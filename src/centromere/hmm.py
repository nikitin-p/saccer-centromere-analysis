"""L4 locator (locator spec phase 3): one nucleotide HMM per centromere, grown round by round
from the strains it finds.

Round k: build one profile per CEN (MAFFT -> `hmmbuild --dna`) from the sequences accepted so
far, `nhmmer` all 16 against every assembly, and accept a new sequence for CEN i only if
  - CEN i's profile is the best one at that locus, by >= `margin` bits over any other
    profile overlapping it (reciprocal check against cross-centromere capture);
  - CDEI and CDEIII, placed with the S288C PWMs inside the hit, pass `cde_ok`;
  - its length is within S288C +- `len_tol` bp (a longer locus that L2/L3 also anchor is
    reported as `insertion` but never enters the profile);
  - the identical sequence is not already in the profile.

Coordinates in the tables here are contig coordinates, 1-based, inclusive (as L1/L2).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd

from centromere.cde import CDEModel, place_cde

TBL_COLS = ["contig", "t_acc", "profile", "q_acc", "hmm_from", "hmm_to", "ali_from", "ali_to",
            "env_from", "env_to", "sq_len", "strand", "evalue", "bits", "bias"]
REASONS = ["ok", "runner_up", "cde", "length", "insertion", "no_hit"]


# ---------------------------------------------------------------- HMMER / MAFFT wrappers

def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def build_profile(name: str, seqs: list[str], out_dir: Path, hmmbuild: str = "hmmbuild",
                  mafft_args: str = "--auto") -> Path:
    """MAFFT-align `seqs` (one sequence is used as is) and build `<out_dir>/<name>.hmm`."""
    from centromere.align import run_mafft
    aln = run_mafft(seqs, mafft_args) if len(seqs) > 1 else [seqs[0]]
    afa = out_dir / f"{name}.afa"
    write_fasta(afa, [(f"{name}_{i}", s.upper()) for i, s in enumerate(aln)])
    hmm = out_dir / f"{name}.hmm"
    subprocess.run([hmmbuild, "--dna", "--informat", "afa", "-n", name, str(hmm), str(afa)],
                   check=True, capture_output=True)
    return hmm


def run_nhmmer(hmm: Path, target: Path, tbl: Path, nhmmer: str = "nhmmer", cpu: int = 4,
               report_evalue: float = 1.0) -> pd.DataFrame:
    """All profiles in `hmm` against an uncompressed FASTA (nhmmer must rewind it per query)."""
    subprocess.run([nhmmer, "--dna", "--noali", "-o", "/dev/null", "--tblout", str(tbl),
                    "-E", str(report_evalue), "--cpu", str(cpu), str(hmm), str(target)],
                   check=True, capture_output=True)
    return parse_tblout(tbl)


def parse_tblout(path) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.split(maxsplit=len(TBL_COLS))[:len(TBL_COLS)])
    df = pd.DataFrame(rows, columns=TBL_COLS)
    for c in ["hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "sq_len"]:
        df[c] = df[c].astype(int)
    for c in ["evalue", "bits", "bias"]:
        df[c] = df[c].astype(float)
    return df.drop(columns=["t_acc", "q_acc", "bias"])


def model_lengths(hmm_path) -> dict[str, int]:
    """NAME -> LENG (match states) for every model in an HMMER3 file."""
    out, name = {}, None
    with open(hmm_path) as fh:
        for line in fh:
            if line.startswith("NAME "):
                name = line.split()[1]
            elif line.startswith("LENG ") and name is not None:
                out[name] = int(line.split()[1])
    return out


# ---------------------------------------------------------------- hits -> loci

def add_margin(hits: pd.DataFrame) -> pd.DataFrame:
    """Per hit: the best other profile overlapping its envelope (either strand) and the margin
    in bits. No overlapping other profile -> runner-up bits 0."""
    h = hits.copy()
    h["lo"] = h[["env_from", "env_to"]].min(axis=1)
    h["hi"] = h[["env_from", "env_to"]].max(axis=1)
    runner, runner_bits = [], []
    for contig, g in h.groupby("contig", sort=False):
        lo, hi, prof, bits = (g[c].to_numpy() for c in ("lo", "hi", "profile", "bits"))
        for i in range(len(g)):
            m = (lo <= hi[i]) & (hi >= lo[i]) & (prof != prof[i])
            if m.any():
                j = bits[m].argmax()
                runner.append((g.index[i], prof[m][j], bits[m][j]))
            else:
                runner.append((g.index[i], None, 0.0))
    r = pd.DataFrame(runner, columns=["idx", "runner_up", "runner_up_bits"]).set_index("idx")
    h = h.join(r)
    h["margin"] = h["bits"] - h["runner_up_bits"]
    return h


def best_hits(hits: pd.DataFrame, profiles: list[str], evalue: float = 1e-10,
              margin: float = 10) -> pd.DataFrame:
    """One row per profile: its best hit where it beats every other profile; if it never does,
    its top hit (margin <= 0). `n_loci` = hits of the profile with E <= evalue and margin >=
    `margin` (more than 1 = an extra CEN-like locus)."""
    h = add_margin(hits) if len(hits) else hits.assign(runner_up=None, runner_up_bits=0.0, margin=0.0)
    h = h[h["evalue"] <= evalue]
    rows = []
    for p in profiles:
        g = h[h["profile"] == p]
        if g.empty:
            rows.append({"profile": p, "n_loci": 0})
            continue
        win = g[g["margin"] > 0]
        best = (win if len(win) else g).sort_values("bits", ascending=False).iloc[0]
        rows.append({**best.drop(["lo", "hi"]).to_dict(), "n_loci": int((g["margin"] >= margin).sum())})
    return pd.DataFrame(rows)


def project_cen(hit: pd.Series, model_len: int) -> tuple[int, int, bool]:
    """Projected CEN span (1-based, inclusive) if the hit covered the whole model, and strand."""
    plus = hit["strand"] == "+"
    head, tail = int(hit["hmm_from"]) - 1, model_len - int(hit["hmm_to"])
    if plus:
        return int(hit["env_from"]) - head, int(hit["env_to"]) + tail, True
    return int(hit["env_to"]) - tail, int(hit["env_from"]) + head, False


def overlaps(contig: str, start: int, end: int, other: pd.DataFrame | None) -> bool:
    """`other` has contig/start/end columns for the same sample x CEN."""
    if other is None or other.empty:
        return False
    return bool(((other["contig"] == contig) & (other["start"] <= end) & (other["end"] >= start)).any())


def classify(hit: dict, ref_len: int, margin: float = 10, len_tol: int = 40,
             anchored: bool = False) -> str:
    """Status of a locus placed from a hit (`no_hit` is set by the caller)."""
    if hit["margin"] < margin:
        return "runner_up"
    if not hit["cde_ok"]:
        return "cde"
    if abs(hit["length"] - ref_len) > len_tol:
        return "insertion" if anchored else "length"
    return "ok"


def locate_in_genome(hits: pd.DataFrame, contigs: dict[str, str], lengths: dict[str, int],
                     model: CDEModel, ref_len: dict[str, int], anchors: pd.DataFrame | None = None,
                     evalue: float = 1e-10, margin: float = 10, len_tol: int = 40,
                     flank: int = 150, radius: int = 20) -> pd.DataFrame:
    """One row per CEN for one assembly: best hit, CDE-placed locus, status.

    `anchors`: L2/L3 loci of this assembly (centromere, contig, start, end), for `insertion`.
    """
    rows = []
    for _, b in best_hits(hits, list(lengths), evalue, margin).iterrows():
        cen = b["profile"]
        row = {"centromere": cen, "n_loci": b["n_loci"]}
        if pd.isna(b.get("bits")):
            rows.append({**row, "status": "no_hit"})
            continue
        lo, hi, plus = project_cen(b, lengths[cen])
        cde = place_cde(contigs[b["contig"]], lo, hi, plus, model, flank, radius)
        row.update(contig=b["contig"], hmm_start=int(cde["cen_start"]), hmm_end=int(cde["cen_end"]),
                   strand="+" if plus else "-", length=len(cde["sequence"]),
                   bits=b["bits"], evalue=b["evalue"], hmm_from=int(b["hmm_from"]),
                   hmm_to=int(b["hmm_to"]), model_len=lengths[cen], runner_up=b["runner_up"],
                   runner_up_bits=b["runner_up_bits"], margin=round(float(b["margin"]), 1),
                   cdei_score=cde["cdei_score"], cdeiii_score=cde["cdeiii_score"],
                   cde_ok=cde["cde_ok"], sequence=cde["sequence"])
        anc = None if anchors is None else anchors[anchors["centromere"] == cen]
        row["anchored"] = overlaps(row["contig"], row["hmm_start"], row["hmm_end"], anc)
        row["status"] = classify(row, ref_len[cen], margin, len_tol, row["anchored"])
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- iteration bookkeeping

def new_members(loci: pd.DataFrame, profile: dict[str, set[str]], frozen: set[str]) -> pd.DataFrame:
    """Loci with status ok whose sequence is not yet in their CEN's profile (first occurrence
    kept when several strains share a new sequence)."""
    if loci.empty:
        return loci
    ok = loci[(loci["status"] == "ok") & ~loci["centromere"].isin(frozen)]
    ok = ok[[s not in profile.get(c, set()) for c, s in zip(ok["centromere"], ok["sequence"])]]
    return ok.drop_duplicates(["centromere", "sequence"])


def drifting(rounds: pd.DataFrame, k: int) -> set[str]:
    """CENs whose `runner_up` rejections rose from round k-1 to round k (iteration stops)."""
    if k == 0:
        return set()
    a = rounds[rounds["round"] == k].set_index("centromere")["rejected_runner_up"]
    b = rounds[rounds["round"] == k - 1].set_index("centromere")["rejected_runner_up"]
    d = a.sub(b, fill_value=0)
    return set(d[d > 0].index)


def round_summary(k: int, loci: pd.DataFrame, profile: dict[str, set[str]], accepted: pd.DataFrame,
                  s288c_bits: dict[str, float], cens: list[str]) -> pd.DataFrame:
    rows = []
    for cen in cens:
        g = loci[loci["centromere"] == cen]
        st = g["status"].value_counts()
        rows.append({"round": k, "centromere": cen, "n_profile": len(profile[cen]),
                     "found_ok": int(st.get("ok", 0)),
                     "accepted": int((accepted["centromere"] == cen).sum()) if len(accepted) else 0,
                     **{f"rejected_{r}": int(st.get(r, 0)) for r in REASONS[1:]},
                     "s288c_bits": s288c_bits.get(cen)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- phase 4: CEN-like scan

def scan_cen_like(hits: pd.DataFrame, contigs: dict[str, str], lengths: dict[str, int],
                  model: CDEModel, known: pd.DataFrame | None = None, evalue: float = 1e-5,
                  margin: float = 10, flank: int = 150, radius: int = 20) -> pd.DataFrame:
    """Profile hits that look like a centromere but lie outside every known locus.

    Kept: E <= `evalue`, best profile at the spot by >= `margin` bits, and no overlap with
    `known` (contig, start, end: the loci of any CEN in this assembly). One row per spot (the
    best profile's hit), CDEI/CDEIII placed inside the projected span; `cde_ok` is the CDE
    gate, applied by the caller (a CEN-like element needs it). `hmm_cov` = fraction of the
    model the hit covers (a CEN fragment covers part of it)."""
    if hits.empty:
        return pd.DataFrame()
    h = add_margin(hits)
    h = h[(h["evalue"] <= evalue) & (h["margin"] >= margin)]
    rows = []
    for _, b in h.sort_values("bits", ascending=False).iterrows():
        lo, hi, plus = project_cen(b, lengths[b["profile"]])
        lo, hi = max(1, lo), min(len(contigs[b["contig"]]), hi)
        if known is not None and overlaps(b["contig"], lo, hi, known):
            continue
        if any(r["contig"] == b["contig"] and r["start"] <= hi and r["end"] >= lo for r in rows):
            continue
        cde = place_cde(contigs[b["contig"]], lo, hi, plus, model, flank, radius)
        rows.append({"profile": b["profile"], "contig": b["contig"], "start": int(cde["cen_start"]),
                     "end": int(cde["cen_end"]), "strand": "+" if plus else "-",
                     "length": len(cde["sequence"]), "bits": b["bits"], "evalue": b["evalue"],
                     "margin": round(float(b["margin"]), 1), "runner_up": b["runner_up"],
                     "hmm_from": int(b["hmm_from"]), "hmm_to": int(b["hmm_to"]),
                     "hmm_cov": round((b["hmm_to"] - b["hmm_from"] + 1) / lengths[b["profile"]], 2),
                     "cdei_score": cde["cdei_score"], "cdeiii_score": cde["cdeiii_score"],
                     "cde_ok": cde["cde_ok"], "sequence": cde["sequence"]})
    return pd.DataFrame(rows)
