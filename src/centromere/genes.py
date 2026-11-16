"""Analysis C: kinetochore genes extracted from assemblies with tblastn.

Per genome and S288C protein: pick the best tblastn locus (HSP cluster on one contig
and strand), find the in-frame ATG nearest the projected protein start, read to the
first in-frame stop, translate, and compare with the S288C protein.
"""
import gzip
import subprocess
from pathlib import Path

import pandas as pd
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

from centromere.io import tool

HSP_COLS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send",
            "evalue", "bitscore", "sframe", "qlen", "slen"]
STOPS = {"TAA", "TAG", "TGA"}
START_WINDOW = 30  # codons searched around the projected start for an ATG
STOP_TOLERANCE = 10  # codons; a stop earlier than this before the projected end is premature


def gene_name(qseqid: str) -> str:
    """`CBF1|P17106` -> `CBF1`."""
    return qseqid.split("|")[0]


# ---------------------------------------------------------------- BLAST
def make_db(genome: Path, db_dir: Path, blast_bin=None) -> Path:
    """Uncompressed copy of the genome + nucleotide BLAST DB; returns the FASTA path (= DB name)."""
    db_dir.mkdir(parents=True, exist_ok=True)
    name = genome.name.removesuffix(".gz").removesuffix(".fasta").removesuffix(".fa")
    fa = db_dir / f"{name}.fa"
    if not fa.exists():
        opener = gzip.open if genome.suffix == ".gz" else open
        with opener(genome, "rt") as src, open(fa, "w") as dst:
            for rec in SeqIO.parse(src, "fasta"):
                dst.write(f">{rec.id}\n{str(rec.seq).upper()}\n")
    if not Path(f"{fa}.nsq").exists() and not Path(f"{fa}.nal").exists():
        subprocess.run([tool(blast_bin, "makeblastdb"), "-in", str(fa), "-dbtype", "nucl"],
                       check=True, capture_output=True)
    return fa


def run_tblastn(proteins: Path, db: Path, blast_bin, evalue: float) -> pd.DataFrame:
    # -seg no: SEG would mask low-complexity stretches (e.g. Cbf1 N-terminus) and split HSPs
    cmd = [tool(blast_bin, "tblastn"), "-query", str(proteins), "-db", str(db),
           "-evalue", str(evalue), "-seg", "no", "-comp_based_stats", "0",
           "-max_target_seqs", "50", "-outfmt", "6 " + " ".join(HSP_COLS)]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    rows = [line.split("\t") for line in res.stdout.splitlines() if line]
    df = pd.DataFrame(rows, columns=HSP_COLS)
    for c in HSP_COLS[2:]:
        df[c] = pd.to_numeric(df[c])
    return df


# ---------------------------------------------------------------- locus
def oriented(hsps: pd.DataFrame) -> pd.DataFrame:
    """Add strand and 0-based half-open coordinates on the strand the gene is read from."""
    h = hsps.copy()
    plus = h["sstart"] <= h["send"]
    h["strand"] = plus.map({True: "+", False: "-"})
    h["o_start"] = (h["sstart"] - 1).where(plus, h["slen"] - h["sstart"])
    h["o_end"] = h["send"].where(plus, h["slen"] - h["send"] + 1)
    return h


def pick_locus(hsps: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """HSPs of one protein in one genome -> (HSPs of the best locus, number of loci scoring
    >= 50 % of the best, i.e. possible extra copies)."""
    h = oriented(hsps).sort_values(["sseqid", "strand", "o_start"])
    qlen = int(h["qlen"].iloc[0])
    cluster, cid, cur_key, cur_end = [], -1, None, 0
    for _, r in h.iterrows():
        key = (r["sseqid"], r["strand"])
        if key != cur_key or r["o_start"] - cur_end > 3 * qlen:
            cid, cur_key, cur_end = cid + 1, key, r["o_end"]
        else:
            cur_end = max(cur_end, r["o_end"])
        cluster.append(cid)
    h["cluster"] = cluster
    score = h.groupby("cluster")["bitscore"].sum()
    best = score.idxmax()
    return h[h["cluster"] == best].sort_values("qstart"), int((score >= 0.5 * score.max()).sum())


# ---------------------------------------------------------------- ORF
def _codon(seq: str, i: int) -> str:
    return seq[i:i + 3] if 0 <= i <= len(seq) - 3 else ""


def extract_cds(seq: str, locus: pd.DataFrame) -> dict:
    """`seq` = contig on the gene's strand; `locus` = oriented HSPs of one locus.

    Returns start/end (0-based, oriented), cds, protein (no stop), frameshift and status flags.
    """
    first = locus.loc[locus["qstart"].idxmin()]
    last = locus.loc[locus["qend"].idxmax()]
    qlen = int(first["qlen"])
    anchor = int(first["o_start"])
    projected = anchor - (int(first["qstart"]) - 1) * 3
    frames = {int(s) % 3 for s in locus["o_start"]}

    # in-frame ATG nearest the projected start, not separated from the HSP by a stop
    start, best_d = None, None
    for k in range(-START_WINDOW, START_WINDOW + 1):
        i = projected + 3 * k
        # downstream ATGs (inside the HSP) allowed: N-terminally shorter protein
        if i < 0 or _codon(seq, i) != "ATG":
            continue
        if any(_codon(seq, j) in STOPS for j in range(i, anchor, 3)):
            continue
        d = (abs(k), k > 0)  # tie -> upstream
        if best_d is None or d < best_d:
            start, best_d = i, d
    start_found = start is not None
    if not start_found:
        start = anchor  # read from the first aligned codon

    end = start
    while end <= len(seq) - 3 and _codon(seq, end) not in STOPS:
        end += 3
    stop_found = end <= len(seq) - 3
    if stop_found:
        end += 3
    cds = seq[start:end]
    protein = str(Seq(cds[: len(cds) - len(cds) % 3]).translate()).rstrip("*")

    projected_end = int(last["o_end"]) + (qlen - int(last["qend"])) * 3 + 3
    if len(frames) > 1:
        status = "frameshift"
    elif stop_found and end < projected_end - 3 * STOP_TOLERANCE:
        status = "premature_stop"
    elif not start_found or not stop_found:
        status = "partial"
    else:
        status = "intact"
    return {"o_start": start, "o_end": end, "cds": cds, "protein": protein,
            "start_shift": (start - projected) // 3, "frames": len(frames), "start_found": start_found, "stop_found": stop_found,
            "projected_end": projected_end, "status": status}


def to_plus(o_start: int, o_end: int, strand: str, slen: int) -> tuple[int, int]:
    """Oriented 0-based half-open -> 1-based inclusive plus-strand coordinates."""
    if strand == "+":
        return o_start + 1, o_end
    return slen - o_end + 1, slen - o_start


# ---------------------------------------------------------------- variants
def _aligner() -> Align.PairwiseAligner:
    a = Align.PairwiseAligner(mode="global")
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -10, -0.5
    # free end gaps: a protein starting at an internal Met aligns as an N-terminal truncation,
    # not as Met1 + internal deletion
    a.open_end_gap_score = a.extend_end_gap_score = 0
    return a


def protein_variants(ref: str, alt: str) -> list[dict]:
    """Substitutions and indels of `alt` vs `ref` in 1-based reference coordinates.

    Gaps running to the end of either sequence are `truncation` / `extension`.
    """
    aln = _aligner().align(ref, alt)[0]
    r, a = aln[0], aln[1]
    out, pos, i, n = [], 0, 0, len(r)
    while i < n:
        if r[i] != "-" and a[i] != "-":
            pos += 1
            if r[i] != a[i]:
                out.append({"ref_pos": pos, "ref": r[i], "alt": a[i], "type": "substitution"})
            i += 1
            continue
        j = i
        gap_in_alt = a[i] == "-"
        while j < n and (a[j] == "-" if gap_in_alt else r[j] == "-"):
            j += 1
        end = "n_" if i == 0 else "c_" if j == n else ""
        if gap_in_alt:
            seg = r[i:j]
            kind = f"{end}truncation" if end else "deletion"
            out.append({"ref_pos": pos + 1, "ref": seg, "alt": "-", "type": kind})
            pos += len(seg)
        else:
            kind = f"{end}extension" if end else "insertion"
            out.append({"ref_pos": pos, "ref": "-", "alt": a[i:j], "type": kind})
        i = j
    return out


def variant_label(v: dict) -> str:
    """Short label: `S12T`, `del45-47`, `ins12_A`, `Ntrunc1-11`, `Ctrunc300-478`, `Cext_5aa`."""
    t = v["type"]
    if t == "substitution":
        return f"{v['ref']}{v['ref_pos']}{v['alt']}"
    span = f"{v['ref_pos']}-{v['ref_pos'] + len(v['ref']) - 1}"
    if t == "deletion":
        return f"del{span}"
    if t.endswith("truncation"):
        return f"{t[0].upper()}trunc{span}"
    if t == "insertion":
        return f"ins{v['ref_pos']}_{v['alt']}"
    return f"{t[0].upper()}ext_{len(v['alt'])}aa"


def longest_run(s: str) -> int:
    best = cur = 0
    for k, ch in enumerate(s):
        cur = cur + 1 if k and ch == s[k - 1] else 1
        best = max(best, cur)
    return best


def frameshift_sites(seq: str, locus: pd.DataFrame, flank: int = 10) -> list[dict]:
    """Frame changes between consecutive HSPs (by query order): approximate query position,
    oriented genome position, and the longest homopolymer within +-`flank` bp — a long run
    points to a possible assembly error rather than a real indel."""
    h = locus.sort_values("qstart").reset_index(drop=True)
    out = []
    for k in range(1, len(h)):
        a, b = h.iloc[k - 1], h.iloc[k]
        if int(a["o_start"]) % 3 == int(b["o_start"]) % 3:
            continue
        g = (int(a["o_end"]) + int(b["o_start"])) // 2
        out.append({"query_pos": (int(a["qend"]) + int(b["qstart"])) // 2, "o_pos": g,
                    "homopolymer": longest_run(seq[max(0, g - flank):g + flank])})
    return out


def assign_haplotypes(genes: pd.DataFrame, reference: dict[str, str]) -> pd.Series:
    """Per gene: `ref` for the S288C protein, else H1, H2, ... by descending frequency."""
    labels = pd.Series(index=genes.index, dtype=object)
    for gene, g in genes.groupby("gene"):
        counts = g.loc[g["protein"] != reference[gene], "protein"].value_counts(sort=False)
        order = sorted(counts.index, key=lambda p: (-counts[p], g.index[g["protein"] == p][0]))
        names = {p: f"H{k}" for k, p in enumerate(order, 1)}
        names[reference[gene]] = "ref"
        labels[g.index] = g["protein"].map(names)
    return labels


# ---------------------------------------------------------------- C.3 domains
def msa_positions(aln: dict[str, str], ref_name: str) -> pd.DataFrame:
    """Variable columns of a protein MSA in reference coordinates.

    `aln` = {name: aligned sequence}, must contain `ref_name`. Gaps count as an allele (`-`).
    Columns where the reference has a gap are insertions, placed after `ref_pos`.
    Columns: col, ref_pos, ref, insertion, alleles (`R:9,C:1`), n_alt (sequences differing
    from the reference), n_sub (of them, carrying a residue rather than a gap), n,
    alt_strains, ref_private (every non-reference sequence differs from the reference).
    """
    ref = aln[ref_name]
    others = {k: v for k, v in aln.items() if k != ref_name}
    rows, pos = [], 0
    for c, r in enumerate(ref):
        ins = r == "-"
        if not ins:
            pos += 1
        col = {k: v[c] for k, v in others.items()}
        alt = sorted(k for k, a in col.items() if a != r)
        if not alt:
            continue
        counts = pd.Series(list(col.values())).value_counts()
        rows.append({"col": c + 1, "ref_pos": pos, "ref": r, "insertion": ins,
                     "alleles": ",".join(f"{a}:{n}" for a, n in counts.items()),
                     "n_alt": len(alt), "n_sub": sum(col[k] != "-" for k in alt),
                     "n": len(col), "alt_strains": ",".join(alt),
                     "ref_private": len(alt) == len(col)})
    return pd.DataFrame(rows, columns=["col", "ref_pos", "ref", "insertion", "alleles", "n_alt",
                                       "n_sub", "n", "alt_strains", "ref_private"])


def domains_at(pos: int, domains: pd.DataFrame) -> str:
    """`;`-joined domains covering reference position `pos` of one gene; `-` if none."""
    hit = domains[(domains["start"] <= pos) & (domains["end"] >= pos)]["domain"]
    return ";".join(dict.fromkeys(hit)) or "-"


def domain_summary(positions: pd.DataFrame, domains: pd.DataFrame, length: int) -> pd.DataFrame:
    """Per domain of one gene (+ `unannotated` = positions in no domain): length, variable
    reference positions (any difference, gaps included), substituted positions, those
    excluding S288C-private alleles, and that last count per 100 aa. The rate uses
    substitutions only: an N-terminal truncation would otherwise count once per residue."""
    pos = positions[~positions["insertion"]]
    covered = set()
    for d in domains.itertuples():
        covered |= set(range(d.start, d.end + 1))
    regions = []
    for name, g in domains.groupby("domain", sort=False):
        regions.append((name, set().union(*(range(s, e + 1) for s, e in zip(g["start"], g["end"])))))
    regions.append(("unannotated", set(range(1, length + 1)) - covered))
    rows = []
    for name, span in regions:
        p = pos[pos["ref_pos"].isin(span)]
        sub = p[p["n_sub"] > 0]
        n_shared = sub.loc[~sub["ref_private"], "ref_pos"].nunique()
        rows.append({"domain": name, "length": len(span), "n_variable": p["ref_pos"].nunique(),
                     "n_substituted": sub["ref_pos"].nunique(),
                     "n_substituted_excl_ref_private": n_shared,
                     "subst_per_100aa": round(100 * n_shared / len(span), 2) if span else None})
    return pd.DataFrame(rows)
