"""Validate our CEN loci against Helsen et al. 2026 (Supplementary Data 4) and PCAn.

Conventions
-----------
Our locus: CDEI start .. CDEIII end (SGD: CDEIII = last 26 bp), genomic 1-based
inclusive `cen_start`..`cen_end` on `sseqid`, `strand` = orientation of the CEN.
Helsen SD4: one sequence per chromosome (column `<n>_A`, second variant in `<n>_B`),
oriented like the CEN. Their sequences run one base past the SGD CDEIII end
(checked on S288C: all 16 equal SGD + 1 bp), so `HELSEN_EXTRA_3P = 1`.

Shifts are reported in CEN orientation: `d5` > 0 means the other call starts
downstream of ours (shorter at CDEI), `d3` > 0 means it ends downstream of ours
(longer at CDEIII), both after removing the known convention offset.
"""
import gzip
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner

HELSEN_EXTRA_3P = 1
_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ---------------------------------------------------------------- Helsen SD4

def load_helsen_sd4(path: str | Path) -> pd.DataFrame:
    """Long table: assembly, centromere (CEN1..CEN16), variant (A/B), sequence."""
    raw = pd.read_excel(path, header=1)
    raw = raw.loc[:, ~raw.columns.astype(str).str.startswith("Unnamed")]
    long = raw.melt(id_vars="Assembly", var_name="col", value_name="sequence").dropna()
    long["centromere"] = "CEN" + long["col"].str.split("_").str[0]
    long["variant"] = long["col"].str.split("_").str[1]
    long["sequence"] = long["sequence"].astype(str).str.strip().str.upper()
    long = long.rename(columns={"Assembly": "assembly"})
    return long[["assembly", "centromere", "variant", "sequence"]].reset_index(drop=True)


def _accession_root(acc: str) -> str:
    return str(acc).split(".")[0]


def locate(seq: str, contig: str) -> list[tuple[int, int, str]]:
    """Exact occurrences of seq on either strand: (start, end, strand), 1-based inclusive."""
    hits = []
    for strand, query in (("+", seq), ("-", revcomp(seq))):
        for m in re.finditer(f"(?={re.escape(query)})", contig):
            hits.append((m.start() + 1, m.start() + len(seq), strand))
    return hits


def interval_shifts(ours: dict, start: int, end: int, extra_3p: int = 0) -> tuple[int, int]:
    """(d5, d3) of another call vs our locus, in CEN orientation (see module doc)."""
    if ours["strand"] == "+":
        d5 = start - ours["cen_start"]
        d3 = (end - extra_3p) - ours["cen_end"]
    else:
        d5 = ours["cen_end"] - end
        d3 = ours["cen_start"] - (start + extra_3p)
    return int(d5), int(d3)


def compare_same_assembly(loci: pd.DataFrame, helsen: pd.DataFrame, genomes: dict) -> pd.DataFrame:
    """Our loci vs Helsen calls on the very same assembly (matched by GCA accession).

    `genomes`: sample -> {contig id -> sequence string}.
    """
    helsen = helsen.assign(acc=helsen["assembly"].map(_accession_root))
    rows = []
    for _, our in loci.iterrows():
        our = our.to_dict()
        calls = helsen[(helsen["acc"] == _accession_root(our["sample"]))
                       & (helsen["centromere"] == our["centromere"])]
        base = {"sample": our["sample"], "strain": our.get("strain"),
                "centromere": our["centromere"], "source": "helsen_same_assembly",
                "our_length": len(our["sequence"]), "our_cdeii": our["cdeii_length"]}
        if calls.empty:
            rows.append({**base, "status": "absent_in_helsen"})
            continue
        contig = genomes[our["sample"]][our["sseqid"]]
        for _, call in calls.iterrows():
            row = {**base, "variant": call["variant"],
                   "helsen_length": len(call["sequence"]) - HELSEN_EXTRA_3P}
            hits = locate(call["sequence"], contig)
            if not hits:
                row["status"] = "not_on_our_contig"
            else:
                start, end, strand = hits[0]
                d5, d3 = interval_shifts(our, start, end, HELSEN_EXTRA_3P)
                row.update(helsen_start=start, helsen_end=end, d5=d5, d3=d3,
                           status="identical" if d5 == d3 == 0 else "boundary_shift")
                if strand != our["strand"]:
                    row["status"] = "opposite_strand"
            row["cdeii_diff"] = row["helsen_length"] - row["our_length"]
            rows.append(row)
    return pd.DataFrame(rows)


def _aligner() -> PairwiseAligner:
    al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                         open_gap_score=-5, extend_gap_score=-2)
    # Free end gaps: boundary choices should not count as indels
    al.end_gap_score = 0
    return al


def seq_differences(a: str, b: str) -> dict:
    """SNVs and internal indel bp between two CEN sequences (end gaps ignored)."""
    aln = _aligner().align(a, b)[0]
    ta, tb = aln[0], aln[1]
    first = max(next(i for i, c in enumerate(ta) if c != "-"), next(i for i, c in enumerate(tb) if c != "-"))
    last = min(max(i for i, c in enumerate(ta) if c != "-"), max(i for i, c in enumerate(tb) if c != "-"))
    snv = indel = 0
    for x, y in zip(ta[first:last + 1], tb[first:last + 1]):
        if x == "-" or y == "-":
            indel += 1
        elif x != y:
            snv += 1
    return {"snv": snv, "indel_bp": indel}


def compare_other_assembly(loci: pd.DataFrame, helsen: pd.DataFrame, key: str = "strain") -> pd.DataFrame:
    """Our loci vs Helsen calls on another assembly of the same strain (e.g. Peter 2018).

    Sequence-only: the genomes differ, so coordinates cannot be compared. Of the
    A/B calls the closer one is kept.
    """
    rows = []
    for _, our in loci.iterrows():
        calls = helsen[(helsen["assembly"] == our[key]) & (helsen["centromere"] == our["centromere"])]
        base = {"sample": our["sample"], "strain": our[key], "centromere": our["centromere"],
                "source": "helsen_other_assembly",
                "our_length": len(our["sequence"]), "our_cdeii": our["cdeii_length"]}
        if calls.empty:
            rows.append({**base, "status": "absent_in_helsen"})
            continue
        best = None
        for _, call in calls.iterrows():
            seq = call["sequence"][:-HELSEN_EXTRA_3P]
            diff = seq_differences(our["sequence"], seq)
            cand = {**base, "variant": call["variant"], "helsen_length": len(seq),
                    "n_variants": len(calls), **diff}
            if best is None or (diff["snv"] + diff["indel_bp"]) < (best["snv"] + best["indel_bp"]):
                best = cand
        best["cdeii_diff"] = best["helsen_length"] - best["our_length"]
        best["status"] = ("identical" if best["snv"] == best["indel_bp"] == 0 and best["cdeii_diff"] == 0
                          else "differs")
        rows.append(best)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------- PCAn

def read_pcan(path: str | Path) -> pd.DataFrame:
    """CENsequences.txt table (the part above the '****' footer)."""
    lines = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("*"):
                break
            lines.append(line.rstrip("\n").split("\t"))
    df = pd.DataFrame(lines[1:], columns=lines[0])
    for col in ("Start", "End", "CDEIIlen"):
        df[col] = df[col].astype(int)
    for col in ("CDEIIAT", "FIMOIscore", "FIMOIIscore", "OverallScore"):
        df[col] = df[col].astype(float)
    return df


def compare_pcan(loci: pd.DataFrame, pcan: dict[str, pd.DataFrame], genomes: dict) -> pd.DataFrame:
    """Match PCAn calls to our loci by contig overlap; unmatched calls on either side are kept.

    PCAn's 3' end is the end of its CDEIII motif hit, not the SGD CDEIII end, so `d3`
    carries a constant convention offset that the summary reports as the mode.
    """
    rows = []
    for sample, ours in loci.groupby("sample"):
        calls = pcan.get(sample)
        if calls is None:
            continue
        used = set()
        for _, our in ours.iterrows():
            our = our.to_dict()
            base = {"sample": sample, "strain": our.get("strain"), "centromere": our["centromere"],
                    "our_length": len(our["sequence"]), "our_cdeii": our["cdeii_length"]}
            hit = calls[(calls["Contig"] == our["sseqid"])
                        & (calls["Start"] <= our["cen_end"]) & (calls["End"] >= our["cen_start"])]
            if hit.empty:
                rows.append({**base, "status": "missed_by_pcan"})
                continue
            call = hit.iloc[0]
            used.add(call.name)
            contig = genomes[sample][our["sseqid"]]
            plus = contig[call["Start"] - 1:call["End"]] == call["Sequence"]
            d5, d3 = interval_shifts(our, call["Start"], call["End"])
            rows.append({**base, "status": "found", "pcan_start": call["Start"], "pcan_end": call["End"],
                         "pcan_strand": "+" if plus else "-", "d5": d5, "d3": d3,
                         "pcan_length": len(call["Sequence"]), "pcan_cdeii": call["CDEIIlen"],
                         "pcan_score": call["OverallScore"]})
        for idx, call in calls.drop(index=list(used)).iterrows():
            rows.append({"sample": sample, "status": "pcan_only", "pcan_contig": call["Contig"],
                         "pcan_start": call["Start"], "pcan_end": call["End"],
                         "pcan_length": len(call["Sequence"]), "pcan_score": call["OverallScore"]})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------- genomes

def read_genome(genome_dir: str | Path, sample: str) -> dict[str, str]:
    path = sorted(Path(genome_dir).glob(f"{sample}.*"))[0]
    handle = gzip.open(path, "rt") if path.suffix == ".gz" else open(path)
    with handle:
        return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(handle, "fasta")}
