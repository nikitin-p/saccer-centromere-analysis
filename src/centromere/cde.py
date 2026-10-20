"""Per-strain CDEI / CDEII / CDEIII annotation of BLAST-located centromeres.

BLAST only says roughly where a centromere is. Here the locus (hit +/- flank)
is cut from the assembly, oriented like the S288C reference, and CDEI and
CDEIII are placed with position weight matrices built from the 16 S288C
elements. CDEII is whatever lies between them, so insertions are kept.

S288C convention (SGD): CDEI = first 8 bp of the CEN, CDEIII = last 26 bp.
"""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

CDEI_LEN = 8
CDEIII_LEN = 26
BASES = "ACGT"


def read_reference(path: str | Path) -> dict[str, str]:
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(path, "fasta")}


def build_pwm(sites: list[str], pseudocount: float = 0.5) -> np.ndarray:
    """Log2-odds matrix (len x 4) against a uniform background."""
    counts = np.full((len(sites[0]), 4), pseudocount)
    for site in sites:
        for i, base in enumerate(site):
            if base in BASES:
                counts[i, BASES.index(base)] += 1
    freqs = counts / counts.sum(axis=1, keepdims=True)
    return np.log2(freqs / 0.25)


def score_positions(seq: str, pwm: np.ndarray) -> np.ndarray:
    """PWM score of every window start in seq (non-ACGT bases score 0)."""
    width = pwm.shape[0]
    idx = np.array([BASES.find(b) for b in seq])
    n = len(seq) - width + 1
    if n <= 0:
        return np.array([])
    scores = np.zeros(n)
    for i in range(width):
        col = idx[i:i + n]
        scores += np.where(col >= 0, pwm[i, np.clip(col, 0, 3)], 0.0)
    return scores


class CDEModel:
    def __init__(self, reference: dict[str, str]):
        self.reference = reference
        self.cdei = build_pwm([s[:CDEI_LEN] for s in reference.values()])
        self.cdeiii = build_pwm([s[-CDEIII_LEN:] for s in reference.values()])
        # Calls scoring well below the weakest S288C element get flagged
        self.cdei_min = min(score_positions(s[:CDEI_LEN], self.cdei)[0] for s in reference.values()) - 2
        self.cdeiii_min = min(score_positions(s[-CDEIII_LEN:], self.cdeiii)[0] for s in reference.values()) - 2

    def annotate(self, window: str, start_guess: int, end_guess: int, radius: int = 20) -> dict:
        """Place CDEI near start_guess and CDEIII ending near end_guess (window coords)."""
        s1 = score_positions(window, self.cdei)
        lo1, hi1 = max(0, start_guess - radius), min(len(s1), start_guess + radius + 1)
        i1 = lo1 + int(np.argmax(s1[lo1:hi1]))

        s3 = score_positions(window, self.cdeiii)
        lo3 = max(i1 + CDEI_LEN, end_guess - CDEIII_LEN - radius)
        hi3 = min(len(s3), end_guess - CDEIII_LEN + radius + 1)
        if hi3 <= lo3:  # projected end clashes with CDEI: search everything after it
            hi3 = len(s3)
        i3 = lo3 + int(np.argmax(s3[lo3:hi3]))

        return {
            "cdei_start": i1,
            "cdeiii_end": i3 + CDEIII_LEN,
            "cdei_score": round(float(s1[i1]), 2),
            "cdeiii_score": round(float(s3[i3]), 2),
            "cdei_seq": window[i1:i1 + CDEI_LEN],
            "cdeii_seq": window[i1 + CDEI_LEN:i3],
            "cdeiii_seq": window[i3:i3 + CDEIII_LEN],
            "cde_ok": bool(s1[i1] >= self.cdei_min and s3[i3] >= self.cdeiii_min),
        }


def _open_genome(genome_dir: Path, sample: str) -> dict:
    matches = sorted(genome_dir.glob(f"{sample}.*"))
    if not matches:
        raise FileNotFoundError(f"No assembly for {sample} in {genome_dir}")
    path = matches[0]
    handle = gzip.open(path, "rt") if path.suffix == ".gz" else open(path)
    with handle:
        return {rec.id: rec for rec in SeqIO.parse(handle, "fasta")}


def _contig_name(record) -> str:
    # ScRAP headers end in "contig: chrI"; fall back to the id elsewhere
    m = re.search(r"contig:\s*(\S+)", record.description)
    return m.group(1) if m else record.id


def annotate_cde(
    df: pd.DataFrame,
    genome_dir: str | Path,
    reference: dict[str, str],
    flank: int = 150,
    radius: int = 20,
) -> pd.DataFrame:
    """Replace BLAST-trimmed sequences with CDEI-start..CDEIII-end loci."""
    model = CDEModel(reference)
    genome_dir = Path(genome_dir)
    rows = []

    for sample, group in df.groupby("sample", sort=False):
        genome = _open_genome(genome_dir, sample)
        for _, row in group.iterrows():
            row = row.to_dict()
            record = genome[row["sseqid"]]
            qlen = len(reference[row["centromere"]])
            plus = row["send"] >= row["sstart"]
            lo, hi = sorted((row["sstart"], row["send"]))  # 1-based, inclusive

            # Where the full reference CEN would start/end if the hit were extended
            if plus:
                cen_lo, cen_hi = lo - (row["qstart"] - 1), hi + (qlen - row["qend"])
            else:
                cen_lo, cen_hi = lo - (qlen - row["qend"]), hi + (row["qstart"] - 1)

            w_lo = max(1, cen_lo - flank)
            w_hi = min(len(record.seq), cen_hi + flank)
            window = record.seq[w_lo - 1:w_hi]
            if plus:
                start_guess, end_guess = cen_lo - w_lo, cen_hi - w_lo + 1
            else:
                window = window.reverse_complement()
                start_guess, end_guess = w_hi - cen_hi, w_hi - cen_lo + 1
            window = str(window).upper()

            cde = model.annotate(window, start_guess, end_guess, radius)
            if plus:
                g_start, g_end = w_lo + cde["cdei_start"], w_lo + cde["cdeiii_end"] - 1
            else:
                g_start, g_end = w_hi - cde["cdeiii_end"] + 1, w_hi - cde["cdei_start"]

            row["blast_sequence"] = row["sequence"]
            row["sequence"] = window[cde["cdei_start"]:cde["cdeiii_end"]]
            row["sequence_length"] = len(row["sequence"])
            row["strand"] = "+" if plus else "-"
            row["contig_name"] = _contig_name(record)
            row["cen_start"], row["cen_end"] = g_start, g_end
            row["cdeii_length"] = len(cde["cdeii_seq"])
            for key in ("cdei_seq", "cdeii_seq", "cdeiii_seq", "cdei_score", "cdeiii_score", "cde_ok"):
                row[key] = cde[key]
            rows.append(row)

    return pd.DataFrame(rows)
