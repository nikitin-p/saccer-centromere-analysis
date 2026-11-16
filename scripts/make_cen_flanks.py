#!/usr/bin/env python3
"""L2 locator input: 2-kb S288C flanks of the 16 CENs (locator spec, phase 1).

Reads the S288C R64 genome (sacCer3 naming, chrI..chrXVI) and the SGD coordinates in the
`data/sacer_centromere.fasta` headers ("Chromosome VII:496920..497038"). The strand of each CEN
is the one on which the extracted interval equals the FASTA sequence; the script fails if neither
strand matches (reference-version guard).

Writes
  data/s288c_cen_flanks.fasta  32 records `CENn|L` / `CENn|R`, each in CEN orientation
                               (CDEI -> CDEIII): L = flank_bp before CDEI, R = flank_bp after CDEIII
  data/s288c_cen.bed           16 CENs, 0-based start, strand = CDEI->CDEIII orientation
                               (if the file exists it must be identical; it is shared with analysis D)

Usage: scripts/make_cen_flanks.py <sacCer3.fa> [--flank-bp 2000]
"""
import argparse
import re
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
HEADER = re.compile(r"Chromosome ([IVX]+):(\d+)\.\.(\d+)")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genome", type=Path, help="S288C R64 genome FASTA (chrI..chrXVI)")
    ap.add_argument("--reference", type=Path, default=ROOT / "data/sacer_centromere.fasta")
    ap.add_argument("--flank-bp", type=int, default=2000)
    ap.add_argument("--out", type=Path, default=ROOT / "data/s288c_cen_flanks.fasta")
    ap.add_argument("--bed", type=Path, default=ROOT / "data/s288c_cen.bed")
    args = ap.parse_args()

    genome = {r.id: str(r.seq).upper() for r in SeqIO.parse(args.genome.expanduser(), "fasta")}
    flank, bed, records, n_ok = args.flank_bp, [], [], 0
    for rec in SeqIO.parse(args.reference, "fasta"):
        m = HEADER.search(rec.description)
        assert m, f"{rec.id}: no SGD coordinates in header"
        chrom, start, end = "chr" + m.group(1), int(m.group(2)) - 1, int(m.group(3))
        seq, cen = genome[chrom], str(rec.seq).upper()
        core = seq[start:end]
        if core == cen:
            strand = "+"
        elif revcomp(core) == cen:
            strand = "-"
        else:
            raise SystemExit(f"{rec.id}: {chrom}:{start + 1}-{end} does not equal the reference FASTA "
                             "on either strand (genome / reference version mismatch?)")
        assert start >= flank and end + flank <= len(seq), f"{rec.id}: flank runs off {chrom}"
        n_ok += 1
        lo, hi = seq[start - flank:start], seq[end:end + flank]
        left, right = (lo, hi) if strand == "+" else (revcomp(hi), revcomp(lo))
        records += [(f"{rec.id}|L", left), (f"{rec.id}|R", right)]
        bed.append(f"{chrom}\t{start}\t{end}\t{rec.id}\t0\t{strand}\n")
    print(f"{n_ok}/16 extracted CENs equal the reference FASTA")

    with open(args.out, "w") as fh:
        for name, s in records:
            fh.write(f">{name}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    print(f"wrote {args.out} ({len(records)} flanks x {flank} bp)")

    text = "".join(bed)
    if args.bed.exists():
        assert args.bed.read_text() == text, f"{args.bed} differs from the CENs derived here"
        print(f"{args.bed} unchanged (identical)")
    else:
        args.bed.write_text(text)
        print(f"wrote {args.bed}")


if __name__ == "__main__":
    main()
