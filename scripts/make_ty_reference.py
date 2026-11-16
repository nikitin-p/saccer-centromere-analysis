#!/usr/bin/env python3
"""Analysis D step 1: Ty/LTR query library, S288C tRNAs and SGD TE coordinates from the SGD GFF.

Writes (committed, small):
  data/ty_reference.fasta   INT|<family>|<name>  internal region of one representative per family
                            LTR|<family>|<name>  every full-length nuclear S288C LTR
  data/s288c_trna.fasta     275 nuclear tRNA genes
  data/s288c_te_sgd.tsv     SGD Ty elements and LTRs (chrom, start, end, strand, type, family, name)
  data/s288c_cen.bed        16 S288C CENs (BED, 0-based start; strand = CDEI->CDEIII orientation)

Asserts the GFF genome equals the local sacCer3 (chrI-XVI) and that the 16 SGD CEN
intervals equal data/sacer_centromere.fasta, so coordinates and sequences agree.
"""
import argparse
import gzip
import io
from pathlib import Path
from urllib.parse import unquote

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
         "XI", "XII", "XIII", "XIV", "XV", "XVI"]
CHROMS = [f"chr{r}" for r in ROMAN]
FULL_LTR_FRACTION = 0.9   # an LTR copy is "full-length" if >= 90 % of its family's modal length


def read_gff(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    with gzip.open(path, "rt") as f:
        text = f.read()
    head, fasta = text.split("##FASTA\n", 1)
    rows = [line.split("\t") for line in head.splitlines() if line and not line.startswith("#")]
    gff = pd.DataFrame([r for r in rows if len(r) == 9],
                       columns=["chrom", "src", "type", "start", "end", "score", "strand", "phase", "attr"])
    gff[["start", "end"]] = gff[["start", "end"]].astype(int)
    attrs = gff["attr"].map(lambda a: dict(kv.split("=", 1) for kv in a.split(";") if "=" in kv))
    for key in ("ID", "Name", "Alias", "Note", "Parent"):
        gff[key] = attrs.map(lambda d: unquote(d.get(key, "")))
    genome = {r.id: str(r.seq).upper() for r in SeqIO.parse(io.StringIO(fasta), "fasta")}
    return gff, genome


def cut(genome: dict[str, str], chrom: str, start: int, end: int, strand: str) -> str:
    s = genome[chrom][start - 1:end]
    return str(Seq(s).reverse_complement()) if strand == "-" else s


def check_reference(gff, genome, sacCer3: Path, cen_fasta: Path) -> None:
    local = {r.id: str(r.seq).upper() for r in SeqIO.parse(sacCer3, "fasta") if r.id in CHROMS}
    bad = [c for c in CHROMS if local.get(c) != genome.get(c)]
    assert not bad, f"SGD GFF genome differs from {sacCer3} on {bad}"
    ref = {r.id: str(r.seq).upper() for r in SeqIO.parse(cen_fasta, "fasta")}
    cen = gff[(gff["type"] == "centromere") & (gff["Parent"] == "")]
    for _, r in cen.iterrows():
        assert cut(genome, r.chrom, r.start, r.end, r.strand) == ref[r.Name], f"{r.Name} mismatch"
    assert len(cen) == 16, len(cen)
    bed = cen.assign(start0=cen["start"] - 1, score=0)[["chrom", "start0", "end", "Name", "score", "strand"]]
    bed.to_csv(cen_fasta.parent / "s288c_cen.bed", sep="\t", header=False, index=False)
    print("reference check: GFF genome == sacCer3 (chrI-XVI); 16/16 CENs == sacer_centromere.fasta")


def te_table(gff: pd.DataFrame) -> pd.DataFrame:
    ty = gff[(gff["type"] == "LTR_retrotransposon") & gff["chrom"].isin(CHROMS)].copy()
    ty["family"] = ty["Alias"]
    ltr = gff[(gff["type"] == "long_terminal_repeat") & gff["chrom"].isin(CHROMS)].copy()
    ltr["family"] = ltr["Note"].str.replace(" LTR", "", regex=False)
    te = pd.concat([ty, ltr])[["chrom", "start", "end", "strand", "type", "family", "Name"]]
    te = te.rename(columns={"Name": "name"}).sort_values(["chrom", "start"]).reset_index(drop=True)
    assert te["family"].isin(["Ty1", "Ty2", "Ty3", "Ty4", "Ty5"]).all(), te["family"].unique()
    return te


def representative(group: pd.DataFrame) -> pd.Series:
    """Element closest to the family's median length (longest when < 3 copies)."""
    length = group["end"] - group["start"] + 1
    if len(group) < 3:
        return group.loc[length.idxmax()]
    return group.loc[(length - length.median()).abs().idxmin()]


def library(te: pd.DataFrame, genome: dict[str, str]) -> list[tuple[str, str]]:
    out = []
    ty = te[te["type"] == "LTR_retrotransposon"]
    ltr = te[te["type"] == "long_terminal_repeat"]
    for fam, group in ty.groupby("family"):
        r = representative(group)
        inside = ltr[(ltr.chrom == r.chrom) & (ltr.start >= r.start) & (ltr.end <= r.end)]
        start = max([r.start] + [e + 1 for e in inside.end if e < r.start + 1000])
        end = min([r.end] + [s - 1 for s in inside.start if s > r.end - 1000])
        out.append((f"INT|{fam}|{r['name']}", cut(genome, r.chrom, start, end, r.strand)))
    ltr_len = ltr["end"] - ltr["start"] + 1
    modal = ltr_len.groupby(ltr["family"]).agg(lambda x: x.round(-1).mode().iloc[0])
    keep = ltr[ltr_len >= FULL_LTR_FRACTION * ltr["family"].map(modal)]
    for _, r in keep.iterrows():
        out.append((f"LTR|{r.family}|{r['name']}", cut(genome, r.chrom, r.start, r.end, r.strand)))
    print("modal LTR length per family:", modal.to_dict())
    print(f"library: {len(ty['family'].unique())} internal regions, "
          f"{len(keep)}/{len(ltr)} full-length LTRs")
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--gff", default="data/external/sgd/saccharomyces_cerevisiae.gff.gz")
    p.add_argument("--sacCer3", default="data/external/s288c/sacCer3.fa")
    p.add_argument("--cen-fasta", default="data/sacer_centromere.fasta")
    p.add_argument("--out-dir", default="data")
    a = p.parse_args()
    out = Path(a.out_dir)
    gff, genome = read_gff(Path(a.gff))
    check_reference(gff, genome, Path(a.sacCer3).expanduser(), Path(a.cen_fasta))

    te = te_table(gff)
    te.to_csv(out / "s288c_te_sgd.tsv", sep="\t", index=False)
    with open(out / "ty_reference.fasta", "w") as f:
        for name, seq in library(te, genome):
            f.write(f">{name}\n{seq}\n")

    trna = gff[(gff["type"] == "tRNA_gene") & gff["chrom"].isin(CHROMS)]
    with open(out / "s288c_trna.fasta", "w") as f:
        for _, r in trna.iterrows():
            f.write(f">{r.Name}\n{cut(genome, r.chrom, r.start, r.end, r.strand)}\n")
    print(f"tRNAs: {len(trna)}; SGD TEs: {te['type'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
