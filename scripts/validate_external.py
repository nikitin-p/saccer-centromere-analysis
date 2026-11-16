#!/usr/bin/env python3
"""Analysis A: our CEN loci vs Helsen et al. 2026 (SD4) and PCAn.

Inputs (config `validation:`): Helsen SD4 xlsx, PCAn output dir from
scripts/run_pcan.sh. Outputs in config `output.dir`:
  validation_helsen.tsv  per locus x source (same GCA assembly / same strain in Peter 2018)
  validation_pcan.tsv    per locus, plus PCAn-only calls
"""
import argparse
from pathlib import Path

import pandas as pd

from centromere.io import load_config, save_tsv
from centromere.validate import (compare_other_assembly, compare_pcan, compare_same_assembly,
                                 load_helsen_sd4, read_genome, read_pcan)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    out = Path(config["output"]["dir"])
    val = config["validation"]
    loci = pd.read_csv(out / "preprocessed.tsv", sep="\t")
    genomes = {s: read_genome(config["cde"]["genome_dir"], s) for s in loci["sample"].unique()}

    helsen = load_helsen_sd4(val["helsen_sd4"])
    same = compare_same_assembly(loci, helsen, genomes)
    other = compare_other_assembly(loci, helsen)
    table = pd.concat([same, other], ignore_index=True)
    save_tsv(table, out / "validation_helsen.tsv")
    print("Helsen SD4:")
    print(table.groupby(["source", "status"]).size().to_string())
    variants = same[same["variant"] == "B"]
    print(f"  second (B) variants on our assemblies: {len(variants)}")

    pcan_dir = Path(val["pcan_dir"])
    pcan = {p.parent.name: read_pcan(p) for p in pcan_dir.glob("*/CENsequences.txt")}
    pc = compare_pcan(loci, pcan, genomes)
    save_tsv(pc, out / "validation_pcan.tsv")
    print(f"PCAn ({len(pcan)} assemblies):")
    print(pc.groupby("status").size().to_string())
    found = pc[pc["status"] == "found"]
    if len(found):
        print("  d5 (CDEI start shift):", found["d5"].value_counts().sort_index().to_dict())
        print("  d3 (3' end shift):    ", found["d3"].value_counts().sort_index().to_dict())


if __name__ == "__main__":
    main()
