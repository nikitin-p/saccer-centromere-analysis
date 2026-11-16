#!/usr/bin/env python3
"""Analysis G: aneuploidy as a centromere phenotype on the Peter 2018 748-strain BLAST set.

Runs 01_preprocess + 02_filter on config.peter748.yaml (skip with --skip-pipeline), then:
  aneuploidy_parsed.tsv        strain, clade, ploidy, raw Aneuploidies, parsed events
  aneuploidy_counts.tsv        per chromosome: gains/losses, all Peter 2018 strains vs analysed strains
  cen_calls.tsv                strain x CEN: loci, complete loci, non-S288C, length variant, delta vs modal
  aneuploidy_association.tsv   per chromosome x predictor: 2x2 counts, Fisher, MH OR, conditional logit
                               by clade, within-clade permutation p, BH q over tested chromosomes

  PYTHONPATH=src .venv/bin/python scripts/peter748.py --config config.peter748.yaml
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

from centromere.aneuploidy import association, chromosome_status, load_peter2018, locus_calls
from centromere.cde import read_reference
from centromere.io import load_config, save_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True)
    parser.add_argument("--skip-pipeline", action="store_true", help="reuse existing filtered.tsv")
    args = parser.parse_args()

    config = load_config(args.config)
    out = Path(config["output"]["dir"])
    cfg = config["aneuploidy"]

    if not args.skip_pipeline:
        for step in ("01_preprocess.py", "02_filter.py"):
            subprocess.run([sys.executable, str(Path(__file__).parent / step), "--config", args.config],
                           check=True)

    meta = load_peter2018(cfg["metadata"])
    parsed = meta.assign(events=meta["events"].map(
        lambda e: "NA" if e is None else ";".join(f"{c}:{n:+d}" for c, n in e)))
    save_tsv(parsed, out / "aneuploidy_parsed.tsv")
    n_eu = int(meta["aneuploidies"].eq("euploid").sum())
    n_g1 = int(meta["aneuploidies"].eq("aneu;+1*1;").sum())
    print(f"Peter 2018: {len(meta)} strains, {n_eu} euploid, {n_g1} 'aneu;+1*1;', "
          f"{meta['events'].isna().sum()} without status")

    loci = pd.read_csv(out / "filtered.tsv", sep="\t")
    loci["strain"] = loci["sample"].str.split("_").str[0]
    calls = locus_calls(loci, read_reference(config["reference"]))
    save_tsv(calls, out / "cen_calls.tsv")

    status = chromosome_status(meta)
    analysed = set(calls["strain"]) & set(status["strain"])
    print(f"BLAST set: {loci['strain'].nunique()} strains, {len(analysed)} with Peter 2018 aneuploidy status")

    def counts(s: pd.DataFrame, label: str) -> pd.DataFrame:
        return (s.assign(gain=s["change"] > 0, loss=s["change"] < 0)
                 .groupby("chrom")[["gain", "loss"]].sum().add_prefix(f"{label}_"))
    table = counts(status, "all").join(counts(status.loc[status["strain"].isin(analysed)], "set"))
    save_tsv(table.reset_index(), out / "aneuploidy_counts.tsv")

    res = association(calls, status, meta, cfg["min_aneuploid"], cfg["n_perm"])
    save_tsv(res, out / "aneuploidy_association.tsv")
    cols = ["centromere", "predictor", "n", "n_aneuploid", "n_variant", "aneu_variant",
            "fisher_or", "mh_or", "clogit_or", "clogit_p", "perm_p", "perm_q", "separation"]
    with pd.option_context("display.width", 200, "display.float_format", "{:.3g}".format):
        print(res.loc[res["tested"], cols].to_string(index=False))
    print(f"Wrote outputs to {out}")


if __name__ == "__main__":
    main()
