#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.cde import annotate_cde, read_reference
from centromere.io import annotate_samples, load_config, load_data, save_tsv
from centromere.locate import apply_reconcile, locator_loci, reconcile, rescue_rows
from centromere.preprocess import assign_loci, find_split_hits, merge_split_hits, primary_loci, remove_gaps


def reconcile_locators(df, config, reference, output_dir):
    """Locator spec phase 4: merge L1 (this table) with the L2/L3/L4 loci already on disk;
    write loci_reconciled.tsv; with locate.rescue, rescued / insertion loci enter the table."""
    loc = config["locate"]
    loci, ran = locator_loci(config, df)
    print(f"Reconciling locators {sorted(t.upper() for t in loci)}...")
    rec = reconcile(loci, ran, {c: len(s) for c, s in reference.items()}, loc.get("insertion_bp", 30))
    save_tsv(rec, output_dir / "loci_reconciled.tsv")
    print("  " + ", ".join(f"{k} {v}" for k, v in rec["category"].value_counts().items()))
    if not loc.get("rescue", True):
        return df.assign(locator="L1").merge(rec[["sample", "centromere", "support"]],
                                             on=["sample", "centromere"], how="left")
    strains = df.dropna(subset=["strain"]).groupby("sample")["strain"].first().to_dict()
    resc = rescue_rows(rec, config["cde"]["genome_dir"], reference, strains,
                       config["cde"].get("flank", 150), config["cde"].get("search_radius", 20))
    if len(resc):
        print(f"  {len(resc)} loci from L2/L3/L4 enter the table: "
              + ", ".join(f"{r.sample} {r.centromere} ({r.locator})" for r in resc.itertuples()))
    return apply_reconcile(df, rec, resc)


def main():
    parser = argparse.ArgumentParser(description="Preprocess centromere data")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(config["output"]["dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data from {config['input']['path']}...")
    df = load_data(config["input"]["path"], config["input"]["format"])
    print(f"  Loaded {len(df)} rows")

    sample_map = config["input"].get("sample_map")
    if sample_map:
        df = annotate_samples(df, sample_map)
        unmapped = df.loc[df["strain"].isna(), "sample"].unique()
        print(f"  Annotated strains from {sample_map} ({len(unmapped)} samples unmapped)")

    print("Removing gaps...")
    df = remove_gaps(df)
    qc = df.groupby("sample").size().rename("raw_hits").to_frame()

    print("Assigning loci to best-matching centromere...")
    n_before = len(df)
    df = assign_loci(df)
    print(f"  Dropped {n_before - len(df)} cross-centromere hits")
    qc["after_locus_assignment"] = df.groupby("sample").size()

    print("Finding split hits...")
    df["blast_pieces"] = df.groupby(["sample", "centromere", "sseqid"])["sseqid"].transform("size")
    split_df, ok_df = find_split_hits(df)
    print(f"  Found {len(split_df)} split hits, {len(ok_df)} ok hits")

    if len(split_df) > 0:
        print("Merging split hits...")
        max_dist = config["preprocess"]["merge_max_distance"]
        merged_df = merge_split_hits(split_df, max_distance=max_dist)
        print(f"  Merged into {len(merged_df)} rows")

        df = pd.concat([ok_df, merged_df], ignore_index=True)
    else:
        df = ok_df

    df, secondary = primary_loci(df)
    qc["secondary_loci"] = secondary.groupby("sample").size()
    qc["secondary_loci"] = qc["secondary_loci"].fillna(0).astype(int)
    save_tsv(secondary, output_dir / "secondary_loci.tsv")
    if len(secondary):
        print(f"  {len(secondary)} extra loci of an already-found centromere set aside "
              f"(secondary_loci.tsv): " + ", ".join(f"{r.sample} {r.centromere} {r.sseqid}"
                                                    for r in secondary.head(10).itertuples()))

    genome_dir = config.get("cde", {}).get("genome_dir")
    if genome_dir:
        print(f"Annotating CDEI/CDEII/CDEIII from assemblies in {genome_dir}...")
        reference = read_reference(config["reference"])
        df = annotate_cde(
            df, genome_dir, reference,
            flank=config["cde"].get("flank", 150),
            radius=config["cde"].get("search_radius", 20),
        )
        print(f"  {df['cde_ok'].sum()}/{len(df)} loci with confident CDEI and CDEIII")

    loc = config.get("locate", {})
    if genome_dir and loc.get("reconcile"):
        df = reconcile_locators(df, config, reference, output_dir)

    qc["final_loci"] = df.groupby("sample").size()
    save_tsv(qc.reset_index(), output_dir / "preprocess_qc.tsv")

    output_path = output_dir / "preprocessed.tsv"
    save_tsv(df, output_path)
    print(f"Saved preprocessed data to {output_path} ({len(df)} rows)")


if __name__ == "__main__":
    main()
