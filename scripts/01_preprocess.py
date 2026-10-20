#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.cde import annotate_cde, read_reference
from centromere.io import annotate_samples, load_config, load_data, save_tsv
from centromere.preprocess import assign_loci, find_split_hits, merge_split_hits, remove_gaps


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

    qc["final_loci"] = df.groupby("sample").size()
    save_tsv(qc.reset_index(), output_dir / "preprocess_qc.tsv")

    output_path = output_dir / "preprocessed.tsv"
    save_tsv(df, output_path)
    print(f"Saved preprocessed data to {output_path} ({len(df)} rows)")


if __name__ == "__main__":
    main()
