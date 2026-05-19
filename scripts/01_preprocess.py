#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.io import load_config, load_data, save_tsv
from centromere.preprocess import find_split_hits, merge_split_hits, remove_gaps


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

    print("Removing gaps...")
    df = remove_gaps(df)

    print("Finding split hits...")
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

    output_path = output_dir / "preprocessed.tsv"
    save_tsv(df, output_path)
    print(f"Saved preprocessed data to {output_path} ({len(df)} rows)")


if __name__ == "__main__":
    main()
