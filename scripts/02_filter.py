#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.filter import by_gc, by_length
from centromere.io import load_config, save_tsv


def main():
    parser = argparse.ArgumentParser(description="Filter centromere data")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(config["output"]["dir"])

    input_path = output_dir / "preprocessed.tsv"
    print(f"Loading preprocessed data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")
    print(f"  Loaded {len(df)} rows")

    min_length = config["filter"]["min_length"]
    print(f"Filtering by length >= {min_length}...")
    df = by_length(df, min_length)
    print(f"  {len(df)} rows remaining")

    min_gc = config["filter"]["min_gc"]
    print(f"Filtering by GC > {min_gc}...")
    df = by_gc(df, min_gc)
    print(f"  {len(df)} rows remaining")

    output_path = output_dir / "filtered.tsv"
    save_tsv(df, output_path)
    print(f"Saved filtered data to {output_path}")


if __name__ == "__main__":
    main()
