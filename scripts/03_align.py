#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.align import align_all_parallel
from centromere.io import load_config, save_tsv


def main():
    parser = argparse.ArgumentParser(description="Align centromere sequences")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(config["output"]["dir"])

    input_path = output_dir / "filtered.tsv"
    print(f"Loading filtered data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")
    print(f"  Loaded {len(df)} rows")

    mafft_args = config["align"]["mafft_args"]
    workers = config["align"]["workers"]
    print(f"Aligning sequences with MAFFT ({workers} workers)...")

    df = align_all_parallel(df, args=mafft_args, workers=workers)
    print(f"  Alignment complete")

    output_path = output_dir / "aligned.tsv"
    save_tsv(df, output_path)
    print(f"Saved aligned data to {output_path}")


if __name__ == "__main__":
    main()
