#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.io import load_config, save_tsv
from centromere.mutations import (
    compute_consensus,
    compute_mutation_patterns,
    count_bases_per_position,
    identify_variable_sites,
)


def main():
    parser = argparse.ArgumentParser(description="Analyze mutation patterns")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(config["output"]["dir"])

    input_path = output_dir / "aligned.tsv"
    print(f"Loading aligned data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")
    print(f"  Loaded {len(df)} rows")

    print("Computing mutation patterns...")
    mutation_df = compute_mutation_patterns(df)

    print("Computing consensus sequences...")
    consensus_data = []
    for centromere, group in df.groupby("centromere"):
        sequences = group["aligned_seq"].tolist()
        consensus = compute_consensus(sequences)
        counts = count_bases_per_position(sequences)
        variable = identify_variable_sites(counts)
        consensus_data.append({
            "centromere": centromere,
            "consensus": consensus,
            "length": len(consensus),
            "variable_sites": len(variable),
            "total_sequences": len(sequences),
        })
    consensus_df = pd.DataFrame(consensus_data)

    mutation_path = output_dir / "mutation_stats.tsv"
    save_tsv(mutation_df, mutation_path)
    print(f"Saved mutation stats to {mutation_path}")

    consensus_path = output_dir / "consensus.tsv"
    save_tsv(consensus_df, consensus_path)
    print(f"Saved consensus data to {consensus_path}")


if __name__ == "__main__":
    main()
