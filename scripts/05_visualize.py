#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.io import load_config
from centromere.viz import (
    plot_gap_histograms,
    plot_gc_violin,
    plot_length_violin,
    plot_mutation_barplots,
    plot_sequence_logos,
)


def main():
    parser = argparse.ArgumentParser(description="Generate visualizations")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(config["output"]["dir"])
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    centromere_order = config["visualize"].get("centromere_order")
    formats = config["visualize"].get("formats", ["png"])

    filtered_path = output_dir / "filtered.tsv"
    aligned_path = output_dir / "aligned.tsv"
    mutation_path = output_dir / "mutation_stats.tsv"

    print("Loading data...")
    filtered_df = pd.read_csv(filtered_path, sep="\t")
    aligned_df = pd.read_csv(aligned_path, sep="\t")
    mutation_df = pd.read_csv(mutation_path, sep="\t")

    for fmt in formats:
        print(f"Generating {fmt.upper()} figures...")

        print("  Length violin plot...")
        plot_length_violin(
            filtered_df,
            figures_dir / f"length_violin.{fmt}",
            centromere_order,
        )

        print("  GC violin plot...")
        plot_gc_violin(
            filtered_df,
            figures_dir / f"gc_violin.{fmt}",
            centromere_order,
        )

    print("Generating sequence logos...")
    plot_sequence_logos(aligned_df, figures_dir / "logos", centromere_order)

    print("Generating gap histograms...")
    plot_gap_histograms(aligned_df, figures_dir / "gap_histograms")

    print("Generating mutation barplots...")
    plot_mutation_barplots(mutation_df, figures_dir / "mutations")

    print(f"All figures saved to {figures_dir}")


if __name__ == "__main__":
    main()
