import os
from pathlib import Path

import logomaker
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import seaborn as sns


def plot_length_violin(
    df: pd.DataFrame,
    output_path: str | Path,
    centromere_order: list[str] | None = None,
) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    df = df.copy()
    if centromere_order:
        df = df[df["centromere"].isin(centromere_order)]
        df["centromere"] = pd.Categorical(
            df["centromere"], categories=centromere_order, ordered=True
        )

    fig = px.violin(
        df,
        y="sequence_length",
        x="centromere",
        color="centromere",
        box=False,
        points="all",
    )
    fig.update_traces(orientation="v", side="positive", width=0.5)
    fig.write_image(str(output_path))
    plt.close("all")


def plot_gc_violin(
    df: pd.DataFrame,
    output_path: str | Path,
    centromere_order: list[str] | None = None,
) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    df = df.copy()
    if centromere_order:
        df = df[df["centromere"].isin(centromere_order)]
        df["centromere"] = pd.Categorical(
            df["centromere"], categories=centromere_order, ordered=True
        )

    fig = px.violin(
        df,
        y="GC",
        x="centromere",
        color="centromere",
        box=False,
        points="all",
    )
    fig.update_traces(orientation="v", side="positive", width=0.5)
    fig.write_image(str(output_path))
    plt.close("all")


def _sequences_to_counts_df(sequences: list[str]) -> pd.DataFrame:
    if not sequences:
        return pd.DataFrame()

    length = len(sequences[0])
    bases = ["A", "C", "G", "T", "-"]

    data = {base: [0] * length for base in bases}

    for seq in sequences:
        for i, char in enumerate(seq.upper()):
            if i < length and char in data:
                data[char][i] += 1

    df = pd.DataFrame(data)
    row_sums = df.sum(axis=1)
    df = df.div(row_sums, axis=0).fillna(0)
    return df


def plot_sequence_logos(
    df: pd.DataFrame,
    output_dir: str | Path,
    centromere_order: list[str] | None = None,
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    centromeres = centromere_order or df["centromere"].unique()

    for centromere in centromeres:
        group = df[df["centromere"] == centromere]
        if group.empty:
            continue

        sequences = group["aligned_seq"].tolist()
        if not sequences:
            continue

        counts_df = _sequences_to_counts_df(sequences)
        if counts_df.empty:
            continue

        fig, ax = plt.subplots(figsize=(max(10, len(counts_df) * 0.1), 3))
        logomaker.Logo(counts_df, ax=ax)
        ax.set_title(centromere)
        fig.savefig(output_dir / f"{centromere}_logo.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


def plot_gap_histograms(
    df: pd.DataFrame,
    output_dir: str | Path,
) -> None:
    from centromere.mutations import count_bases_per_position

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for centromere, group in df.groupby("centromere"):
        sequences = group["aligned_seq"].tolist()
        counts = count_bases_per_position(sequences)

        gaps = [c.get("-", 0) for c in counts]
        positions = list(range(len(gaps)))

        fig, ax = plt.subplots(figsize=(max(10, len(gaps) * 0.05), 4))
        ax.bar(positions, gaps, width=1.0)
        ax.set_xlabel("Position")
        ax.set_ylabel("Gap count")
        ax.set_title(f"{centromere} - Gap distribution")
        fig.savefig(output_dir / f"{centromere}_gaps.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


def plot_mutation_barplots(
    mutation_df: pd.DataFrame,
    output_dir: str | Path,
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pivot = mutation_df.pivot_table(
        index="mutation",
        columns="centromere",
        values="count",
        fill_value=0,
    )

    fig, ax = plt.subplots(figsize=(12, 6))
    pivot.plot(kind="bar", stacked=True, ax=ax)
    ax.set_xlabel("Mutation type")
    ax.set_ylabel("Count")
    ax.set_title("Mutation patterns by centromere")
    ax.legend(title="Centromere", bbox_to_anchor=(1.05, 1))
    fig.savefig(
        output_dir / "stacked_mutation_barplot.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig)
