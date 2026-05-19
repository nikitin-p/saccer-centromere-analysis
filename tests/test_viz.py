import os
import tempfile

import pandas as pd


def test_plot_length_violin():
    from centromere.viz import plot_length_violin

    df = pd.DataFrame({
        "centromere": ["CEN1"] * 10 + ["CEN2"] * 10,
        "sequence_length": list(range(80, 90)) + list(range(90, 100)),
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "length_violin.png")
        plot_length_violin(df, output_path)
        assert os.path.exists(output_path)


def test_plot_gc_violin():
    from centromere.viz import plot_gc_violin

    df = pd.DataFrame({
        "centromere": ["CEN1"] * 10 + ["CEN2"] * 10,
        "GC": [0.3 + i * 0.01 for i in range(20)],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "gc_violin.png")
        plot_gc_violin(df, output_path)
        assert os.path.exists(output_path)


def test_plot_sequence_logos():
    from centromere.viz import plot_sequence_logos

    df = pd.DataFrame({
        "centromere": ["CEN1"] * 5,
        "aligned_seq": ["ATCG", "ATCG", "ATGG", "ATCG", "ATCG"],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        plot_sequence_logos(df, tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "CEN1_logo.png"))


def test_plot_gap_histograms():
    from centromere.viz import plot_gap_histograms

    df = pd.DataFrame({
        "centromere": ["CEN1"] * 3,
        "aligned_seq": ["AT-G", "ATCG", "AT-G"],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        plot_gap_histograms(df, tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "CEN1_gaps.png"))


def test_plot_mutation_barplots():
    from centromere.viz import plot_mutation_barplots

    mutation_df = pd.DataFrame({
        "centromere": ["CEN1", "CEN1", "CEN2"],
        "mutation": ["A>T", "A>G", "A>T"],
        "count": [10, 5, 8],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        plot_mutation_barplots(mutation_df, tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "stacked_mutation_barplot.png"))
