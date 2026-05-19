import pandas as pd


def test_remove_gaps():
    from centromere.preprocess import remove_gaps

    df = pd.DataFrame({
        "sequence_gapped": ["ATCG--AT", "GC-TTA--A"],
    })

    result = remove_gaps(df)

    assert "sequence" in result.columns
    assert "sequence_length" in result.columns
    assert result.iloc[0]["sequence"] == "ATCGAT"
    assert result.iloc[0]["sequence_length"] == 6
    assert result.iloc[1]["sequence"] == "GCTTAA"
    assert result.iloc[1]["sequence_length"] == 6


def test_find_split_hits():
    from centromere.preprocess import find_split_hits

    df = pd.DataFrame({
        "sample": ["s1", "s1", "s1", "s2"],
        "centromere": ["CEN1", "CEN1", "CEN2", "CEN1"],
        "sseqid": ["c1", "c1", "c1", "c1"],
        "sstart": [100, 200, 300, 100],
    })

    split_df, ok_df = find_split_hits(df)

    assert len(split_df) == 2  # s1/CEN1/c1 appears twice
    assert len(ok_df) == 2     # s1/CEN2/c1 and s2/CEN1/c1


def test_merge_split_hits_overlapping():
    from centromere.preprocess import merge_split_hits

    df = pd.DataFrame({
        "sample": ["s1", "s1"],
        "centromere": ["CEN1", "CEN1"],
        "sseqid": ["c1", "c1"],
        "qstart": [0, 20],
        "qend": [30, 50],
        "sstart": [100, 120],
        "send": [130, 150],
        "sequence": ["ATCACATGATATATTTTTTATTTTTAATTTT", "TTTAATTTTTTAATTATAAAAATAATTTTTTT"],
        "sequence_length": [31, 32],
        "sequence_gapped": ["ATCACATGATATATTTTTTATTTTTAATTTT", "TTTAATTTTTTAATTATAAAAATAATTTTTTT"],
        "evalue": [1e-10, 1e-10],
    })

    result = merge_split_hits(df, max_distance=5000)

    assert len(result) == 1
    assert result.iloc[0]["sample"] == "s1"
    assert len(result.iloc[0]["sequence"]) > 31  # merged sequence is longer
