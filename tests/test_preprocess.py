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


def test_assign_loci_drops_cross_centromere_hits():
    from centromere.preprocess import assign_loci

    # CEN15 query weakly hits the CEN14 locus on c14 (and vice versa);
    # CEN16 is split into two overlapping pieces, both must be kept
    df = pd.DataFrame({
        "sample": ["s1"] * 6,
        "centromere": ["CEN14", "CEN15", "CEN15", "CEN14", "CEN16", "CEN16"],
        "sseqid": ["c14", "c14", "c15", "c15", "c16", "c16"],
        "sstart": [609978, 609992, 311728, 311707, 544678, 544724],
        "send": [610095, 610089, 311610, 311613, 544733, 544816],
        "evalue": [1e-57, 1e-6, 1e-58, 1e-6, 1e-23, 1e-43],
    })

    result = assign_loci(df)

    kept = set(zip(result["centromere"], result["sseqid"]))
    assert kept == {("CEN14", "c14"), ("CEN15", "c15"), ("CEN16", "c16")}
    assert (result["centromere"] == "CEN16").sum() == 2


def test_merge_split_hits_minus_strand():
    from centromere.preprocess import merge_split_hits

    # Same layout as the overlapping test, but on the reverse strand
    df = pd.DataFrame({
        "sample": ["s1", "s1"],
        "centromere": ["CEN1", "CEN1"],
        "sseqid": ["c1", "c1"],
        "qstart": [0, 20],
        "qend": [30, 50],
        "sstart": [150, 130],
        "send": [120, 100],
        "sequence": ["ATCACATGATATATTTTTTATTTTTAATTTT", "TTTAATTTTTTAATTATAAAAATAATTTTTTT"],
        "sequence_length": [31, 32],
        "sequence_gapped": ["ATCACATGATATATTTTTTATTTTTAATTTT", "TTTAATTTTTTAATTATAAAAATAATTTTTTT"],
        "evalue": [1e-10, 1e-10],
    })

    result = merge_split_hits(df, max_distance=0)

    assert len(result) == 1
    assert len(result.iloc[0]["sequence"]) > 31


def test_merge_split_hits_distant_piece_does_not_stretch_coordinates():
    from centromere.preprocess import merge_split_hits

    df = pd.DataFrame({
        "sample": ["s1", "s1"],
        "centromere": ["CEN14", "CEN14"],
        "sseqid": ["c1", "c1"],
        "qstart": [1, 15],
        "qend": [118, 112],
        "sstart": [28457, 287558],
        "send": [28574, 287649],
        "sequence": ["A" * 118, "T" * 92],
        "sequence_length": [118, 92],
        "sequence_gapped": ["A" * 118, "T" * 92],
        "evalue": [1e-57, 1e-6],
    })

    result = merge_split_hits(df, max_distance=5000).iloc[0]

    assert result["send"] == 28574
    assert result["qend"] == 118
