import pandas as pd


def test_count_bases_per_position():
    from centromere.mutations import count_bases_per_position

    sequences = ["ATCG", "ATCG", "ATGG"]

    counts = count_bases_per_position(sequences)

    assert len(counts) == 4
    assert counts[0] == {"A": 3}
    assert counts[1] == {"T": 3}
    assert counts[2] == {"C": 2, "G": 1}
    assert counts[3] == {"G": 3}


def test_compute_consensus():
    from centromere.mutations import compute_consensus

    sequences = ["ATCG", "ATCG", "ATGG"]

    consensus = compute_consensus(sequences)

    assert consensus == "ATCG"  # majority at each position


def test_identify_variable_sites():
    from centromere.mutations import identify_variable_sites

    counts = [
        {"A": 3},           # position 0: not variable
        {"T": 2, "C": 1},   # position 1: variable
        {"G": 3},           # position 2: not variable
        {"A": 1, "T": 1, "G": 1},  # position 3: variable
    ]

    variable = identify_variable_sites(counts)

    assert variable == [1, 3]


def test_compute_mutation_patterns():
    from centromere.mutations import compute_mutation_patterns

    df = pd.DataFrame({
        "centromere": ["CEN1", "CEN1"],
        "aligned_seq": ["ATCG", "ATGG"],
    })

    result = compute_mutation_patterns(df)

    assert "mutation" in result.columns
    assert "count" in result.columns
    assert len(result) > 0
