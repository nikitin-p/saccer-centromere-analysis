import pandas as pd


def test_by_length():
    from centromere.filter import by_length

    df = pd.DataFrame({
        "sequence_length": [50, 80, 100, 79],
    })

    result = by_length(df, min_length=80)

    assert len(result) == 2
    assert list(result["sequence_length"]) == [80, 100]


def test_by_gc():
    from centromere.filter import by_gc

    df = pd.DataFrame({
        "sequence": ["AAAA", "ATGC", "GCGC", "ATGC"],
    })

    result = by_gc(df, min_gc=0.05)

    assert len(result) == 3  # AAAA filtered out (GC=0)
    assert "GC" in result.columns


def test_gc_content_calculation():
    from centromere.filter import compute_gc

    assert compute_gc("AAAA") == 0.0
    assert compute_gc("GCGC") == 1.0
    assert compute_gc("ATGC") == 0.5
    assert compute_gc("atgc") == 0.5  # case insensitive
