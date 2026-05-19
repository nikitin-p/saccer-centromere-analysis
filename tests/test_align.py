import shutil

import pandas as pd
import pytest


@pytest.mark.skipif(
    shutil.which("mafft") is None,
    reason="mafft not installed"
)
def test_run_mafft():
    from centromere.align import run_mafft

    sequences = ["ATCGATCG", "ATCGATCG", "ATCGTTCG"]

    aligned = run_mafft(sequences, args="--auto")

    assert len(aligned) == 3
    assert all(len(s) == len(aligned[0]) for s in aligned)


@pytest.mark.skipif(
    shutil.which("mafft") is None,
    reason="mafft not installed"
)
def test_align_centromere_group():
    from centromere.align import align_centromere_group

    df = pd.DataFrame({
        "sample": ["s1", "s2", "s3"],
        "centromere": ["CEN1", "CEN1", "CEN1"],
        "sequence": ["ATCGATCG", "ATCGATCG", "ATCGTTCG"],
    })

    result = align_centromere_group(df, "CEN1", "--auto")

    assert "aligned_seq" in result.columns
    assert len(result) == 3


@pytest.mark.skipif(
    shutil.which("mafft") is None,
    reason="mafft not installed"
)
def test_align_all_parallel():
    from centromere.align import align_all_parallel

    df = pd.DataFrame({
        "sample": ["s1", "s2", "s3", "s4"],
        "centromere": ["CEN1", "CEN1", "CEN2", "CEN2"],
        "sequence": ["ATCGATCG", "ATCGTTCG", "GCTAGCTA", "GCTTGCTA"],
    })

    result = align_all_parallel(df, args="--auto", workers=2)

    assert "aligned_seq" in result.columns
    assert len(result) == 4
