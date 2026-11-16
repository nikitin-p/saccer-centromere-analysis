from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from centromere.aneuploidy import (association, chromosome_status, load_peter2018, locus_calls,
                                   normalize_clade, parse_aneuploidies, stratified_permutation)

PETER = Path(__file__).parents[1] / "data" / "41586_2018_30_MOESM3_ESM.xls"


def test_parse_aneuploidies():
    assert parse_aneuploidies("euploid") == []
    assert parse_aneuploidies("aneu;+1*9;") == [(9, 1)]
    assert parse_aneuploidies("aneu;+2*2;-1*3;+1*16;") == [(2, 2), (3, -1), (16, 1)]
    assert parse_aneuploidies(np.nan) is None


def test_normalize_clade():
    assert normalize_clade("1. Wine/European (subclade 3)") == "1. Wine/European"
    assert normalize_clade("M3. Mosaic region 3 ") == "M3. Mosaic region 3"
    assert normalize_clade(np.nan) == "unassigned"


@pytest.mark.skipif(not PETER.exists(), reason="Peter 2018 table not present")
def test_parser_reproduces_peter2018_counts():
    meta = load_peter2018(PETER)
    assert (meta["events"].map(lambda e: e == [])).sum() == 793
    assert (meta["events"].map(lambda e: e == [(1, 1)])).sum() == 31


def test_chromosome_status():
    meta = pd.DataFrame({"strain": ["A", "B", "C"], "events": [[], [(9, 1), (3, -1)], None]})
    st = chromosome_status(meta).set_index(["strain", "chrom"])["change"]
    assert st[("B", 9)] == 1 and st[("B", 3)] == -1 and st[("A", 9)] == 0
    assert "C" not in st.index.get_level_values(0)


def _loci(rows):
    return pd.DataFrame(rows, columns=["strain", "centromere", "qstart", "qend", "sequence"]).assign(
        sequence_length=lambda d: d["sequence"].str.len())


def test_locus_calls_trimmed_end_and_variants():
    ref = {"CEN1": "GAAAATTTTTCC"}
    loci = _loci([
        ("S288C", "CEN1", 1, 12, "GAAAATTTTTCC"),   # reference
        ("A", "CEN1", 2, 12, "AAAATTTTTCC"),        # base 1 mismatch trimmed: modal, non-S288C
        ("B", "CEN1", 2, 12, "AAAATTTTTCC"),
        ("C", "CEN1", 2, 12, "AAAATTATTTTCC"),      # +2 bp
        ("D", "CEN1", 5, 12, "ATTTTTCC"),           # truncated: not called
    ])
    c = locus_calls(loci, ref).set_index("strain")
    assert c.loc["S288C", "non_s288c"] == 0 and c.loc["A", "non_s288c"] == 1
    assert c.loc["A", "non_modal"] == 0 and c.loc["S288C", "non_modal"] == 1
    assert c.loc["A", "length_variant"] == 0 and c.loc["C", "delta_len"] == 2
    assert c.loc["D", "n_complete"] == 0 and np.isnan(c.loc["D", "non_modal"])


def test_permutation_detects_within_clade_association_only():
    rng = np.random.default_rng(1)
    strata = np.repeat(["a", "b"], 100)
    x = rng.integers(0, 2, 200)
    y_assoc = (x & (rng.random(200) < 0.6)).astype(int)
    assert stratified_permutation(y_assoc, x, strata, 2000) < 0.01
    # confounded by clade: x and y both 1 in clade b only, no association inside a clade
    x_conf = np.r_[np.zeros(100, int), np.ones(100, int)]
    y_conf = np.r_[rng.random(100) < 0.05, rng.random(100) < 0.5].astype(int)
    assert stratified_permutation(y_conf, x_conf, strata, 2000) == 1.0


def test_association_runs_and_flags_tested():
    rng = np.random.default_rng(2)
    n = 200
    strains = [f"S{i}" for i in range(n)]
    x = rng.integers(0, 2, n)
    y = ((x & (rng.random(n) < 0.5)) | (rng.random(n) < 0.05)).astype(int)
    calls = pd.DataFrame({"strain": strains, "centromere": "CEN9", "chrom": 9, "n_complete": 1,
                          "non_modal": x, "length_variant": x, "non_s288c": x})
    status = pd.DataFrame({"strain": strains, "chrom": 9, "change": y})
    meta = pd.DataFrame({"strain": strains, "clade": np.where(np.arange(n) % 2, "a", "b")})
    res = association(calls, status, meta, min_aneuploid=10, n_perm=2000)
    row = res.loc[res["predictor"] == "non_modal"].iloc[0]
    assert row["tested"] and row["perm_p"] < 0.01 and row["clogit_or"] > 1
