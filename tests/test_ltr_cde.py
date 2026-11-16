import numpy as np
import pandas as pd
import pytest

from centromere.cde import CDEModel, read_reference
from centromere.ltr_cde import (CDEScanner, count_contained, dinuc_counts, dinuc_shuffle, extract,
                                free_mask, ltr_intervals, merge, revcomp, sample_windows)

REF = read_reference("data/sacer_centromere.fasta")
CEN1 = REF["CEN1"]
CDEI, CDEII, CDEIII = CEN1[:8], CEN1[8:-26], CEN1[-26:]


@pytest.fixture(scope="module")
def scanner():
    return CDEScanner(CDEModel(REF))


def gc_filler(n, seed=0):
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("GC"), size=n))


def test_unit_found_in_insert(scanner):
    seq = gc_filler(100) + CEN1 + gc_filler(100, 1)
    h = scanner.hits(seq)
    assert len(h.units) == 1
    u = h.units.iloc[0]
    assert (u["strand"], u["lo"], u["hi"]) == (1, 100, 100 + len(CEN1))
    row = scanner.summary(seq)
    assert row["n_units"] == 1 and row["unit_cdei"] == CDEI and row["unit_cdeiii"] == CDEIII


def test_unit_on_reverse_strand(scanner):
    seq = gc_filler(50) + revcomp(CEN1) + gc_filler(50, 1)
    h = scanner.hits(seq)
    assert len(h.units) == 1
    assert (h.units.iloc[0]["strand"], h.units.iloc[0]["lo"]) == (-1, 50)
    row = scanner.summary(seq)
    assert row["unit_strand"] == "-" and row["unit_cdeii"] == CDEII


def test_short_or_gc_rich_cdeii_rejected(scanner):
    short = gc_filler(50) + CDEI + CDEII[:60] + CDEIII + gc_filler(50, 1)
    assert len(scanner.hits(short).units) == 0
    rng = np.random.default_rng(3)
    at60 = "".join(rng.choice(list("ATGC"), p=[0.3, 0.3, 0.2, 0.2], size=len(CDEII)))
    gc_rich = gc_filler(50) + CDEI + at60 + CDEIII + gc_filler(50, 1)
    h = scanner.hits(gc_rich)
    assert len(h.cdei) >= 1 and len(h.cdeiii) >= 1 and len(h.units) == 0


def test_all_s288c_cens_are_units(scanner):
    for name, s in REF.items():
        assert len(scanner.hits(s).units) >= 1, name


def test_dinuc_shuffle_preserves_counts():
    rng = np.random.default_rng(0)
    seq = CEN1 * 3
    for _ in range(20):
        sh = dinuc_shuffle(seq, rng)
        assert dinuc_counts(sh) == dinuc_counts(seq)
        assert sh[0] == seq[0] and sh[-1] == seq[-1]
    assert len({dinuc_shuffle(seq, rng) for _ in range(10)}) > 1


def test_ltr_intervals_full_and_solo():
    el = pd.DataFrame([
        {"genome": "G", "element_id": "G:e0", "contig": "c", "start": 1001, "end": 7000,
         "strand": "+", "family": "Ty1", "class": "full"},
        {"genome": "G", "element_id": "G:e1", "contig": "c", "start": 9001, "end": 9330,
         "strand": "-", "family": "Ty3", "class": "solo_LTR"},
        {"genome": "G", "element_id": "G:e2", "contig": "d", "start": 101, "end": 6000,
         "strand": "-", "family": "Ty2", "class": "full"},
    ])
    hits = pd.DataFrame([  # left LTR of e0 found by an HSP; right end has none
        {"kind": "LTR", "sseqid": "c", "pident": 95, "qcov": 0.95, "s_lo": 1003, "s_hi": 1336,
         "bitscore": 500},
        {"kind": "LTR", "sseqid": "c", "pident": 70, "qcov": 0.95, "s_lo": 6660, "s_hi": 7000,
         "bitscore": 100},   # below min_pident -> ignored
    ])
    iv = ltr_intervals(el, hits, {"Ty1": 334, "Ty2": 332, "Ty3": 340})
    e0 = iv[iv["element_id"] == "G:e0"].set_index("ltr")
    assert (e0.loc["5p", "start"], e0.loc["5p", "end"], e0.loc["5p", "source"]) == (1003, 1336, "hsp")
    assert (e0.loc["3p", "start"], e0.loc["3p", "end"], e0.loc["3p", "source"]) == (6667, 7000, "median_length")
    e1 = iv[iv["element_id"] == "G:e1"].iloc[0]
    assert (e1["ltr"], e1["start"], e1["end"]) == ("solo", 9001, 9330)
    e2 = iv[iv["element_id"] == "G:e2"].set_index("ltr")   # minus strand: left end is 3'
    assert e2.loc["3p", "start"] == 101 and e2.loc["5p", "end"] == 6000


def test_extract_orients():
    assert extract("aaCCGTtt", 3, 6, "+") == "CCGT"
    assert extract("aaCCGTtt", 3, 6, "-") == "ACGG"


def test_windows_avoid_blocked():
    rng = np.random.default_rng(1)
    blocked = {"a": merge(np.array([[100, 900], [850, 950]]))}
    w = sample_windows({"a": 1000, "b": 500}, blocked, np.array([50, 80]), 200, rng)
    assert len(w) == 400 and set(w["rep"]) == set(range(200))
    a = w[w["contig"] == "a"]
    assert free_mask(a["lo"].to_numpy(), a["hi"].to_numpy(), blocked["a"]).all()
    assert (w["hi"] <= w["contig"].map({"a": 1000, "b": 500})).all()


def test_count_contained():
    lo, hi = np.array([10, 20, 30]), np.array([18, 28, 60])
    assert count_contained(lo, hi, np.array([0, 15, 0]), np.array([30, 30, 59])).tolist() == [2, 1, 2]


def test_windows_too_long_dropped_when_allowed():
    rng = np.random.default_rng(2)
    blocked = {"a": np.array([[100, 1000]])}   # only a 100-bp gap is free
    with pytest.raises(RuntimeError):
        sample_windows({"a": 1000}, blocked, np.array([50, 300]), 5, rng, max_rounds=20)
    w = sample_windows({"a": 1000}, blocked, np.array([50, 300]), 5, rng, max_rounds=500, allow_missing=True)
    assert len(w) == 5 and (w["hi"] - w["lo"] == 50).all() and w.attrs["n_missing"] == 5
