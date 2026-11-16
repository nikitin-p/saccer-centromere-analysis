import shutil
import subprocess

import numpy as np
import pytest
from scipy.stats import kstest

from centromere.grammar import (GenomeWindows, aatt_power, atracts, empirical_p, parse_shape_file,
                                percentile, phase_score, run_dnashape, shape_distance, spectrum_null,
                                substitutions)
from centromere.ltr_cde import dinuc_shuffle


def test_atracts_hand_cases():
    assert atracts("GAAAATTTTG") == [(1, 8)]
    assert atracts("GTTTTAAAAG") == [(1, 4), (5, 4)]
    assert atracts("AAATAAA") == [(0, 4)]          # AAAT, then AAA is too short
    assert atracts("AAAGAAAA") == [(4, 4)]          # G breaks
    assert atracts("ATATATAT") == []                # TpA every other step
    assert atracts("AAAAAA", min_len=7) == []


def test_phase_score():
    assert phase_score([(0, 5), (21, 5)]) == 1.0     # 21 = 2 x 10.5
    assert phase_score([(0, 5), (15, 5)]) == 0.0
    assert np.isnan(phase_score([(0, 5)]))


def _phased(rng, n=90):
    s = list(rng.choice(list("GCATGC"), n))
    for k in range(0, n - 6, 10):
        s[k:k + 5] = "AAAAA"
    return "".join(s)


def test_phasing_beats_shuffles():
    rng = np.random.default_rng(0)
    seq = _phased(rng)
    obs = phase_score(atracts(seq, 4), 10.0)
    null = [phase_score(atracts(dinuc_shuffle(seq, rng), 4), 10.0) for _ in range(200)]
    assert obs > np.nanpercentile(null, 95)
    assert aatt_power(seq, (10.0,)) > 5


def test_empirical_p_and_percentile():
    null = np.arange(100.0)
    assert empirical_p(null, 1000) == pytest.approx(2 / 101)
    assert percentile(null, 49.5) == 0.5
    assert np.isnan(percentile(np.array([np.nan]), 1.0))


def test_genome_windows_match_at_and_skip_blocked():
    rng = np.random.default_rng(1)
    genome = {"c1": "A" * 50 + "GC" * 50 + "A" * 50, "c2": "N" * 30 + "AT" * 40}
    gw = GenomeWindows(genome, {"c1": [(0, 50)]})
    seqs, tol = gw.sample(10, 1.0, 0.0, 1000, rng)
    assert tol == 0.0 and len(seqs) > 0
    assert all(set(s) <= set("AT") for s in seqs)
    # c1 0-50 blocked, c2 has N in 0-30: all-AT windows come from c1 tail or c2 >= 30
    assert len(seqs) == (150 - 100 - 10 + 1) + (80 - 10 + 1)
    _, tol = gw.sample(10, 0.5, 0.0, 5, rng, tol_wide=0.2, min_n=10 ** 6)
    assert tol == 0.2


def test_spectrum_null_keeps_substitutions():
    rng = np.random.default_rng(2)
    ref = "AATTGCAATTAAGCTTAA"
    alt = "AGTTGCAATTAAGCTTTA"
    subs = substitutions(ref, alt)
    assert subs == [(1, "A", "G"), (16, "A", "T")]
    for m in spectrum_null(ref, subs, 50, rng):
        s = sorted((a, b) for _, a, b in substitutions(ref, m))
        assert s == [("A", "G"), ("A", "T")]


def test_null_calibration_uniform():
    """Check 3: an 'observed' allele drawn from the null gets a uniform percentile."""
    rng = np.random.default_rng(0)
    ref = "".join(rng.choice(list("AAATTTGC"), 85))

    def toy_shape(s):  # pentamer AT count, a stand-in for MGW
        a = np.array([b in "AT" for b in s], float)
        return np.convolve(a, np.ones(5), "valid")

    pct = []
    for _ in range(300):
        n = int(rng.integers(1, 6))
        idx = rng.choice(len(ref), n, replace=False)
        subs = [(int(i), ref[i], rng.choice([b for b in "ACGT" if b != ref[i]])) for i in idx]
        obs, *null = spectrum_null(ref, subs, 101, rng)
        d0 = shape_distance(toy_shape(ref), toy_shape(obs))
        dn = [shape_distance(toy_shape(ref), toy_shape(m)) for m in null]
        dn = np.array(dn)  # randomised tie-breaking: the toy distances are very discrete
        pct.append((np.sum(dn < d0) + rng.uniform() * (np.sum(dn == d0) + 1)) / (len(dn) + 1))
    assert kstest(pct, "uniform").pvalue > 0.01


def test_parse_shape_file(tmp_path):
    p = tmp_path / "x.fa.Roll"
    p.write_text(">s0\nNA,1.5,2.0,\n-3.0,NA\n>s1\nNA,4.0,NA\n")
    r = parse_shape_file(p)
    assert np.allclose(r["s0"][1:4], [1.5, 2.0, -3.0]) and len(r["s0"]) == 5
    assert len(r["s1"]) == 3


def _has_dnashaper():
    if not shutil.which("Rscript"):
        return False
    r = subprocess.run(["Rscript", "-e", "library(DNAshapeR)"], capture_output=True)
    return r.returncode == 0


@pytest.mark.skipif(not _has_dnashaper(), reason="DNAshapeR not installed")
def test_shape_sanity_atract_narrow_groove(tmp_path):
    """Check 2: an A-tract narrows the minor groove relative to GC DNA."""
    sh = run_dnashape({"at": "CGCAAAAAAGCG", "gc": "CGCGCGCGCGCG"}, tmp_path)
    assert sh["MGW"]["at"][6] < sh["MGW"]["gc"][6] - 1.0
    assert len(sh["Roll"]["at"]) == 11


def test_event_null_positions():
    from centromere.grammar import event_null_positions, mutate
    ref = "AATAGATAC"
    assert event_null_positions(ref, 0) == [1, 3, 5, 7]
    assert event_null_positions(ref, 2, context=True) == [6]   # ATA at 1-3 and 5-7
    assert event_null_positions(ref, 0, context=True) == []    # terminal base: no context
    assert mutate(ref, 2, "G") == "AAGAGATAC"
