import numpy as np
import pandas as pd
import pytest

from centromere.pheno import (cen_predictors, cis_genes, eigen, fit_null, genomic_lambda,
                              inverse_normal, lmm_test, loco_kinships, ols_test, rare_burden,
                              read_plink, write_plink)


# ------------------------------------------------------------- plink round trip

def test_plink_roundtrip(tmp_path):
    rng = np.random.default_rng(1)
    geno = rng.choice([-1, 0, 1, 2], size=(7, 5), p=[0.1, 0.3, 0.3, 0.3]).astype(np.int8)
    bim = pd.DataFrame({"chrom": [1, 1, 2, 2, 17], "id": [f"rs{i}" for i in range(5)], "cm": 0,
                        "pos": [10, 20, 30, 40, 50], "a1": "A", "a2": "G"})
    samples = [f"S{i}" for i in range(7)]  # 7 samples -> padding in the last byte
    write_plink(tmp_path / "x", samples, bim, geno)
    ids, bim2, g2 = read_plink(tmp_path / "x")
    assert ids == samples
    assert bim2["id"].tolist() == bim["id"].tolist()
    np.testing.assert_array_equal(g2, geno)


def test_plink_byte_decode(tmp_path):
    # one marker, 4 samples: codes 00 (hom A1), 01 (missing), 10 (het), 11 (hom A2) in one byte
    (tmp_path / "y.fam").write_text("".join(f"s{i} s{i} 0 0 0 -9\n" for i in range(4)))
    (tmp_path / "y.bim").write_text("1\trs1\t0\t1\tA\tG\n")
    byte = 0b11_10_01_00
    (tmp_path / "y.bed").write_bytes(bytes([0x6C, 0x1B, 0x01, byte]))
    _, _, g = read_plink(tmp_path / "y")
    assert g[:, 0].tolist() == [2, -1, 1, 0]


# ------------------------------------------------------------- LMM calibration

def _structured_population(rng, n=600, n_clades=6, m=1500):
    clade = np.repeat(np.arange(n_clades), n // n_clades)
    freq = rng.uniform(0.05, 0.95, size=(n_clades, m))
    geno = rng.binomial(2, freq[clade]).astype(np.int8)
    k = loco_kinships(geno, np.ones(m, dtype=int))["all"]
    return clade, k


def _simulate(rng, clade, k, effect):
    n = len(clade)
    # predictor tracks clade, trait has a clade-driven polygenic component
    pred = (rng.uniform(size=n) < np.linspace(0.1, 0.9, clade.max() + 1)[clade]).astype(float)
    g = np.linalg.cholesky(k + 1e-6 * np.eye(n)) @ rng.standard_normal(n)
    y = 1.5 * g + rng.standard_normal(n) + effect * pred
    return y, pred


@pytest.fixture(scope="module")
def population():
    rng = np.random.default_rng(7)
    clade, k = _structured_population(rng)
    return clade, k, eigen(k)


def test_lmm_null_calibrated_ols_inflated(population):
    clade, k, (s, u) = population
    rng = np.random.default_rng(11)
    x = np.ones((len(clade), 1))
    p_lmm, p_ols = [], []
    for _ in range(200):
        y, pred = _simulate(rng, clade, k, effect=0.0)
        p_lmm.append(lmm_test(fit_null(y, x, s, u), pred)["p"])
        p_ols.append(ols_test(y, x, pred))
    fpr_lmm = np.mean(np.array(p_lmm) < 0.05)
    fpr_ols = np.mean(np.array(p_ols) < 0.05)
    assert 0.03 <= fpr_lmm <= 0.07, fpr_lmm
    assert fpr_ols > 0.2, fpr_ols


def test_lmm_recovers_effect(population):
    clade, k, (s, u) = population
    rng = np.random.default_rng(3)
    x = np.ones((len(clade), 1))
    hits = 0
    for _ in range(20):
        y, pred = _simulate(rng, clade, k, effect=0.5)
        r = lmm_test(fit_null(y, x, s, u), pred)
        hits += abs(r["beta"] - 0.5) < 2 * r["se"]
    assert hits >= 17


def test_fit_null_h2_without_kinship_signal():
    rng = np.random.default_rng(5)
    n = 300
    k = np.eye(n)
    s, u = eigen(k + 0.0)
    y = rng.standard_normal(n)
    nl = fit_null(y, np.ones((n, 1)), s, u)
    assert 0 <= nl["h2"] <= 1


# ------------------------------------------------------------- helpers

def test_inverse_normal_and_lambda():
    z = inverse_normal(np.array([3.0, 1.0, 2.0, 10.0]))
    assert np.all(np.diff(z[[1, 2, 0, 3]]) > 0)
    rng = np.random.default_rng(0)
    assert genomic_lambda(pd.Series(rng.uniform(size=20000))) == pytest.approx(1, abs=0.05)


def _cen(cdeii):
    return "ATCACGTG" + cdeii + "TGTTTTGTTTTCCGAAGATGTAAAA" + "A" + "A"  # CDEIII 26 + 1 extra


def test_cen_predictors():
    modal = "A" * 80
    rows = [("S1", "A", modal), ("S2", "A", modal), ("S3", "A", modal), ("S4", "A", "A" * 90),
            ("S5", "A", "A" * 85), ("S5", "B", modal)]
    h = pd.DataFrame([{"assembly": a, "centromere": "CEN1", "variant": v, "sequence": _cen(c)}
                      for a, v, c in rows])
    p = cen_predictors(h, rare_freq=0.2).set_index("strain")
    assert p.loc["S1", "non_modal"] == 0
    assert p.loc["S4", "non_modal"] == 1 and p.loc["S4", "cdeii_delta"] == 10
    assert p.loc["S4", "off_step"] == 0          # +10 = one helical turn
    assert p.loc["S5", "non_modal"] == 0.5 and p.loc["S5", "het"] == 1
    assert p.loc["S5", "off_step"] == 1          # +5 = half a turn
    assert p.loc["S4", "rare"] == 1 and p.loc["S1", "rare"] == 0
    assert rare_burden(p.reset_index(), n_cen=1).loc["S4"] == 1


def test_cis_genes():
    genes = pd.DataFrame({"chrom": ["chrI"] * 3, "start": [100, 15000, 40000],
                          "end": [900, 16000, 41000], "orf": ["A", "B", "C"], "name": ["A", "B", "C"]})
    cens = pd.DataFrame({"chrom": ["chrI"], "start": [20000], "end": [20118], "centromere": ["CEN1"]})
    assert cis_genes(genes, cens, 10000)["orf"].tolist() == ["B"]
