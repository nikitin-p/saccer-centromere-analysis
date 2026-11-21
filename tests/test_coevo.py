import pandas as pd
import pytest

from centromere.coevo import (cell_label, centromere_variants, haplotype_by_centromere,
                              reference_elements, strain_matrix, substitutions)

REF = {"CENA": "GTCACATG" + "A" * 10 + "TGTTTATGATTTCCGAAAAAGAAAAA",
       "CENB": "ATCACGTG" + "T" * 12 + "TGATTTATGTTTCCGAAATAAAAATT"}


def test_substitutions():
    assert substitutions("ATCACATG", "GTCACATG") == "G1A"
    assert substitutions("ACGT", "ACGT") == ""
    with pytest.raises(ValueError):
        substitutions("ACG", "ACGT")


def test_reference_elements_ccg():
    r = reference_elements(REF).set_index("centromere")
    assert r.loc["CENA", "ref_cdei"] == "GTCACATG" and r.loc["CENA", "ref_cdeii_length"] == 10
    assert r.loc["CENA", "ccg_start"] == 13 and r.loc["CENB", "ccg_start"] == 13


def _inputs():
    ref = reference_elements(REF)
    a3, b3 = REF["CENA"][-26:], REF["CENB"][-26:]
    pre = pd.DataFrame({
        "strain": ["s1", "s1", "s2", "s2"], "centromere": ["CENA", "CENB"] * 2,
        "cdei_seq": ["ATCACATG", "ATCACGTG", "GTCACATG", "NNNNNNNN"],
        "cdeiii_seq": [a3[:13] + "A" + a3[14:], b3, a3, b3],   # s1 CENA: C14A inside CCG
        "cde_ok": [True, True, True, False]})
    steps = pd.DataFrame({"strain": ["s1", "s1", "s2", "s2"], "centromere": ["CENA", "CENB"] * 2,
                          "cdeii_length": [10, 14, 10, 12], "delta_s288c": [0, 2, 0, 0],
                          "class_s288c": ["within", "off_step", "within", "within"]})
    genes = pd.DataFrame({"strain": ["s1", "s2"], "gene": ["CBF1", "CBF1"],
                          "haplotype": ["H1", "ref"], "status": ["intact", "intact"]})
    return centromere_variants(pre, steps, ref), genes


def test_centromere_variants_and_labels():
    cen, _ = _inputs()
    c = cen.set_index(["strain", "centromere"])
    assert c.loc[("s1", "CENA"), "cdei_variants"] == "G1A"
    assert c.loc[("s1", "CENA"), "cdeiii_variants"] == "C14A" and c.loc[("s1", "CENA"), "ccg_hit"]
    assert c.loc[("s1", "CENB"), "cdei_core"] == "CACGTG" and not c.loc[("s1", "CENB"), "ccg_hit"]
    assert pd.isna(c.loc[("s2", "CENB"), "cdei_variants"])           # failed CDE call -> NA
    assert cell_label(c.loc[("s2", "CENB")], "cdei") == "?"
    assert cell_label(c.loc[("s1", "CENB")], "cdeii") == "+2"
    assert cell_label(c.loc[("s2", "CENA")], "cdeiii") == "."


def test_haplotype_table_and_strain_matrix():
    cen, genes = _inputs()
    t = haplotype_by_centromere(cen, genes, "CBF1", "cdei", ["CENA", "CENB"])
    assert list(t.index) == ["s2", "s1"]                                # ref first
    assert t.loc["s1", "CENA"] == "G1A" and t.loc["s2", "CENB"] == "?"
    m = strain_matrix(cen, genes)
    assert m.loc["s1", "n_cdei_alt"] == 1 and m.loc["s1", "n_ccg_hit"] == 1
    assert m.loc["s1", "n_cdeii_off_step"] == 1 and m.loc["s2", "n_no_cde_call"] == 1
    assert m.loc["s2", "n_cdei_alt"] == 0 and m.loc["s1", "n_cdei_cacatg"] == 1


def test_missing_cde_calls_and_loci():
    """No sequence (NaN) for a failed call, no CDEII step, and a strain without a locus row."""
    ref = reference_elements(REF)
    pre = pd.DataFrame({"strain": ["s1", "s1", "s2"], "centromere": ["CENA", "CENB", "CENA"],
                        "cdei_seq": ["ATCACATG", float("nan"), "GTCACATG"],
                        "cdeiii_seq": [REF["CENA"][-26:], float("nan"), REF["CENA"][-26:]],
                        "cde_ok": [True, False, True]})
    steps = pd.DataFrame({"strain": ["s1", "s2"], "centromere": ["CENA", "CENA"],
                          "cdeii_length": [10, 10], "delta_s288c": [0, 0],
                          "class_s288c": ["within", "within"]})
    genes = pd.DataFrame({"strain": ["s1", "s2"], "gene": ["CBF1", "CBF1"],
                          "haplotype": ["H1", "ref"], "status": ["intact", "intact"]})
    cen = centromere_variants(pre, steps, ref)
    c = cen.set_index(["strain", "centromere"])
    assert pd.isna(c.loc[("s1", "CENB"), "cdei_variants"]) and pd.isna(c.loc[("s1", "CENB"), "ccg_hit"])
    assert cell_label(c.loc[("s1", "CENB")], "cdeii") == "?"         # no CDEII length
    t = haplotype_by_centromere(cen, genes, "CBF1", "cdei", ["CENA", "CENB"])
    assert t.loc["s2", "CENB"] == "?"                                  # s2 has no CENB row
    m = strain_matrix(cen, genes)
    assert m.loc["s1", "n_no_cde_call"] == 1 and m.loc["s1", "n_cdei_alt"] == 1
