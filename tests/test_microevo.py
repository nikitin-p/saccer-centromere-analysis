import numpy as np
import pandas as pd
import pytest

from centromere.microevo import (call_steps, describe_extra, helsen_cdeii, indel_events,
                                 reference_alleles, step_class)

# A CDEII-like, non-repetitive AT-rich sequence (S288C CEN1 CDEII)
CDEII = "ACGTAGTAATAAATAATTTTAAAAATATTTTATATTTTAAATAGTTTTCAAATATTTTACAGTTTATTTTTTAAATTTATTTAT"
CDEI, CDEIII = "ATCACGTG", "TGTTTTTGTTTTCCGAAGCAGTCAAA"


@pytest.mark.parametrize("delta,cls", [(0, "within"), (-2, "within"), (3, "off_step"), (-5, "off_step"),
                                       (8, "one_turn"), (-10, "one_turn"), (12, "one_turn"),
                                       (15, "off_step"), (20, "two_turns"), (-22, "two_turns"),
                                       (23, "large")])
def test_step_class(delta, cls):
    assert step_class(delta) == cls


def test_ten_bp_tandem_duplication_is_one_turn():
    unit = CDEII[30:40]
    alt = CDEII[:40] + unit + CDEII[40:]
    loci = pd.DataFrame({"strain": ["s1", "s2", "s3"], "centromere": "CEN1",
                         "cdeii_seq": [CDEII, CDEII, alt]})
    steps, events, _ = call_steps(loci, {"CEN1": CDEI + CDEII + CDEIII}, ["strain"])
    row = steps[steps["strain"] == "s3"].iloc[0]
    assert row["delta_median"] == 10 and row["class_median"] == "one_turn"
    ev = events[events["strain"] == "s3"]
    assert len(ev) == 1
    e = ev.iloc[0]
    assert e["type"] == "insertion" and e["length"] == 10
    assert e["tandem"] == "exact" and e["orientation"] == "direct" and e["unit_length"] == 10


def test_tandem_deletion_and_microhomology():
    # ref has XUUY; deleting one U is a tandem deletion; junction microhomology
    # comes from the base shared by the end of X and the end of U
    ref = "GCGCA" + "TTCAGA" * 2 + "CGGCC"
    alt = "GCGCA" + "TTCAGA" + "CGGCC"
    (e,) = indel_events(ref, alt)
    assert e["type"] == "deletion" and e["length"] == -6
    assert e["tandem"] == "exact" and e["microhomology"] == 1  # 'A' before TTCAGA = its last base
    assert e["ref_pos"] == 4  # left-normalised


def test_non_tandem_insertion_and_inverted_copy():
    seq = "GGGGG" + "ACCTGA" + "CCCCC"
    d = describe_extra(seq, 5, 6)
    assert d["tandem"] == "none"
    inv = "GGGG" + "ACCTGA" + "TCAGGT" + "CCCC"  # extra TCAGGT = revcomp(ACCTGA)
    d = describe_extra(inv, 10, 6)
    assert d["tandem"] == "inverted" and d["orientation"] == "inverted"


def test_multi_copy_expansion_unit():
    ref = "GGCGC" + "ATTAT" * 2 + "CCGCG"
    alt = "GGCGC" + "ATTAT" * 4 + "CCGCG"
    (e,) = indel_events(ref, alt)
    assert e["tandem"] == "exact" and e["unit_length"] == 5 and e["copies_added"] == 2


def test_reference_allele_is_modal_of_median_length():
    loci = pd.DataFrame({"centromere": "CEN1",
                         "cdeii_seq": ["AAAA", "AATA", "AATA", "AAAAAAAA", "AA"]})
    ref = reference_alleles(loci).iloc[0]
    assert ref["median_length"] == 4 and ref["ref_seq"] == "AATA"


def test_helsen_cdeii_strips_cde_and_extra_base():
    seq = CDEI + CDEII + CDEIII + "T"
    h = pd.DataFrame({"assembly": ["a", "b"], "centromere": "CEN1", "variant": "A",
                      "sequence": [seq, seq]})
    out = helsen_cdeii(h)
    assert out["cdeii_seq"].tolist() == [CDEII] and out["n_assemblies"].tolist() == [2]


def test_imperfect_tandem_with_one_base_shift():
    # 12-bp unit duplicated with one extra base: insertion is 13 bp, not an exact repeat
    unit = "TTTATTAATTAA"
    ref = "GCGCGACGCG" + unit + "GGCCGG"
    alt = "GCGCGACGCG" + unit + "C" + unit + "GGCCGG"
    (e,) = [x for x in indel_events(ref, alt) if abs(x["length"]) >= 3]
    assert e["length"] == 13 and e["tandem"] == "imperfect" and e["orientation"] == "direct"


def test_divergent_flag_uses_gap_free_identity():
    unit = CDEII[30:52]
    alt = CDEII[:52] + unit + CDEII[52:]
    loci = pd.DataFrame({"strain": ["a", "b", "c"], "centromere": "CEN1",
                         "cdeii_seq": [CDEII, CDEII, alt]})
    steps, _, _ = call_steps(loci, {"CEN1": CDEI + CDEII + CDEIII}, ["strain"])
    row = steps[steps["strain"] == "c"].iloc[0]
    assert row["identity_to_ref"] == 1.0 and not row["divergent"]


# --- Analysis F: pi (spec check 6) ---
from centromere.microevo import diversity_windows, pairwise_pi, region_map


def test_pairwise_pi_exact():
    seqs = ["ACGTACGTAC", "ACGTTCGTAC", "AC--ACGTAC"]
    r = pairwise_pi(seqs, np.array(["g"] * 10)).iloc[0]
    # pairs: (0,1) 1 sub / 10; (0,2) 0 sub, 1 indel / 8; (1,2) 1 sub, 1 indel / 8
    assert r["n_pairs"] == 3
    assert r["pi_sub"] == pytest.approx((0.1 + 0 + 0.125) / 3)
    assert r["indel_rate"] == pytest.approx((0 + 0.125 + 0.125) / 3)
    assert r["pi"] == pytest.approx(r["pi_sub"] + r["indel_rate"])


def test_pairwise_pi_terminal_gaps_are_missing():
    # a truncated sequence is not an indel; adjacent gaps in different sequences are 2 events
    r = pairwise_pi(["ACGTACGT", "--GTACGT"], np.array(["g"] * 8)).iloc[0]
    assert r["pi"] == 0 and r["indel_rate"] == 0
    r = pairwise_pi(["AC-GTA", "ACG-TA"], np.array(["g"] * 6)).iloc[0]
    assert r["pi_sub"] == 0 and r["indel_rate"] == pytest.approx(2 / 4)


def test_pairwise_pi_na_when_too_few_pairs():
    seqs = ["ACGTACGT", "ACGTACGT", "ACGTACGT"]
    g = np.array(["a"] * 4 + ["b"] * 4)
    seqs[1] = "ACGT----"
    seqs[2] = "ACGT----"            # group b: two sequences missing -> 0 of 3 pairs usable
    r = pairwise_pi(seqs, g).set_index("group")
    assert r.loc["a", "n_pairs"] == 3 and r.loc["a", "pi"] == 0
    assert r.loc["b", "n_pairs"] == 0 and np.isnan(r.loc["b", "pi"])


def test_region_map_and_windows():
    # S288C window: 3 bp left flank, CEN of 8+4+26 = 38 bp, 3 bp right flank; window = 2 bp
    cen = "TCACGTGA" + "AAAA" + "C" * 26
    s = "GGG" + cen + "TTT"
    s_aln = s[:5] + "--" + s[5:]                 # S288C gap inside CDEI
    rm = region_map(s_aln, cen_start=3, cen_len=38, window=2)
    assert list(rm["region"][:3]) == ["left_2", "left_1", "left_1"]
    assert list(rm["region"][5:7]) == ["CDEI", "CDEI"]   # gap columns take the preceding position
    assert (rm["region"] == "CDEII").sum() == 4
    assert list(rm["region"][-3:]) == ["right_1", "right_1", "right_2"]
    other = s[:5] + "AA" + s[5:13] + "G" + s[14:]  # 2-bp insertion in CDEI, 1 sub in CDEII
    d = diversity_windows({"S288C": s_aln, "x": other}, "S288C", 3, 38, window=2).set_index("region")
    assert d.loc["CDEI", "indel_rate"] == pytest.approx(1 / 8)
    assert d.loc["CDEII", "pi_sub"] == pytest.approx(1 / 4)
    assert d.loc["left_1", "offset"] == -1 and d.loc["right_2", "offset"] == 2
    assert list(d.index[:2]) == ["left_2", "left_1"]
