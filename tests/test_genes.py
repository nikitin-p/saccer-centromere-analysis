import pandas as pd

from centromere.genes import (assign_haplotypes, extract_cds, frameshift_sites, longest_run,
                              pick_locus, protein_variants, to_plus, variant_label)

# 10-codon ORF: M K L V D E R S T W *
ORF = "ATGAAACTGGTGGATGAACGTTCTACCTGG" + "TAA"
UP, DOWN = "CCCTAGCC", "GGGCCC"  # upstream in-frame stop (TAG at offset 3 -> frame differs; fine)


def hsp(qstart, qend, o_start, o_end, qlen=10, contig="c1", strand="+", bits=50.0, slen=1000):
    sstart, send = (o_start + 1, o_end) if strand == "+" else (slen - o_start, slen - o_end + 1)
    return {"qseqid": "G|X", "sseqid": contig, "qstart": qstart, "qend": qend, "sstart": sstart,
            "send": send, "bitscore": bits, "qlen": qlen, "slen": slen}


def test_extract_intact_from_partial_hsp():
    seq = UP + ORF + DOWN
    o = len(UP)
    locus = pd.DataFrame([hsp(3, 10, o + 6, o + 30)])
    locus = pick_locus(locus)[0]
    cds = extract_cds(seq, locus)
    assert cds["protein"] == "MKLVDERSTW"
    assert cds["status"] == "intact" and cds["start_shift"] == 0
    assert cds["o_start"] == o and cds["o_end"] == o + len(ORF)


def test_downstream_start_when_projected_start_is_lost():
    # start codon ATG -> ACG; codon 4 GTG -> ATG becomes the start
    seq = UP + "ACGAAACTGATGGATGAACGTTCTACCTGGTAA" + DOWN
    o = len(UP)
    locus = pick_locus(pd.DataFrame([hsp(1, 10, o, o + 30)]))[0]
    cds = extract_cds(seq, locus)
    assert cds["protein"] == "MDERSTW" and cds["start_shift"] == 3
    v = protein_variants("MKLVDERSTW", cds["protein"])
    assert [variant_label(x) for x in v] == ["Ntrunc1-3", "V4M"]


def test_premature_stop():
    seq = UP + "ATGAAACTGTAG" + "GCT" * 26 + "TAA" + DOWN  # 30 codons, stop at codon 4
    o = len(UP)
    locus = pick_locus(pd.DataFrame([hsp(1, 3, o, o + 9, qlen=30),
                                     hsp(5, 30, o + 12, o + 90, qlen=30, bits=30)]))[0]
    cds = extract_cds(seq, locus)
    assert cds["protein"] == "MKL" and cds["status"] == "premature_stop"


def test_frameshift_flag_and_site():
    # 1-bp insertion (extra A after codon 4, inside an A-run) shifts the frame
    seq = UP + "ATGAAACTGGTGAAAAAGGATGAACGTTCTACCTGGTAA" + DOWN
    o = len(UP)
    locus = pick_locus(pd.DataFrame([hsp(1, 5, o, o + 15), hsp(6, 10, o + 16, o + 31)]))[0]
    cds = extract_cds(seq, locus)
    assert cds["status"] == "frameshift" and cds["frames"] == 2
    fs = frameshift_sites(seq, locus)
    assert len(fs) == 1 and fs[0]["homopolymer"] >= 5


def test_pick_locus_prefers_highest_score_and_counts_copies():
    h = pd.DataFrame([hsp(1, 10, 100, 130, contig="a", bits=40), hsp(1, 10, 500, 530, contig="b", bits=50),
                      hsp(1, 5, 900, 915, contig="b", strand="-", bits=10)])
    locus, n = pick_locus(h)
    assert set(locus["sseqid"]) == {"b"} and locus["o_start"].iloc[0] == 500
    assert n == 2  # contig a scores >= 50 % of the best


def test_minus_strand_coordinates():
    h = pd.DataFrame([hsp(1, 10, 100, 130, strand="-", slen=1000)])
    locus, _ = pick_locus(h)
    assert locus["sstart"].iloc[0] == 900 and locus["send"].iloc[0] == 871
    assert to_plus(100, 130, "-", 1000) == (871, 900)
    assert to_plus(100, 130, "+", 1000) == (101, 130)


def test_protein_variants_labels():
    ref = "MKLVDERSTWPQ"
    assert [variant_label(v) for v in protein_variants(ref, "MKLVDQRSTWPQ")] == ["E6Q"]
    assert [variant_label(v) for v in protein_variants(ref, "MKLVDERSTW")] == ["Ctrunc11-12"]
    v = protein_variants("MKLVDERSTWPQGHIKLM", "MKLVDERSTGHIKLM")
    assert [x["type"] for x in v] == ["deletion"] and len(v[0]["ref"]) == 3


def test_haplotypes_ref_and_frequency_order():
    g = pd.DataFrame({"gene": "G", "strain": list("abcde"), "protein": ["AAA", "AAC", "AAD", "AAD", "AAA"]})
    lab = assign_haplotypes(g, {"G": "AAA"})
    assert list(lab) == ["ref", "H2", "H1", "H1", "ref"]


def test_longest_run():
    assert longest_run("ACGTTTTTAC") == 5 and longest_run("") == 0


# ---------------------------------------------------------------- C.3
def test_msa_positions_ref_coordinates_insertions_and_private():
    from centromere.genes import msa_positions
    aln = {"REF": "MK-LVD", "a": "MKALVE", "b": "MK-LVE", "c": "-K-LVE"}
    p = msa_positions(aln, "REF").set_index("col")
    assert list(p.index) == [1, 3, 6]
    assert p.loc[1, "ref_pos"] == 1 and p.loc[1, "n_alt"] == 1 and p.loc[1, "n_sub"] == 0
    assert p.loc[3, "insertion"] and p.loc[3, "ref_pos"] == 2 and p.loc[3, "alt_strains"] == "a"
    assert p.loc[6, "ref_pos"] == 5 and p.loc[6, "ref_private"] and p.loc[6, "alleles"] == "E:3"
    assert not p.loc[3, "ref_private"]


def test_domains_and_summary():
    from centromere.genes import domain_summary, domains_at, msa_positions
    dom = pd.DataFrame({"domain": ["tail", "core", "core"], "start": [1, 3, 6], "end": [4, 4, 6]})
    assert domains_at(3, dom) == "tail;core" and domains_at(5, dom) == "-"
    # pos1 gap in a; pos3 L>M in b only; pos4 V>A and pos6 E>Q in both (ref-private); pos5 same
    aln = {"REF": "MKLVDE", "a": "-KLADQ", "b": "MKMADQ"}
    s = domain_summary(msa_positions(aln, "REF"), dom, 6).set_index("domain")
    assert s.loc["tail", "length"] == 4 and s.loc["tail", "n_variable"] == 3
    assert s.loc["tail", "n_substituted"] == 2 and s.loc["tail", "n_substituted_excl_ref_private"] == 1
    assert s.loc["tail", "subst_per_100aa"] == 25.0
    assert s.loc["core", "length"] == 3 and s.loc["core", "n_substituted"] == 3
    assert s.loc["core", "n_substituted_excl_ref_private"] == 1
    assert s.loc["unannotated", "length"] == 1 and s.loc["unannotated", "n_variable"] == 0
