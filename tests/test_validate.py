import pandas as pd

from centromere.validate import (compare_other_assembly, compare_pcan, compare_same_assembly,
                                 interval_shifts, locate, read_pcan, revcomp, seq_differences)

CEN = "ATCACGTG" + "AT" * 40 + "TGTTTTTGTTTTCCGAAAAAGAAAAA"   # 8 + 80 + 26 = 114 bp
LEFT, RIGHT = "GGGCCCGGGC" * 3, "CCGGGCCCGG" * 3


def _locus(strand, contig):
    if strand == "+":
        start = contig.find(CEN) + 1
    else:
        start = contig.find(revcomp(CEN)) + 1
    return {"sample": "GCA_1.1", "strain": "XYZ", "centromere": "CEN1", "sseqid": "c1",
            "strand": strand, "cen_start": start, "cen_end": start + len(CEN) - 1,
            "sequence": CEN, "cdeii_length": 80}


def test_locate_both_strands():
    contig = LEFT + CEN + RIGHT
    assert locate(CEN, contig) == [(31, 30 + len(CEN), "+")]
    assert locate(CEN, revcomp(contig))[0][2] == "-"


def test_interval_shifts_are_in_cen_orientation():
    ours = {"strand": "-", "cen_start": 100, "cen_end": 213}
    # other call starts 2 bp later at CDEI (= lower genomic end on minus) and ends 1 bp later
    assert interval_shifts(ours, 99, 211) == (2, 1)
    plus = {"strand": "+", "cen_start": 100, "cen_end": 213}
    assert interval_shifts(plus, 102, 214, extra_3p=1) == (2, 0)


def test_same_assembly_handles_helsen_extra_base_on_minus_strand():
    contig = LEFT + revcomp(CEN + "T") + RIGHT  # Helsen call = CEN + 1 bp
    loci = pd.DataFrame([_locus("-", contig)])
    helsen = pd.DataFrame([{"assembly": "GCA_1.2", "centromere": "CEN1", "variant": "A",
                            "sequence": CEN + "T"}])
    res = compare_same_assembly(loci, helsen, {"GCA_1.1": {"c1": contig}})
    row = res.iloc[0]
    assert row["status"] == "identical" and row["d5"] == 0 and row["d3"] == 0
    assert row["cdeii_diff"] == 0


def test_seq_differences_ignores_end_gaps():
    assert seq_differences(CEN, CEN[3:]) == {"snv": 0, "indel_bp": 0}
    mutated = CEN[:40] + ("C" if CEN[40] != "C" else "G") + CEN[41:]
    assert seq_differences(CEN, mutated) == {"snv": 1, "indel_bp": 0}
    duplicated = CEN[:50] + "ATATATATAT" + CEN[50:]
    assert seq_differences(CEN, duplicated)["indel_bp"] == 10


def test_other_assembly_keeps_closest_variant():
    loci = pd.DataFrame([_locus("+", LEFT + CEN + RIGHT)])
    far = CEN[:50] + "GGGGG" + CEN[55:]
    helsen = pd.DataFrame([
        {"assembly": "XYZ", "centromere": "CEN1", "variant": "A", "sequence": far + "A"},
        {"assembly": "XYZ", "centromere": "CEN1", "variant": "B", "sequence": CEN + "A"},
    ])
    row = compare_other_assembly(loci, helsen).iloc[0]
    assert row["variant"] == "B" and row["status"] == "identical" and row["n_variants"] == 2


def test_pcan_match_and_unmatched(tmp_path):
    contig = LEFT + CEN + RIGHT + LEFT * 5
    loci = pd.DataFrame([_locus("+", contig)])
    path = tmp_path / "CENsequences.txt"
    header = "Contig\tContigHit\tStart\tEnd\tSequence\tCDEIIlen\tCDEIIAT\tFIMOIscore\tFIMOIIscore\tOverallScore"
    s = loci.iloc[0]
    line1 = f"c1\tc1_1\t{s.cen_start}\t{s.cen_end - 8}\t{CEN[:-8]}\t80\t95.0\t20\t10\t4000"
    line2 = "c1\tc1_2\t200\t300\tACGT\t80\t90.0\t15\t9\t100"
    path.write_text("\n".join([header, line1, line2, "", "****************", "PCAn v1.0"]) + "\n")
    pc = read_pcan(path)
    assert len(pc) == 2 and pc["Start"].dtype.kind == "i"
    res = compare_pcan(loci, {"GCA_1.1": pc}, {"GCA_1.1": {"c1": contig}})
    found = res[res["status"] == "found"].iloc[0]
    assert (found["d5"], found["d3"], found["pcan_strand"]) == (0, -8, "+")
    assert (res["status"] == "pcan_only").sum() == 1
