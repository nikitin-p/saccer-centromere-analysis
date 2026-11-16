import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from centromere.locate import Block, lift_interval, lift_window, liftover_paf, parse_cs, parse_paf, ref_alleles
from centromere.microevo import revcomp

# target 0..19, block covers target 2..18 and query 10..26 (+) :
#   :3  t2-4  -> q10-12
#   *ag t5    -> q13      (substitution, query base g)
#   +tt       -> q14-15   (insertion in query)
#   :4  t6-9  -> q16-19
#   -cc t10-11 deleted
#   :6  t12-17 -> q20-25
CS = ":3*ag+tt:4-cc:6"
PLUS = Block("ctg", 40, 10, 26, "+", "chr", 2, 18, 60, CS)
MINUS = Block("ctg", 40, 10, 26, "-", "chr", 2, 18, 60, CS)


def test_parse_cs_spans():
    o = parse_cs(CS)
    assert list(o["op"]) == ["=", "*", "+", "=", "-", "="]
    assert o["tlen"].sum() == 16 and o["qlen"].sum() == 16
    with pytest.raises(ValueError):
        parse_cs(":3x")


def test_lift_plus():
    assert liftover_paf([PLUS], "chr", 2) == ("ctg", 10, "+")
    assert liftover_paf([PLUS], "chr", 5) == ("ctg", 13, "+")      # substitution
    assert liftover_paf([PLUS], "chr", 6) == ("ctg", 16, "+")      # after the insertion
    assert liftover_paf([PLUS], "chr", 12) == ("ctg", 20, "+")     # after the deletion
    assert liftover_paf([PLUS], "chr", 17) == ("ctg", 25, "+")
    assert liftover_paf([PLUS], "chr", 10) == "in_insertion_gap"
    assert liftover_paf([PLUS], "chr", 18) == "unaligned"
    assert liftover_paf([PLUS], "chr", 1) == "unaligned"
    assert liftover_paf([PLUS], "chrX", 5) == "unaligned"


def test_lift_minus():
    # k-th aligned query base is qe-1-k = 25-k
    assert liftover_paf([MINUS], "chr", 2) == ("ctg", 25, "-")
    assert liftover_paf([MINUS], "chr", 6) == ("ctg", 19, "-")
    assert liftover_paf([MINUS], "chr", 17) == ("ctg", 10, "-")
    assert lift_interval([MINUS], "chr", 2, 18) == ("ctg", 10, 26, "-")


def test_lift_interval_reasons():
    assert lift_interval([PLUS], "chr", 3, 9) == ("ctg", 11, 19, "+")
    assert lift_interval([PLUS], "chr", 3, 11) == "in_insertion_gap"
    other = Block("ctg2", 40, 0, 4, "+", "chr", 18, 22, 60, ":4")
    assert lift_interval([PLUS, other], "chr", 3, 20) == "split_contigs"
    flip = Block("ctg", 40, 30, 34, "-", "chr", 18, 22, 60, ":4")
    assert lift_interval([PLUS, flip], "chr", 3, 20) == "strand_flip"
    dup = Block("ctg3", 40, 0, 4, "+", "chr", 4, 8, 60, ":4")
    assert liftover_paf([PLUS, dup], "chr", 5) == "multiple"


def test_lift_window_shortens_to_block():
    w = lift_window([PLUS], "chr", 8, 10, flank=10, step=1)
    assert w["status"] == "shortened"
    # left: target 2 is the first aligned base; right: 18 is unaligned, 17 lifts (the deletion
    # at 10-11 lies inside the span and is kept)
    assert (w["left_kept"], w["right_kept"]) == (6, 8)
    assert (w["qstart"], w["qend"]) == (10, 26)


def test_ref_alleles():
    ref = "AAACGTACGTCCAAAAAATT"
    a = ref_alleles([PLUS], "chr", ref)
    assert "".join("-ACGTN"[x] for x in a) == "--ACGGACGTNNAAAAAA--"
    dup = Block("ctg3", 40, 0, 4, "+", "chr", 4, 8, 60, ":4")
    a2 = ref_alleles([PLUS, dup], "chr", ref)
    assert (a2[4:8] == 0).all() and a2[3] == 2


MM2 = Path("~/miniconda3/envs/minimap2/bin/minimap2").expanduser()
S288C = Path("~/Desktop/pyselectal-analysis/data/sacCer3.fa").expanduser()


@pytest.mark.skipif(not (MM2.exists() and S288C.exists()), reason="minimap2 or sacCer3 missing")
def test_self_alignment_identity(tmp_path):
    from Bio import SeqIO
    chr3 = next(str(r.seq).upper() for r in SeqIO.parse(S288C, "fasta") if r.id == "chrIII")
    region = chr3[94000:135000]                      # CEN3 = chrIII:114384-114501
    (tmp_path / "t.fa").write_text(f">chrIII\n{region}\n")
    (tmp_path / "q.fa").write_text(f">fwd\n{region}\n>rev\n{revcomp(region)}\n")
    paf = tmp_path / "a.paf"
    with open(paf, "w") as fh:
        subprocess.run([str(MM2), "-x", "asm5", "-c", "--cs", str(tmp_path / "t.fa"), str(tmp_path / "q.fa")],
                       stdout=fh, stderr=subprocess.DEVNULL, check=True)
    blocks = parse_paf(paf)
    fwd = [b for b in blocks if b.qname == "fwd"]
    rev = [b for b in blocks if b.qname == "rev"]
    s, e = 114383 - 94000, 114501 - 94000
    assert lift_interval(fwd, "chrIII", s, e) == ("fwd", s, e, "+")
    n = len(region)
    assert lift_interval(rev, "chrIII", s, e) == ("rev", n - e, n - s, "-")
    a = ref_alleles(fwd, "chrIII", region)
    assert "".join("-ACGTN"[x] for x in a[s:e]) == region[s:e]


def test_bridge_collinear_gap():
    # target 0-9 -> query 100-109, gap target 10-14 / query 110-116, target 15-24 -> query 117-126
    L = Block("c", 300, 100, 110, "+", "chr", 0, 10, 60, ":10")
    R = Block("c", 300, 117, 127, "+", "chr", 15, 25, 60, ":10")
    assert lift_interval([L, R], "chr", 11, 14) == "unaligned"
    assert lift_interval([L, R], "chr", 11, 14, bridge=50) == ("c", 110, 117, "+")
    assert lift_interval([L, R], "chr", 5, 20, bridge=50) == ("c", 105, 122, "+")
    assert lift_interval([L, R], "chr", 11, 14, bridge=3) == "unaligned"          # gap too long
    other = Block("d", 300, 117, 127, "+", "chr", 15, 25, 60, ":10")
    assert lift_interval([L, other], "chr", 11, 14, bridge=50) == "unaligned"     # not collinear
    # reverse strand: target 0-9 -> query 190-199 (descending), target 15-24 -> query 173-182
    Lm = Block("c", 300, 190, 200, "-", "chr", 0, 10, 60, ":10")
    Rm = Block("c", 300, 173, 183, "-", "chr", 15, 25, 60, ":10")
    assert lift_interval([Lm, Rm], "chr", 11, 14, bridge=50) == ("c", 183, 190, "-")
    w = lift_window([L, R], "chr", 11, 14, flank=8, step=1, bridge=50)
    assert w["status"] == "ok" and (w["qstart"], w["qend"]) == (103, 124)


# ---- L2: flank locator -------------------------------------------------------------------
import pandas as pd  # noqa: E402

from centromere.locate import locate_by_flanks  # noqa: E402


def _hit(side, sseqid, sstart, send, qstart=1, qend=100, pident=99.0, bitscore=180.0,
         cen="CEN1", sample="s1"):
    return {"sample": sample, "qseqid": f"{cen}|{side}", "centromere": cen, "side": side,
            "qstart": qstart, "qend": qend, "qlen": 100, "sseqid": sseqid,
            "sstart": sstart, "send": send, "pident": pident, "length": qend - qstart + 1,
            "evalue": 0.0, "bitscore": bitscore}


def _locate(*hits, **kw):
    return locate_by_flanks(pd.DataFrame(list(hits)), **kw).iloc[0]


def test_flanks_plus_strand():
    # contig: L flank 101-200, CEN 201-317, R flank 318-417
    r = _locate(_hit("L", "c", 101, 200), _hit("R", "c", 318, 417))
    assert (r["status"], r["contig"], r["l2_start"], r["l2_end"], r["strand"]) == ("ok", "c", 201, 317, "+")
    assert r["l2_length"] == 117


def test_flanks_minus_strand():
    # CEN in reverse orientation: R flank (rc) 101-200, CEN 201-317, L flank (rc) 318-417
    r = _locate(_hit("L", "c", 417, 318), _hit("R", "c", 200, 101))
    assert (r["status"], r["l2_start"], r["l2_end"], r["strand"]) == ("ok", 201, 317, "-")


def test_flanks_partial_hsp_extrapolated():
    # L HSP stops 10 bp short of the CEN (qend 90); R HSP starts 5 bp in (qstart 6)
    r = _locate(_hit("L", "c", 101, 190, qend=90), _hit("R", "c", 323, 417, qstart=6))
    assert (r["l2_start"], r["l2_end"], r["ext_left"], r["ext_right"]) == (201, 317, 10, 5)
    r = _locate(_hit("L", "c", 417, 328, qend=90), _hit("R", "c", 195, 101, qstart=6))
    assert (r["l2_start"], r["l2_end"], r["strand"]) == (201, 317, "-")


def test_flanks_insertion_and_too_far():
    r = _locate(_hit("L", "c", 101, 200), _hit("R", "c", 6318, 6417))
    assert (r["status"], r["l2_length"]) == ("ok", 6117)
    assert _locate(_hit("L", "c", 101, 200), _hit("R", "c", 6318, 6417), max_locus=5000)["status"] == "too_far"


def test_flanks_best_pair_not_best_single_hit():
    # strongest L hit is a paralog on another contig; the pair on contig c wins
    r = _locate(_hit("L", "p", 101, 200, bitscore=190), _hit("L", "c", 101, 200, bitscore=170),
                _hit("R", "c", 318, 417))
    assert (r["status"], r["contig"], r["n_left"]) == ("ok", "c", 2)


@pytest.mark.parametrize("hits,status", [
    ([_hit("R", "c", 318, 417)], "no_left"),
    ([_hit("L", "c", 101, 200)], "no_right"),
    ([_hit("L", "c", 101, 200, pident=80), _hit("R", "c", 318, 417, qend=40)], "no_flanks"),
    ([_hit("L", "a", 101, 200), _hit("R", "b", 318, 417)], "split_contigs"),
    ([_hit("L", "c", 101, 200), _hit("R", "c", 417, 318)], "strand_flip"),
    ([_hit("L", "c", 318, 417), _hit("R", "c", 101, 200)], "order"),
])
def test_flanks_failure_reasons(hits, status):
    assert _locate(*hits)["status"] == status


def test_flanks_split_hsp_refines_boundary():
    # 50-bp insertion in the L flank 20 bp from the CEN: main HSP q1-80 -> 101-180, then
    # insertion 181-230, then q81-100 -> 231-250; CEN 251-367; R flank 368-467
    hits = [_hit("L", "c", 101, 180, qend=80), _hit("L", "c", 231, 250, qstart=81, qend=100, bitscore=40),
            _hit("R", "c", 368, 467)]
    r = _locate(*hits, min_refine_len=20)
    assert (r["l2_start"], r["l2_end"], r["ext_left"], r["split_left"]) == (251, 367, 0, True)
    # same on the minus strand: R (rc) 101-200, CEN 201-317, L tail q81-100 -> 337-318,
    # insertion 338-387, L main q1-80 -> 467-388
    hits = [_hit("L", "c", 467, 388, qend=80), _hit("L", "c", 337, 318, qstart=81, qend=100, bitscore=40),
            _hit("R", "c", 200, 101)]
    r = _locate(*hits, min_refine_len=20)
    assert (r["l2_start"], r["l2_end"], r["strand"], r["split_left"]) == (201, 317, "-", True)
    # a short hit elsewhere on the contig (beyond the R flank) is not collinear -> ignored
    hits = [_hit("L", "c", 101, 180, qend=80), _hit("L", "c", 900, 919, qstart=81, qend=100),
            _hit("R", "c", 368, 467)]
    r = _locate(*hits, min_refine_len=20)
    assert (r["l2_start"], r["split_left"]) == (201, False)


# ---------------------------------------------------------------- phase 4: reconcile

from centromere.locate import apply_reconcile, reconcile  # noqa: E402

REF_LEN = {"CEN1": 118, "CEN2": 117}


def _loci(rows):
    """rows: (tag, sample, cen, contig, start, end[, exact])"""
    import pandas as pd
    out = {}
    for r in rows:
        tag, smp, cen, ctg, s, e, *ex = r
        out.setdefault(tag, []).append({"sample": smp, "centromere": cen, "contig": ctg, "start": s,
                                        "end": e, "strand": "+", "exact": ex[0] if ex else True})
    return {t: pd.DataFrame(v) for t, v in out.items()}


def _rec(rows, ran=("l1", "l2", "l3", "l4"), samples=("s",), **kw):
    loci = _loci(rows)
    return reconcile(loci, {t: set(samples) for t in ran}, REF_LEN, **kw).set_index("centromere")


def test_reconcile_agree_and_missing():
    r = _rec([(t, "s", "CEN1", "c1", 100, 217) for t in ("l1", "l2", "l3", "l4")])
    assert r.loc["CEN1", "category"] == "agree" and r.loc["CEN1", "support"] == "L1,L2,L3,L4"
    assert r.loc["CEN1", "source"] == "L1" and r.loc["CEN1", "start"] == 100
    assert r.loc["CEN2", "category"] == "missing"


def test_reconcile_rescued_prefers_l4_coordinates():
    r = _rec([("l2", "s", "CEN1", "c1", 98, 219), ("l3", "s", "CEN1", "c1", 95, 220),
              ("l4", "s", "CEN1", "c1", 100, 217)])
    assert r.loc["CEN1", "category"] == "rescued" and r.loc["CEN1", "source"] == "L4"
    assert (r.loc["CEN1", "start"], r.loc["CEN1", "end"]) == (100, 217)
    r = _rec([("l2", "s", "CEN1", "c1", 98, 219)])
    assert r.loc["CEN1", "category"] == "rescued" and r.loc["CEN1", "source"] == "L2"


def test_reconcile_sequence_only_partial_conflict():
    r = _rec([("l1", "s", "CEN1", "c1", 100, 217), ("l4", "s", "CEN1", "c1", 100, 217)])
    assert r.loc["CEN1", "category"] == "sequence_only"
    # L2/L3 never ran -> not sequence_only: everything that ran found it
    r = _rec([("l1", "s", "CEN1", "c1", 100, 217)], ran=("l1",))
    assert r.loc["CEN1", "category"] == "agree"
    r = _rec([("l1", "s", "CEN1", "c1", 100, 217), ("l2", "s", "CEN1", "c1", 100, 217)])
    assert r.loc["CEN1", "category"] == "partial"
    r = _rec([("l1", "s", "CEN1", "c1", 100, 217), ("l2", "s", "CEN1", "c2", 5, 122),
              ("l3", "s", "CEN1", "c2", 5, 122)])
    assert r.loc["CEN1", "category"] == "conflict" and r.loc["CEN1", "n_clusters"] == 2
    assert r.loc["CEN1", "source"] == "L1" and "L2,L3:c2:5-122" in r.loc["CEN1", "clusters"]


def test_reconcile_insertion_needs_exact_anchor():
    big = [("l1", "s", "CEN1", "c1", 100, 160), ("l2", "s", "CEN1", "c1", 100, 6217),
           ("l3", "s", "CEN1", "c1", 100, 6217)]
    r = _rec(big)
    assert r.loc["CEN1", "category"] == "insertion" and r.loc["CEN1", "insertion_bp"] == 6000
    assert r.loc["CEN1", "source"] == "L2"
    # a bridged (inexact) L3 interval 60 bp too long is the alignment gap, not an insertion
    r = _rec([("l1", "s", "CEN1", "c1", 100, 217), ("l2", "s", "CEN1", "c1", 100, 217),
              ("l3", "s", "CEN1", "c1", 80, 257, False), ("l4", "s", "CEN1", "c1", 100, 217)])
    assert r.loc["CEN1", "category"] == "agree"
    r = _rec([("l1", "s", "CEN1", "c1", 100, 217), ("l2", "s", "CEN1", "c1", 100, 247)])
    assert r.loc["CEN1", "category"] == "partial"     # +30 bp: not above the threshold


def test_apply_reconcile_replaces_insertion_and_adds_rescued():
    import pandas as pd
    df = pd.DataFrame({"sample": ["s", "s"], "centromere": ["CEN1", "CEN2"], "x": [1, 2]})
    rec = pd.DataFrame({"sample": ["s", "s"], "centromere": ["CEN1", "CEN2"],
                        "category": ["insertion", "agree"], "support": ["L1,L2", "L1,L2"]})
    resc = pd.DataFrame({"sample": ["s"], "centromere": ["CEN1"], "x": [9], "locator": ["L2"]})
    out = apply_reconcile(df, rec, resc).set_index("centromere")
    assert len(out) == 2 and out.loc["CEN1", "x"] == 9 and out.loc["CEN1", "locator"] == "L2"
    assert out.loc["CEN2", "locator"] == "L1" and out.loc["CEN1", "locus_category"] == "insertion"
