import random
import shutil
from pathlib import Path

import pandas as pd
import pytest

from centromere.cde import CDEModel, read_reference
from centromere.hmm import (add_margin, best_hits, build_profile, classify, drifting,
                            locate_in_genome, model_lengths, new_members, parse_tblout,
                            project_cen, run_nhmmer, write_fasta)
from centromere.locate import compare_loci
from centromere.microevo import revcomp

REF = read_reference(Path(__file__).parent.parent / "data" / "sacer_centromere.fasta")

TBL = """\
# target name        accession  query name           accession  hmmfrom hmm to alifrom  ali to envfrom  env to  sq len strand   E-value  score  bias  description of target
c1                   -          CEN1                 -                1     118     101     218     101     218    1000    +     1e-30   95.0   0.5  some desc
c1                   -          CEN2                 -               10      40     150     120     150     120    1000    -     1e-2    12.0   0.1  some desc
#
"""


def hit(profile, contig, ef, et, bits, ev=1e-20, strand="+", hf=1, ht=100):
    return {"contig": contig, "profile": profile, "hmm_from": hf, "hmm_to": ht, "ali_from": ef,
            "ali_to": et, "env_from": ef, "env_to": et, "sq_len": 10**6, "strand": strand,
            "evalue": ev, "bits": bits}


def test_parse_tblout_and_model_lengths(tmp_path):
    p = tmp_path / "x.tbl"
    p.write_text(TBL)
    t = parse_tblout(p)
    assert list(t["profile"]) == ["CEN1", "CEN2"]
    assert t.loc[1, "env_from"] == 150 and t.loc[1, "strand"] == "-" and t.loc[0, "bits"] == 95.0
    h = tmp_path / "m.hmm"
    h.write_text("HMMER3/f\nNAME  CEN1\nLENG  118\n//\nHMMER3/f\nNAME  CEN2\nLENG  117\n//\n")
    assert model_lengths(h) == {"CEN1": 118, "CEN2": 117}


def test_margin_uses_overlapping_other_profiles_on_either_strand():
    h = pd.DataFrame([hit("CEN1", "c1", 100, 220, 95), hit("CEN2", "c1", 200, 150, 30, strand="-"),
                      hit("CEN2", "c1", 5000, 5100, 80), hit("CEN1", "c2", 1, 100, 50)])
    m = add_margin(h)
    assert m.loc[0, "runner_up"] == "CEN2" and m.loc[0, "margin"] == 65
    assert pd.isna(m.loc[2, "runner_up"]) and m.loc[2, "margin"] == 80
    assert m.loc[3, "runner_up_bits"] == 0      # other contig


def test_best_hit_skips_loci_won_by_another_profile():
    # CEN15's top hit sits on CEN14's locus (where CEN14 scores higher); its own locus is weaker
    h = pd.DataFrame([hit("CEN14", "c1", 100, 220, 100), hit("CEN15", "c1", 100, 220, 60),
                      hit("CEN15", "c2", 300, 420, 55), hit("CEN16", "c3", 1, 100, 40, ev=1e-5)])
    b = best_hits(h, ["CEN14", "CEN15", "CEN16"], evalue=1e-10, margin=10).set_index("profile")
    assert b.loc["CEN15", "contig"] == "c2" and b.loc["CEN15", "margin"] == 55
    assert b.loc["CEN14", "margin"] == 40 and b.loc["CEN14", "n_loci"] == 1
    assert b.loc["CEN16", "n_loci"] == 0 and pd.isna(b.loc["CEN16", "bits"])   # E too weak


def test_best_hit_falls_back_to_a_lost_locus():
    h = pd.DataFrame([hit("CEN14", "c1", 100, 220, 100), hit("CEN15", "c1", 100, 220, 95)])
    b = best_hits(h, ["CEN15"]).iloc[0]
    assert b["contig"] == "c1" and b["margin"] == -5 and b["n_loci"] == 0


def test_project_cen_extends_partial_hits():
    plus = pd.Series(hit("CEN1", "c", 110, 200, 50, hf=11, ht=91))
    assert project_cen(plus, 100) == (100, 209, True)
    minus = pd.Series(hit("CEN1", "c", 200, 110, 50, strand="-", hf=11, ht=91))
    assert project_cen(minus, 100) == (101, 210, False)


def test_classify():
    base = {"margin": 50, "cde_ok": True, "length": 118}
    assert classify(base, 118) == "ok"
    assert classify({**base, "margin": 9}, 118) == "runner_up"
    assert classify({**base, "cde_ok": False}, 118) == "cde"
    assert classify({**base, "length": 300}, 118) == "length"
    assert classify({**base, "length": 300}, 118, anchored=True) == "insertion"


def test_new_members_dedupes_and_respects_frozen():
    loci = pd.DataFrame({"centromere": ["CEN1", "CEN1", "CEN1", "CEN2", "CEN3"],
                         "status": ["ok", "ok", "ok", "ok", "cde"],
                         "sequence": ["AAA", "CCC", "CCC", "GGG", "TTT"]})
    n = new_members(loci, {"CEN1": {"AAA"}, "CEN2": set()}, frozen={"CEN2"})
    assert list(n["sequence"]) == ["CCC"]


def test_drifting():
    r = pd.DataFrame({"round": [0, 0, 1, 1], "centromere": ["CEN1", "CEN2"] * 2,
                      "rejected_runner_up": [0, 2, 1, 2]})
    assert drifting(r, 0) == set() and drifting(r, 1) == {"CEN1"}


def test_compare_loci_ok_mask_is_positional():
    df = pd.DataFrame({"sample": ["s", "s"], "centromere": ["CEN1", "CEN2"], "contig": ["c", "c"],
                       "l4_start": [10, 500], "l4_end": [20, 600]}, index=[7, 3])
    other = pd.DataFrame({"sample": ["s", "s"], "centromere": ["CEN1", "CEN2"], "contig": ["c", "c"],
                          "start": [12, 700], "end": [22, 800]})
    m = compare_loci(df, other, "l1", "l4", ok=pd.Series([True, True], index=[7, 3]))
    assert list(m["agree_l1"]) == [True, False] and list(m["d_start_l1"]) == [-2, -200]
    m = compare_loci(df, other, "l1", "l4", ok=pd.Series([False, True], index=[7, 3]))
    assert pd.isna(m.loc[0, "agree_l1"]) and m.loc[1, "agree_l1"] is False


@pytest.mark.skipif(not (shutil.which("hmmbuild") and shutil.which("nhmmer") and shutil.which("mafft")),
                    reason="needs HMMER and MAFFT")
def test_hmm_finds_cens_on_both_strands(tmp_path):
    rng = random.Random(1)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    cens = ["CEN3", "CEN5"]
    mut = REF["CEN5"][:40] + ("T" if REF["CEN5"][40] != "T" else "A") + REF["CEN5"][41:]
    paths = [build_profile("CEN3", [REF["CEN3"]], tmp_path),
             build_profile("CEN5", [REF["CEN5"], mut], tmp_path)]
    allh = tmp_path / "all.hmm"
    allh.write_text("".join(p.read_text() for p in paths))
    left3, left5 = rand(3000), rand(5000)
    contigs = {"ctgA": left3 + REF["CEN3"] + rand(3000),
               "ctgB": left5 + revcomp(mut) + rand(4000)}
    genome = tmp_path / "g.fa"
    write_fasta(genome, list(contigs.items()))
    hits = run_nhmmer(allh, genome, tmp_path / "g.tbl", cpu=1)
    loci = locate_in_genome(hits, contigs, model_lengths(allh), CDEModel(REF),
                            {c: len(REF[c]) for c in cens}).set_index("centromere")
    assert (loci["status"] == "ok").all()
    assert loci.loc["CEN3", ["contig", "hmm_start", "hmm_end", "strand"]].tolist() == \
        ["ctgA", 3001, 3000 + len(REF["CEN3"]), "+"]
    assert loci.loc["CEN5", ["contig", "hmm_start", "strand"]].tolist() == ["ctgB", 5001, "-"]
    assert loci.loc["CEN5", "sequence"] == mut


def test_scan_cen_like_skips_known_and_runner_up_spots():
    from centromere.hmm import scan_cen_like
    rng = random.Random(1)
    bg = "".join(rng.choice("ACGT") for _ in range(3000))
    cen = REF["CEN3"]
    contig = bg[:500] + cen + bg[500:1500] + cen + bg[1500:]
    p1, p2 = 501, 500 + len(cen) + 1000 + 1          # 1-based starts of the two copies
    L = len(cen)
    hits = pd.DataFrame([hit("CEN3", "c", p1, p1 + L - 1, 90, hf=1, ht=L),
                         hit("CEN3", "c", p2, p2 + L - 1, 90, hf=1, ht=L),
                         hit("CEN5", "c", p2, p2 + L - 1, 85, hf=1, ht=L),   # margin 5 at copy 2
                         hit("CEN7", "c", 2800, 2900, 40, ev=1e-3)])         # too weak
    known = pd.DataFrame({"contig": ["c"], "start": [p1], "end": [p1 + L - 1]})
    model = CDEModel(REF)
    lengths = {"CEN3": L, "CEN5": L, "CEN7": 119}
    assert scan_cen_like(hits, {"c": contig}, lengths, model, known).empty
    s = scan_cen_like(hits, {"c": contig}, lengths, model, known, margin=3)
    assert len(s) == 1 and s.iloc[0]["start"] == p2 and s.iloc[0]["end"] == p2 + L - 1
    assert bool(s.iloc[0]["cde_ok"]) and s.iloc[0]["sequence"] == cen and s.iloc[0]["hmm_cov"] == 1.0
