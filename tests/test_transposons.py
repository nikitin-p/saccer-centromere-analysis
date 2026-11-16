import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio import SeqIO

from centromere.transposons import (HIT_COLS, UnionFind, call_elements, candidate_positions, empirical_p,
                                    flank_call, interval_distance, matched_pools, merge_intervals,
                                    near_cen, parse_hits, replicate_sums, unique_best, window_counts)

ROOT = Path(__file__).resolve().parents[1]
BLAST = Path("~/miniconda3/envs/blast/bin").expanduser()


def hsp(q, s, qs, qe, ss, se, qlen, pident=95.0, bits=None):
    return {"qseqid": q, "sseqid": "c1", "pident": pident, "length": abs(se - ss) + 1, "qstart": qs,
            "qend": qe, "sstart": ss, "send": se, "evalue": 0.0, "bitscore": bits or abs(se - ss),
            "qlen": qlen, "slen": 100000}


def test_merge_intervals_gap():
    ids = merge_intervals(np.array([1, 105, 500]), np.array([100, 200, 600]), gap=10)
    assert ids[0] == ids[1] != ids[2]


def test_interval_distance():
    assert interval_distance(1, 10, 11, 20) == 0      # abutting
    assert interval_distance(1, 10, 5, 20) == 0       # overlapping
    assert interval_distance(1, 10, 21, 30) == 10
    assert interval_distance(21, 30, 1, 10) == 10


def test_call_elements_full_solo_fragment():
    rows = [  # full Ty1: LTR | INT | LTR, abutting
        hsp("LTR|Ty1|a", 1000, 1, 330, 1000, 1329, 330),
        hsp("INT|Ty1|x", 1000, 1, 5000, 1330, 6329, 5250),
        hsp("LTR|Ty1|b", 1000, 1, 330, 6330, 6659, 330),
        # solo Ty3 LTR on the minus strand
        hsp("LTR|Ty3|s", 2000, 1, 340, 20340, 20001, 340),
        # fragment: 100 bp of an LTR (< 70 %) -> dropped
        hsp("LTR|Ty1|f", 3000, 1, 100, 40000, 40099, 330),
        # truncated: LTR + 30 % of the internal region
        hsp("LTR|Ty2|t", 4000, 1, 330, 60000, 60329, 330),
        hsp("INT|Ty2|y", 4000, 1, 1600, 60330, 61929, 5295),
    ]
    el, dropped = call_elements(parse_hits(pd.DataFrame(rows)))
    assert dropped == 1
    assert el[["class", "family"]].values.tolist() == [["full", "Ty1"], ["solo_LTR", "Ty3"],
                                                       ["truncated", "Ty2"]]
    assert el.loc[0, ["start", "end"]].tolist() == [1000, 6659]
    assert el.loc[1, "strand"] == "-"


def test_low_identity_hits_ignored():
    el, _ = call_elements(parse_hits(pd.DataFrame([hsp("LTR|Ty1|a", 0, 1, 330, 1, 330, 330, pident=75)])))
    assert el.empty


def elements_df(rows):
    return pd.DataFrame(rows, columns=["genome", "contig", "start", "end", "strand", "class", "family",
                                       "element_id"])


def test_near_cen_side_and_windows():
    el = elements_df([["g", "c1", 900, 950, "+", "solo_LTR", "Ty1", "e0"],     # 49 bp left of CEN
                      ["g", "c1", 1500, 1800, "-", "full", "Ty1", "e1"],       # 300 bp right
                      ["g", "c1", 30000, 30300, "+", "solo_LTR", "Ty3", "e2"]])
    cens = pd.DataFrame([{"genome": "g", "centromere": "CEN1", "contig": "c1", "cen_start": 1000,
                          "cen_end": 1199, "strand": "-"}])
    near = near_cen(el, cens, 20000)
    assert near["element_id"].tolist() == ["e0", "e1"]
    # CEN on minus strand: CDEI is at the right (high-coordinate) end
    assert near["side"].tolist() == ["CDEIII", "CDEI"]
    assert near["rel_strand"].tolist() == ["opposite", "same"]
    assert near.loc[0, ["rel_lo", "rel_hi"]].tolist() == [99.5 - 950 + 1000, 1099.5 - 900]  # right of centre
    assert near.loc[1, "rel_hi"] < 0
    assert window_counts(el, "c1", [1000], [1199], [10, 100, 1000, 50000])[0].tolist() == [0, 1, 2, 3]


def test_matched_pools_within_tolerance():
    rng = np.random.default_rng(0)
    cand = pd.DataFrame({"d_trna": rng.uniform(0, 50000, 20000)})
    pools = matched_pools(cand, np.array([3000.0, 20000.0]), 0.1, 100, rng)
    for d, p in zip([3000, 20000], pools):
        assert len(p) == 100
        assert np.all(np.abs(cand["d_trna"].to_numpy()[p] - d) <= 0.1 * d)


def test_candidate_positions_respect_edges_and_cens():
    rng = np.random.default_rng(1)
    cens = pd.DataFrame([{"contig": "a", "cen_start": 50000, "cen_end": 50120}])
    cand = candidate_positions({"a": 100000, "tiny": 1000}, cens, 5000, 120, 20000, rng)
    assert set(cand["contig"]) == {"a"}
    assert cand["lo"].min() > 20000 and cand["hi"].max() <= 80000
    assert (interval_distance(cand["lo"], cand["hi"], 50000, 50120) >= 20000).all()


def test_replicate_sums_and_p():
    rng = np.random.default_rng(2)
    pools = [np.ones((10, 1)), np.zeros((10, 1))]
    assert np.all(replicate_sums(pools, 100, rng) == 1)
    null = np.arange(100)
    assert empirical_p(null, 1000) == pytest.approx(2 / 101)
    assert empirical_p(null, 50) == 1.0


def flank(contig, lo, hi, strand):
    return pd.Series({"sseqid": contig, "s_lo": lo, "s_hi": hi, "strand": strand})


def test_flank_call_states():
    tgt = pd.DataFrame([{"contig": "t", "start": 1301, "end": 1634, "element_id": "T1"}])
    # empty site with a 5-bp target-site duplication: flanks overlap by 5 bp
    assert flank_call(flank("t", 1001, 1300, "+"), flank("t", 1296, 1595, "+"), tgt)[:2] == ("absent", -5)
    # insertion present between the flanks
    assert flank_call(flank("t", 1001, 1300, "+"), flank("t", 1635, 1934, "+"), tgt) == ("present", 334, "T1")
    # minus-strand orientation: left flank maps to higher coordinates
    assert flank_call(flank("t", 1635, 1934, "-"), flank("t", 1001, 1300, "-"), tgt)[0] == "present"
    assert flank_call(flank("t", 1, 300, "+"), flank("t", 501, 800, "+"), tgt)[0] == "gap_no_te"
    assert flank_call(flank("t", 1, 300, "+"), flank("u", 301, 600, "+"), tgt)[0] == "rearranged"
    assert flank_call(None, flank("t", 301, 600, "+"), tgt)[0] == "unresolved"


def test_unique_best_rejects_ambiguous():
    h = pd.DataFrame([{"qseqid": "q1", "bitscore": 500, "pident": 99, "qcov": 1.0},
                      {"qseqid": "q1", "bitscore": 480, "pident": 99, "qcov": 1.0},
                      {"qseqid": "q2", "bitscore": 500, "pident": 99, "qcov": 1.0},
                      {"qseqid": "q2", "bitscore": 100, "pident": 90, "qcov": 0.3}])
    assert unique_best(h, 90, 0.8)["qseqid"].tolist() == ["q2"]


def test_union_find():
    uf = UnionFind()
    uf.union("a", "b")
    uf.union("c", "b")
    assert uf.find("a") == uf.find("c") and uf.find("d") == "d"


# ---------------------------------------------------------------- BLAST integration (check 2)
@pytest.mark.skipif(not (BLAST / "blastn").exists(), reason="BLAST+ not installed")
def test_synthetic_insertions_blast(tmp_path):
    lib = {r.id: str(r.seq) for r in SeqIO.parse(ROOT / "data/ty_reference.fasta", "fasta")}
    int_ty1 = next(v for k, v in lib.items() if k.startswith("INT|Ty1"))
    ltr = next(v for k, v in lib.items() if k.startswith("LTR|Ty1"))
    rng = np.random.default_rng(3)
    bg = "".join(rng.choice(list("ACGT"), size=40000, p=[0.3, 0.2, 0.2, 0.3]))
    tsd1, tsd2 = bg[10000:10005], bg[25000:25005]
    with_te = bg[:10005] + ltr + int_ty1 + ltr + tsd1 + bg[10005:25005] + ltr + tsd2 + bg[25005:]
    for name, seq in (("src", with_te), ("empty", bg)):
        (tmp_path / f"{name}.fa").write_text(f">c1\n{seq}\n")
        subprocess.run([str(BLAST / "makeblastdb"), "-in", str(tmp_path / f"{name}.fa"), "-dbtype", "nucl"],
                       check=True, capture_output=True)

    def blast(query, db, task):
        res = subprocess.run([str(BLAST / "blastn"), "-task", task, "-query", str(query), "-db", str(db),
                              "-evalue", "1e-10", "-dust", "no", "-outfmt", "6 " + " ".join(HIT_COLS)],
                             check=True, capture_output=True, text=True)
        rows = [line.split("\t") for line in res.stdout.splitlines()]
        df = pd.DataFrame(rows, columns=HIT_COLS)
        for c in HIT_COLS[2:]:
            df[c] = pd.to_numeric(df[c])
        return parse_hits(df)

    el, _ = call_elements(blast(ROOT / "data/ty_reference.fasta", tmp_path / "src.fa", "blastn"))
    assert el[["class", "family"]].values.tolist() == [["full", "Ty1"], ["solo_LTR", "Ty1"]]
    full_len = 2 * len(ltr) + len(int_ty1)
    assert abs(el.loc[0, "start"] - 10006) <= 5 and abs(el.loc[0, "end"] - (10005 + full_len)) <= 5

    # junction flanks of the solo LTR: absent in the empty genome (TSD -> gap ~ -5), present in itself
    s, e = int(el.loc[1, "start"]), int(el.loc[1, "end"])
    (tmp_path / "fl.fa").write_text(f">L\n{with_te[s - 301:s - 1]}\n>R\n{with_te[e:e + 300]}\n")
    el = el.assign(element_id=["E0", "E1"], contig="c1")
    for db, expect in (("empty", "absent"), ("src", "present")):
        best = unique_best(blast(tmp_path / "fl.fa", tmp_path / f"{db}.fa", "megablast"), 90, 0.8)
        best = {r["qseqid"]: r for _, r in best.iterrows()}
        call, gap, tid = flank_call(best.get("L"), best.get("R"), el)
        assert call == expect
        if expect == "absent":
            assert -6 <= gap <= 0
        else:
            assert tid == "E1"
