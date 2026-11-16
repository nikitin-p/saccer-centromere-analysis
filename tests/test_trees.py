import numpy as np
import pandas as pd
import pytest

from centromere.trees import (concat_cen, discordance, dist_frame, encode, k80, swap_blocks, swap_check)


def test_k80_hand_computed():
    d, L = k80(encode(["AAAAAAAAAA", "GAAAAAAAAC"]))
    # 1 transition (A/G), 1 transversion (A/C) in 10 sites: P = Q = 0.1
    assert L[0, 1] == 10
    assert d[0, 1] == pytest.approx(-0.5 * np.log(0.7) - 0.25 * np.log(0.8))
    # gaps are pairwise-deleted; extra invariant sites enlarge L
    d2, L2 = k80(encode(["AA-AAAAAAAA", "GAAAAAAAAAC"]), extra_sites=0)
    assert L2[0, 1] == 10 and d2[0, 1] == pytest.approx(d[0, 1])
    d3, L3 = k80(encode(["GA", "AA"]), extra_sites=8)   # 1 transition in 10
    assert L3[0, 1] == 10 and d3[0, 1] == pytest.approx(-0.5 * np.log(0.8) - 0.25 * np.log(1.0))


def test_concat_cen_blocks_and_missing():
    aligned = pd.DataFrame({"centromere": ["CEN1", "CEN1", "CEN2"], "strain": ["X", "Y", "X"],
                            "ref_aligned_seq": ["AC-G", "AC-G", "TT"], "aligned_seq": ["ACTG", "AC-G", "TA"]})
    seqs, blocks = concat_cen(aligned, ["CEN1", "CEN2"])
    assert seqs == {"S288C": "AC-GTT", "X": "ACTGTA", "Y": "AC-G--"}
    assert blocks.to_dict("records") == [{"centromere": "CEN1", "start": 0, "end": 4},
                                         {"centromere": "CEN2", "start": 4, "end": 6}]


def _mutate(seq, rate, rng):
    s = seq.copy()
    m = rng.random(len(s)) < rate
    s[m] = (s[m] + rng.integers(1, 4, m.sum())) % 4
    return s


def _caterpillar(n_taxa, length, rng):
    """t0..t{n-1}: each lineage branches off a backbone that keeps evolving, so every
    taxon's distances to the others have a distinct order."""
    backbone = rng.integers(0, 4, length)
    out = []
    for k in range(n_taxa):
        backbone = _mutate(backbone, 0.03, rng)
        out.append(_mutate(backbone, 0.01, rng))
    return ["".join("ACGT"[x] for x in s) for s in out]


def test_swap_detected():
    rng = np.random.default_rng(3)
    taxa = [f"t{k}" for k in range(8)]
    full = _caterpillar(8, 16 * 150 + 20000, rng)
    cen = {t: s[:16 * 150] for t, s in zip(taxa, full)}
    blocks = pd.DataFrame({"centromere": [f"CEN{k + 1}" for k in range(16)],
                           "start": np.arange(16) * 150, "end": np.arange(1, 17) * 150})
    genome = dist_frame(k80(encode([s[16 * 150:] for s in full]))[0], taxa)

    base = discordance(dist_frame(k80(encode(list(cen.values())))[0], taxa), genome)
    assert base.max() < 0.2     # CENs share the genome's history

    r = swap_check(cen, blocks, genome, strain="t2", donor="t7", n=4)
    assert r["detected"], r
    assert r["most_discordant"] == "t2"
    assert set(r["top_contributors"]) == set(r["swapped"])
    sw = swap_blocks(cen, blocks, "t2", "t7", r["swapped"])
    b = blocks.set_index("centromere").loc[r["swapped"][0]]
    assert sw["t2"][b["start"]:b["end"]] == cen["t7"][b["start"]:b["end"]]
