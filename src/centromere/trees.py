"""Analysis H: centromere tree vs genome tree.

Both distance matrices come from the same K80 code: a CEN matrix from the concatenated
per-CEN MSAs (aligned.tsv) and a genome matrix from single-copy core SNPs in the minimap2
alignments to S288C. Trees, RF and Mantel are done in R (scripts/trees_nj.R).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr



def encode(seqs: list[str]) -> np.ndarray:
    """1-4 = ACGT, 0 = anything else (gap, N)."""
    lut = np.zeros(256, dtype=np.int8)
    for i, b in enumerate("ACGT", 1):
        lut[ord(b)] = lut[ord(b.lower())] = i
    return np.stack([lut[np.frombuffer(s.encode(), dtype=np.uint8)] for s in seqs])


def concat_cen(aligned: pd.DataFrame, order: list[str]) -> tuple[dict[str, str], pd.DataFrame]:
    """Concatenate the per-CEN MSAs (S288C = ref_aligned_seq). A strain without a locus gets
    gaps for that CEN. Returns {taxon: seq} and the block table (centromere, start, end)."""
    strains = sorted(aligned["strain"].unique())
    if "S288C" in strains:
        raise ValueError("strain S288C collides with the reference taxon; rename it in the sample map")
    seqs = {t: [] for t in ["S288C"] + strains}
    blocks, pos = [], 0
    for cen in order:
        a = aligned[aligned["centromere"] == cen]
        if a.empty:
            continue
        L = len(a["ref_aligned_seq"].iloc[0])
        by = dict(zip(a["strain"], a["aligned_seq"]))
        seqs["S288C"].append(a["ref_aligned_seq"].iloc[0])
        for s in strains:
            seqs[s].append(by.get(s, "-" * L).upper())
        blocks.append({"centromere": cen, "start": pos, "end": pos + L})
        pos += L
    return {t: "".join(v) for t, v in seqs.items()}, pd.DataFrame(blocks)


def k80(codes: np.ndarray, extra_sites: int = 0, chunk: int = 200_000) -> tuple[np.ndarray, np.ndarray]:
    """Kimura 2-parameter distance per pair (pairwise deletion of 0 codes).

    `extra_sites` = invariant sites comparable in every pair but not passed in `codes`
    (the genome matrix stores variable core sites only). Returns (distance, comparable sites);
    distance is NaN where the K80 log argument is not positive. All pairs at once as
    indicator-matrix products (column chunks), so 1000+ taxa x 1e5+ sites stay cheap.
    """
    n = len(codes)
    L = np.zeros((n, n))
    same = np.zeros((n, n))
    same_class = np.zeros((n, n))
    for lo in range(0, codes.shape[1], chunk):
        c = codes[:, lo:lo + chunk]
        ok = (c > 0).astype(np.float32)
        L += ok @ ok.T
        for b in (1, 2, 3, 4):
            x = (c == b).astype(np.float32)
            same += x @ x.T
        for cls in ((1, 3), (2, 4)):          # purines A/G, pyrimidines C/T
            x = np.isin(c, cls).astype(np.float32)
            same_class += x @ x.T
    ts = same_class - same
    tv = (L - same) - ts
    L = L + extra_sites
    with np.errstate(divide="ignore", invalid="ignore"):
        P, Q = ts / L, tv / L
        x, y = 1 - 2 * P - Q, 1 - 2 * Q
        d = np.where((x > 0) & (y > 0), -0.5 * np.log(np.where(x > 0, x, 1)) - 0.25 * np.log(np.where(y > 0, y, 1)),
                     np.nan)
    np.fill_diagonal(d, 0.0)
    np.fill_diagonal(L, 0.0)
    return d, L


def dist_frame(d: np.ndarray, taxa: list[str]) -> pd.DataFrame:
    return pd.DataFrame(d, index=taxa, columns=taxa)


def discordance(d1: pd.DataFrame, d2: pd.DataFrame) -> pd.Series:
    """1 - Spearman rho between each taxon's row in the two matrices (self excluded)."""
    taxa = list(d1.index)
    d2 = d2.loc[taxa, taxa]
    out = {}
    for t in taxa:
        others = [x for x in taxa if x != t]
        out[t] = 1 - spearmanr(d1.loc[t, others], d2.loc[t, others]).statistic
    return pd.Series(out, name="discordance")


def cen_matrix(seqs: dict[str, str], blocks: pd.DataFrame, drop: list[str] = ()) -> pd.DataFrame:
    """K80 matrix of the concatenation, optionally without some CEN blocks."""
    taxa = list(seqs)
    keep = np.ones(len(next(iter(seqs.values()))), dtype=bool)
    for _, b in blocks[blocks["centromere"].isin(drop)].iterrows():
        keep[b["start"]:b["end"]] = False
    codes = encode(list(seqs.values()))[:, keep]
    return dist_frame(k80(codes)[0], taxa)


def cen_contribution(seqs: dict[str, str], blocks: pd.DataFrame, genome: pd.DataFrame) -> pd.DataFrame:
    """Per strain x CEN: discordance(full) - discordance(without CEN k); positive = CEN k
    pushes the strain away from its genome-tree position."""
    full = discordance(cen_matrix(seqs, blocks), genome)
    rows = []
    for cen in blocks["centromere"]:
        wo = discordance(cen_matrix(seqs, blocks, [cen]), genome)
        for t in full.index:
            rows.append({"strain": t, "centromere": cen, "discordance": full[t],
                         "discordance_without": wo[t], "contribution": full[t] - wo[t]})
    return pd.DataFrame(rows)


def swap_blocks(seqs: dict[str, str], blocks: pd.DataFrame, strain: str, donor: str,
                cens: list[str]) -> dict[str, str]:
    """Copy of the concatenation where `strain` carries `donor`'s sequence at `cens`."""
    s = list(seqs[strain])
    for _, b in blocks[blocks["centromere"].isin(cens)].iterrows():
        s[b["start"]:b["end"]] = seqs[donor][b["start"]:b["end"]]
    return {**seqs, strain: "".join(s)}


def swap_targets(seqs: dict[str, str], blocks: pd.DataFrame, strain: str, donor: str, n: int) -> list[str]:
    """The n CENs where strain and donor differ most (so the swap changes something)."""
    codes = encode([seqs[strain], seqs[donor]])
    diff = (codes[0] > 0) & (codes[1] > 0) & (codes[0] != codes[1])
    score = {b["centromere"]: diff[b["start"]:b["end"]].sum() for _, b in blocks.iterrows()}
    return sorted(score, key=lambda c: (-score[c], c))[:n]


def swap_check(seqs: dict[str, str], blocks: pd.DataFrame, genome: pd.DataFrame,
               strain: str, donor: str, n: int) -> dict:
    """Spec check 8: after the swap, is `strain` the most discordant taxon and do the
    swapped CENs top its contribution list?"""
    cens = swap_targets(seqs, blocks, strain, donor, n)
    sw = swap_blocks(seqs, blocks, strain, donor, cens)
    disc = discordance(cen_matrix(sw, blocks), genome).sort_values(ascending=False)
    contrib = cen_contribution(sw, blocks, genome)
    top = (contrib[contrib["strain"] == strain].sort_values("contribution", ascending=False)
           ["centromere"].head(n).tolist())
    return {"strain": strain, "donor": donor, "swapped": cens, "most_discordant": disc.index[0],
            "discordance": disc.round(3).to_dict(), "top_contributors": top,
            "detected": bool(disc.index[0] == strain and set(top) == set(cens))}


def swap_pairs(taxa: list[str], max_pairs: int | None, seed: int = 1) -> list[tuple[str, str]]:
    """Every ordered (strain, donor) pair, or a seeded random sample of max_pairs of them."""
    pairs = [(x, y) for x in taxa for y in taxa if x != y]
    if max_pairs is not None and len(pairs) > max_pairs:
        rng = np.random.default_rng(seed)
        pairs = [pairs[i] for i in sorted(rng.choice(len(pairs), max_pairs, replace=False))]
    return pairs


def default_swap_pair(cen: pd.DataFrame, genome: pd.DataFrame) -> tuple[str, str]:
    """The a-priori easiest case: the most concordant strain (S288C excluded) and its most
    genome-distant donor. Used when the configured pair is not in the data."""
    disc = discordance(cen, genome).drop("S288C", errors="ignore")
    strain = disc.idxmin()
    donor = genome.loc[strain].drop([strain, "S288C"], errors="ignore").idxmax()
    return strain, donor


def swap_scan(seqs: dict[str, str], blocks: pd.DataFrame, genome: pd.DataFrame, n: int,
              max_pairs: int | None = None, seed: int = 1) -> pd.DataFrame:
    """swap_check for every ordered (strain, donor) pair (or max_pairs sampled ones): how often
    the spec-8 criteria hold on this data (power of the discordance / contribution read-out)."""
    rows = []
    for x, y in swap_pairs(list(seqs), max_pairs, seed):
        r = swap_check(seqs, blocks, genome, x, y, n)
        rank = sorted(r["discordance"], key=lambda t: -r["discordance"][t]).index(x) + 1
        rows.append({"strain": x, "donor": y, "genome_dist": genome.loc[x, y],
                     "swapped": ",".join(r["swapped"]), "discordance_rank": rank,
                     "most_discordant": rank == 1,
                     "n_swapped_in_top": len(set(r["top_contributors"]) & set(r["swapped"])),
                     "detected": r["detected"]})
    return pd.DataFrame(rows)
