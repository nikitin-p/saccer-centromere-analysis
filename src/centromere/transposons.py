"""Analysis D: Ty elements and solo LTRs around centromeres.

Pure functions; BLAST is run by scripts/10_transposons.py. Coordinates are 1-based,
inclusive, on the contig's plus strand throughout.
"""
import numpy as np
import pandas as pd

HIT_COLS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send",
            "evalue", "bitscore", "qlen", "slen"]
MERGE_GAP = 10          # bp: HSPs closer than this belong to one element (LTR|INT|LTR abut)
FULL_INT_COV = 0.7      # internal-region coverage for a "full" element
SOLO_INT_COV = 0.1      # at most this much internal coverage for a "solo_LTR"


# ---------------------------------------------------------------- hits -> elements
def parse_hits(df: pd.DataFrame) -> pd.DataFrame:
    """Add kind (INT/LTR), family, strand, plus-strand s_lo/s_hi and query coverage."""
    df = df.copy()
    parts = df["qseqid"].astype(str).str.split("|")
    df["kind"] = parts.str[0]
    df["family"] = parts.str[1]   # NaN for queries without a family (tRNAs, flanks)
    df["strand"] = np.where(df["sstart"] <= df["send"], "+", "-")
    df["s_lo"] = df[["sstart", "send"]].min(axis=1)
    df["s_hi"] = df[["sstart", "send"]].max(axis=1)
    df["qcov"] = (df["qend"] - df["qstart"] + 1) / df["qlen"]
    return df


def merge_intervals(lo: np.ndarray, hi: np.ndarray, gap: int = 0) -> np.ndarray:
    """Cluster id per interval; intervals closer than `gap` bp share a cluster (input any order)."""
    order = np.argsort(lo, kind="stable")
    ids = np.empty(len(lo), dtype=int)
    cid, end = -1, -np.inf
    for i in order:
        if lo[i] > end + gap + 1:
            cid += 1
            end = hi[i]
        else:
            end = max(end, hi[i])
        ids[i] = cid
    return ids


def _union_len(intervals: list[tuple[int, int]]) -> int:
    total, cur_lo, cur_hi = 0, None, None
    for lo, hi in sorted(intervals):
        if cur_hi is None or lo > cur_hi + 1:
            if cur_hi is not None:
                total += cur_hi - cur_lo + 1
            cur_lo, cur_hi = lo, hi
        else:
            cur_hi = max(cur_hi, hi)
    return total + (cur_hi - cur_lo + 1 if cur_hi is not None else 0)


def call_elements(hits: pd.DataFrame, min_pident: float = 80, min_qcov: float = 0.7,
                  min_hsp: int = 50) -> tuple[pd.DataFrame, int]:
    """Merge HSPs into elements and classify them.

    An element is kept if it contains an LTR HSP covering >= min_qcov of its query, or
    its INT HSPs cover >= min_qcov of one family's internal region. Class by internal
    coverage: full (>= 0.7), solo_LTR (<= 0.1), truncated (between). Returns the
    elements and the number of dropped fragments.
    """
    h = hits[(hits["pident"] >= min_pident) & (hits["length"] >= min_hsp)]
    rows, dropped = [], 0
    for contig, g in h.groupby("sseqid", sort=False):
        g = g.reset_index(drop=True)
        g["cid"] = merge_intervals(g["s_lo"].to_numpy(), g["s_hi"].to_numpy(), MERGE_GAP)
        for _, e in g.groupby("cid"):
            ints = e[e["kind"] == "INT"]
            ltrs = e[(e["kind"] == "LTR") & (e["qcov"] >= min_qcov)]
            int_cov, int_fam = 0.0, None
            if not ints.empty:
                cov = {f: _union_len(list(zip(x["qstart"], x["qend"]))) / x["qlen"].iloc[0]
                       for f, x in ints.groupby("family")}
                int_fam = max(cov, key=cov.get)
                int_cov = cov[int_fam]
            if int_cov < min_qcov and ltrs.empty:
                dropped += 1
                continue
            if int_cov >= FULL_INT_COV:
                cls = "full"
            elif int_cov <= SOLO_INT_COV:
                cls = "solo_LTR"
            else:
                cls = "truncated"
            if cls == "solo_LTR":
                best = ltrs.loc[ltrs["bitscore"].idxmax()]
            else:
                pool = ints[ints["family"] == int_fam]
                best = pool.loc[pool["bitscore"].idxmax()]
            rows.append({"contig": contig, "start": int(e["s_lo"].min()), "end": int(e["s_hi"].max()),
                         "strand": best["strand"], "class": cls, "family": best["family"],
                         "int_cov": round(int_cov, 3), "n_ltr_hsps": len(ltrs),
                         "best_query": best["qseqid"], "pident": best["pident"]})
    cols = ["contig", "start", "end", "strand", "class", "family", "int_cov", "n_ltr_hsps",
            "best_query", "pident"]
    out = pd.DataFrame(rows, columns=cols).sort_values(["contig", "start"]).reset_index(drop=True)
    return out, dropped


def call_trnas(hits: pd.DataFrame, min_pident: float = 90, min_qcov: float = 0.9) -> pd.DataFrame:
    h = hits[(hits["pident"] >= min_pident) & (hits["qcov"] >= min_qcov)]
    rows = []
    for contig, g in h.groupby("sseqid", sort=False):
        g = g.reset_index(drop=True)
        g["cid"] = merge_intervals(g["s_lo"].to_numpy(), g["s_hi"].to_numpy(), 0)
        for _, e in g.groupby("cid"):
            best = e.loc[e["bitscore"].idxmax()]
            rows.append({"contig": contig, "start": int(e["s_lo"].min()), "end": int(e["s_hi"].max()),
                         "strand": best["strand"], "trna": best["qseqid"]})
    return pd.DataFrame(rows, columns=["contig", "start", "end", "strand", "trna"])


# ---------------------------------------------------------------- distances / windows
def interval_distance(lo1, hi1, lo2, hi2):
    """Gap in bp between intervals (0 if they overlap); works elementwise on arrays."""
    return np.maximum(0, np.maximum(np.asarray(lo2) - np.asarray(hi1), np.asarray(lo1) - np.asarray(hi2)) - 1)


def near_cen(elements: pd.DataFrame, cens: pd.DataFrame, max_dist: int = 20000) -> pd.DataFrame:
    """Elements within max_dist of each CEN (same genome, same contig).

    side: 'CDEI' if the element lies on the CDEI side of the CEN (in CEN orientation),
    'CDEIII' otherwise, 'overlap' if it overlaps the CEN. rel_strand: same/opposite.
    rel_lo/rel_hi: element ends in bp from the CEN centre, oriented so the CDEI side is negative.
    """
    rows = []
    for _, c in cens.iterrows():
        e = elements[(elements["genome"] == c["genome"]) & (elements["contig"] == c["contig"])]
        if e.empty:
            continue
        d = interval_distance(e["start"], e["end"], c["cen_start"], c["cen_end"])
        e = e[d <= max_dist].assign(distance=d[d <= max_dist])
        mid, sign = (c["cen_start"] + c["cen_end"]) / 2, 1 if c["strand"] == "+" else -1
        for _, x in e.iterrows():
            a, b = sign * (x["start"] - mid), sign * (x["end"] - mid)   # bp from CEN centre, CDEI side < 0
            if x["distance"] == 0:
                side = "overlap"
            else:
                left = x["end"] < c["cen_start"]
                side = "CDEI" if left == (c["strand"] == "+") else "CDEIII"
            rows.append({"genome": c["genome"], "centromere": c["centromere"], "contig": c["contig"],
                         "start": x["start"], "end": x["end"], "strand": x["strand"],
                         "class": x["class"], "family": x["family"], "distance": int(x["distance"]),
                         "side": side,
                         "rel_strand": "same" if x["strand"] == c["strand"] else "opposite",
                         "rel_lo": min(a, b), "rel_hi": max(a, b),
                         "element_id": x["element_id"]})
    return pd.DataFrame(rows)


def window_counts(elements: pd.DataFrame, contig: str, lo: np.ndarray, hi: np.ndarray,
                  windows: list[int]) -> np.ndarray:
    """Counts of elements within each window of each query interval -> (n_queries, n_windows)."""
    e = elements[elements["contig"] == contig]
    out = np.zeros((len(lo), len(windows)), dtype=int)
    if e.empty:
        return out
    d = interval_distance(e["start"].to_numpy()[None, :], e["end"].to_numpy()[None, :],
                          np.asarray(lo)[:, None], np.asarray(hi)[:, None])
    for j, w in enumerate(windows):
        out[:, j] = (d <= w).sum(axis=1)
    return out


def nearest_trna(trnas: pd.DataFrame, contig: str, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    t = trnas[trnas["contig"] == contig]
    if t.empty:
        return np.full(len(lo), np.inf)
    d = interval_distance(t["start"].to_numpy()[None, :], t["end"].to_numpy()[None, :],
                          np.asarray(lo)[:, None], np.asarray(hi)[:, None])
    return d.min(axis=1).astype(float)


# ---------------------------------------------------------------- null
def candidate_positions(contig_len: dict[str, int], cens: pd.DataFrame, n: int, length: int,
                        edge: int, rng: np.random.Generator) -> pd.DataFrame:
    """n random intervals of `length` bp, >= edge from contig ends and from every CEN."""
    names = [c for c, L in contig_len.items() if L > 2 * edge + length]
    lens = np.array([contig_len[c] for c in names], dtype=float)
    pick = rng.choice(len(names), size=n, p=lens / lens.sum())
    lo = (edge + 1 + rng.random(n) * (lens[pick] - 2 * edge - length)).astype(int)
    df = pd.DataFrame({"contig": np.array(names)[pick], "lo": lo, "hi": lo + length - 1})
    keep = np.ones(n, dtype=bool)
    for _, c in cens.iterrows():
        on = df["contig"].to_numpy() == c["contig"]
        d = interval_distance(df["lo"], df["hi"], c["cen_start"], c["cen_end"])
        keep &= ~(on & (d < edge))
    return df[keep].reset_index(drop=True)


def matched_pools(cand: pd.DataFrame, cen_d: np.ndarray, tol: float, per_cen: int,
                  rng: np.random.Generator, min_tol: float = 100) -> list[np.ndarray]:
    """For each CEN, up to per_cen candidate row indices with nearest-tRNA distance within
    max(tol * d, min_tol) of that CEN's distance (sampled without replacement)."""
    pools = []
    cd = cand["d_trna"].to_numpy()
    for d in cen_d:
        if np.isinf(d):
            idx = np.flatnonzero(np.isinf(cd))
        else:
            idx = np.flatnonzero(np.abs(cd - d) <= max(tol * d, min_tol))
        if len(idx) > per_cen:
            idx = rng.choice(idx, size=per_cen, replace=False)
        pools.append(idx)
    return pools


def replicate_sums(pool_counts: list[np.ndarray], n_rep: int, rng: np.random.Generator) -> np.ndarray:
    """Null distribution of the total count: one random draw from each CEN's pool, summed.
    pool_counts[i] is (n_pool_i, n_windows). Returns (n_rep, n_windows)."""
    total = np.zeros((n_rep, pool_counts[0].shape[1]))
    for pc in pool_counts:
        total += pc[rng.integers(0, len(pc), size=n_rep)]
    return total


def empirical_p(null: np.ndarray, obs: float) -> float:
    """Two-sided empirical p with +1 correction."""
    n = len(null)
    hi = (np.sum(null >= obs) + 1) / (n + 1)
    lo = (np.sum(null <= obs) + 1) / (n + 1)
    return float(min(1.0, 2 * min(hi, lo)))


# ---------------------------------------------------------------- orthology by junction flanks
def te_clusters(elements: pd.DataFrame, gap: int) -> pd.DataFrame:
    """Group a genome's elements into clusters (gap <= `gap` bp) so flanks are TE-free."""
    out = []
    for contig, g in elements.groupby("contig", sort=False):
        g = g.copy()
        g["cluster"] = merge_intervals(g["start"].to_numpy(), g["end"].to_numpy(), gap)
        for _, c in g.groupby("cluster"):
            out.append({"contig": contig, "start": int(c["start"].min()), "end": int(c["end"].max()),
                        "element_ids": ",".join(c["element_id"]),
                        "classes": ",".join(c["class"]), "families": ",".join(c["family"])})
    return pd.DataFrame(out)


def unique_best(hits: pd.DataFrame, min_pident: float, min_qcov: float, ratio: float = 0.9) -> pd.DataFrame:
    """Per query, the best hit if it passes the thresholds and the runner-up scores < ratio x best."""
    h = hits.sort_values("bitscore", ascending=False, kind="stable")
    rank = h.groupby("qseqid", sort=False).cumcount()
    best = h[rank == 0]
    second = h[rank == 1].set_index("qseqid")["bitscore"]
    runner = best["qseqid"].map(second)
    keep = ((best["pident"] >= min_pident) & (best["qcov"] >= min_qcov)
            & ~(runner >= ratio * best["bitscore"]))
    return best[keep.to_numpy()]


def flank_call(left: pd.Series | None, right: pd.Series | None, target_elements: pd.DataFrame,
               max_empty_gap: int = 50) -> tuple[str, int | None, str | None]:
    """Presence call at an orthologous site from the two junction-flank hits in a target genome.

    Flanks are oriented with the source contig: the left flank ends where the insertion
    starts. Calls: absent (flanks abut, |gap| <= max_empty_gap; a Ty target-site
    duplication shows as gap ~ -5), present (gap holds a called element), gap_no_te,
    rearranged (different contigs/strands or overlap beyond max_empty_gap), unresolved
    (a flank without a unique hit). Returns (call, gap, target element_ids)."""
    if left is None or right is None:
        return "unresolved", None, None
    if left["sseqid"] != right["sseqid"] or left["strand"] != right["strand"]:
        return "rearranged", None, None
    if left["strand"] == "+":
        gap_lo, gap_hi = left["s_hi"] + 1, right["s_lo"] - 1
    else:
        gap_lo, gap_hi = right["s_hi"] + 1, left["s_lo"] - 1
    gap = int(gap_hi - gap_lo + 1)
    if -max_empty_gap <= gap <= max_empty_gap:
        return "absent", gap, None
    if gap < -max_empty_gap:
        return "rearranged", gap, None
    t = target_elements[target_elements["contig"] == left["sseqid"]]
    ov = (np.minimum(t["end"], gap_hi) - np.maximum(t["start"], gap_lo) + 1).clip(lower=0)
    hit = t[ov >= 0.5 * (t["end"] - t["start"] + 1)]
    if hit.empty:   # filled by something that is not a called element (often an LTR fragment)
        return "gap_no_te", gap, None
    return "present", gap, ",".join(hit["element_id"])


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)
