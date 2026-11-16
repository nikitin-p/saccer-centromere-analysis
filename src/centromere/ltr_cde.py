"""Analysis E: CDE-like motifs (CDEI, CDEIII, full CDEI-CDEII-CDEIII units) in Ty LTRs.

A *unit* is a CDEI PWM hit, then a CDEII-like spacer (length and AT criteria), then a
CDEIII PWM hit, all on one strand -- the same layout `cde.CDEModel` annotates at real
centromeres, with the same thresholds (weakest S288C element - 2).

Coordinates: 0-based, half-open on the forward strand of the scanned sequence.
Spec: docs/specs/2026-11-13-ltr-cde-motifs.md
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from centromere.cde import CDEI_LEN, CDEIII_LEN, CDEModel, score_positions

_COMP = str.maketrans("ACGTN", "TGCAN")
EBOX = ("CACGTG", "CACATG", "CATGTG")   # CDEI cores (CACRTG, both strands)


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ---------------------------------------------------------------- scanning
@dataclass
class Hits:
    """Hits on both strands, forward-strand intervals [lo, hi)."""
    cdei: np.ndarray      # (n, 3): lo, hi, strand(+1/-1)
    cdeiii: np.ndarray
    units: pd.DataFrame   # strand, lo, hi, cdei_score, cdeiii_score, cdeii_len, cdeii_at


class CDEScanner:
    def __init__(self, model: CDEModel, cdeii_min: int = 70, cdeii_max: int = 95,
                 cdeii_min_at: float = 0.80):
        self.m = model
        self.cdeii_min, self.cdeii_max, self.cdeii_min_at = cdeii_min, cdeii_max, cdeii_min_at

    def _strand(self, seq: str):
        """Hits on the given strand, coordinates in that strand's own frame."""
        s1 = score_positions(seq, self.m.cdei)
        s3 = score_positions(seq, self.m.cdeiii)
        i1 = np.flatnonzero(s1 >= self.m.cdei_min)
        i3 = np.flatnonzero(s3 >= self.m.cdeiii_min)
        units = []
        if len(i1) and len(i3):
            at = np.concatenate([[0], np.cumsum([b in "AT" for b in seq])])
            for a in i1:
                lo, hi = a + CDEI_LEN + self.cdeii_min, a + CDEI_LEN + self.cdeii_max
                for b in i3[np.searchsorted(i3, lo):np.searchsorted(i3, hi, side="right")]:
                    L = b - a - CDEI_LEN
                    frac = (at[b] - at[a + CDEI_LEN]) / L
                    if frac >= self.cdeii_min_at:
                        units.append((a, b + CDEIII_LEN, s1[a], s3[b], L, frac))
        return i1, s1, i3, s3, units

    def hits(self, seq: str) -> Hits:
        seq = seq.upper()
        n = len(seq)
        c1, c3, us = [], [], []
        for strand, s in ((1, seq), (-1, revcomp(seq))):
            i1, _, i3, _, units = self._strand(s)
            if strand == 1:
                c1 += [(a, a + CDEI_LEN, 1) for a in i1]
                c3 += [(b, b + CDEIII_LEN, 1) for b in i3]
                us += [(1, a, e, x, y, L, f) for a, e, x, y, L, f in units]
            else:
                c1 += [(n - a - CDEI_LEN, n - a, -1) for a in i1]
                c3 += [(n - b - CDEIII_LEN, n - b, -1) for b in i3]
                us += [(-1, n - e, n - a, x, y, L, f) for a, e, x, y, L, f in units]
        cols = ["strand", "lo", "hi", "cdei_score", "cdeiii_score", "cdeii_len", "cdeii_at"]
        return Hits(np.array(c1, dtype=int).reshape(-1, 3), np.array(c3, dtype=int).reshape(-1, 3),
                    pd.DataFrame(us, columns=cols))

    def counts(self, seq: str) -> tuple[int, int, int]:
        h = self.hits(seq)
        return len(h.cdei), len(h.cdeiii), len(h.units)

    def summary(self, seq: str) -> dict:
        """Per-sequence row for ltr_cde_scores.tsv."""
        seq = seq.upper()
        rc = revcomp(seq)
        best1 = max(score_positions(seq, self.m.cdei).max(initial=-np.inf),
                    score_positions(rc, self.m.cdei).max(initial=-np.inf))
        best3 = max(score_positions(seq, self.m.cdeiii).max(initial=-np.inf),
                    score_positions(rc, self.m.cdeiii).max(initial=-np.inf))
        h = self.hits(seq)
        core = [(seq[lo:hi] if st == 1 else revcomp(seq[lo:hi]))[2:8] for lo, hi, st in h.cdei]
        row = {"length": len(seq), "at": round(sum(b in "AT" for b in seq) / max(len(seq), 1), 3),
               "best_cdei": round(float(best1), 2), "best_cdeiii": round(float(best3), 2),
               "n_cdei": len(h.cdei), "n_cdeiii": len(h.cdeiii), "n_units": len(h.units),
               "n_ebox": sum(seq.count(e) for e in EBOX),
               "cdei_core_cacrtg": sum(c in ("CACGTG", "CACATG") for c in core)}
        if len(h.units):
            u = h.units.assign(tot=h.units["cdei_score"] + h.units["cdeiii_score"])
            u = u.loc[u["tot"].idxmax()]
            unit = seq[int(u["lo"]):int(u["hi"])]
            if u["strand"] == -1:
                unit = revcomp(unit)
            row.update({"unit_strand": "+" if u["strand"] == 1 else "-", "unit_lo": int(u["lo"]),
                        "unit_hi": int(u["hi"]), "unit_cdei_score": round(u["cdei_score"], 2),
                        "unit_cdeiii_score": round(u["cdeiii_score"], 2),
                        "unit_cdeii_len": int(u["cdeii_len"]), "unit_cdeii_at": round(u["cdeii_at"], 3),
                        "unit_cdei": unit[:CDEI_LEN], "unit_cdeii": unit[CDEI_LEN:-CDEIII_LEN],
                        "unit_cdeiii": unit[-CDEIII_LEN:]})
        return row


# ---------------------------------------------------------------- dinucleotide shuffle
def dinuc_shuffle(seq: str, rng: np.random.Generator) -> str:
    """Altschul-Erickson shuffle: same dinucleotide counts, same first and last base."""
    if len(seq) < 3:
        return seq
    edges: dict[str, list[str]] = {}
    for a, b in zip(seq, seq[1:]):
        edges.setdefault(a, []).append(b)
    last = seq[-1]
    # Wilson's algorithm: random arborescence into `last` from loop-erased random walks;
    # final[v] is v's last-used out-edge in the Eulerian path
    in_tree, final = {last}, {}
    for v in edges:
        u = v
        while u not in in_tree:
            final[u] = edges[u][rng.integers(len(edges[u]))]
            u = final[u]
        u = v
        while u not in in_tree:
            in_tree.add(u)
            u = final[u]
    order: dict[str, list[str]] = {}
    for v, outs in edges.items():
        outs = list(outs)
        if v != last:
            outs.remove(final[v])
        rng.shuffle(outs)
        if v != last:
            outs.append(final[v])
        order[v] = outs
    out, cur, pos = [seq[0]], seq[0], {v: 0 for v in order}
    for _ in range(len(seq) - 1):
        nxt = order[cur][pos[cur]]
        pos[cur] += 1
        out.append(nxt)
        cur = nxt
    return "".join(out)


def dinuc_counts(seq: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for a, b in zip(seq, seq[1:]):
        out[a + b] = out.get(a + b, 0) + 1
    return out


# ---------------------------------------------------------------- LTR extraction
def ltr_intervals(elements: pd.DataFrame, hits: pd.DataFrame, ltr_len: dict[str, int],
                  min_pident: float = 80, min_qcov: float = 0.7, end_tol: int = 50) -> pd.DataFrame:
    """One row per LTR (1-based inclusive contig coords).

    solo_LTR -> the element itself. full -> two LTRs at the element ends: the best
    LTR HSP starting (ending) within `end_tol` bp of the element start (end), else a
    terminal segment of the family's median S288C LTR length.
    """
    h = hits[(hits["kind"] == "LTR") & (hits["pident"] >= min_pident) & (hits["qcov"] >= min_qcov)]
    by_contig = {c: g for c, g in h.groupby("sseqid")}
    rows = []
    for _, e in elements.iterrows():
        base = {"genome": e["genome"], "element_id": e["element_id"], "contig": e["contig"],
                "strand": e["strand"], "family": e["family"], "class": e["class"]}
        if e["class"] == "solo_LTR":
            rows.append({**base, "ltr": "solo", "start": e["start"], "end": e["end"], "source": "element"})
            continue
        if e["class"] != "full":
            continue
        g = by_contig.get(e["contig"])
        L = ltr_len[e["family"]]
        for side in ("left", "right"):
            src, lo, hi = "median_length", None, None
            if g is not None:
                if side == "left":
                    c = g[(abs(g["s_lo"] - e["start"]) <= end_tol) & (g["s_hi"] <= e["end"])]
                else:
                    c = g[(abs(g["s_hi"] - e["end"]) <= end_tol) & (g["s_lo"] >= e["start"])]
                if not c.empty:
                    best = c.loc[c["bitscore"].idxmax()]
                    lo, hi, src = int(best["s_lo"]), int(best["s_hi"]), "hsp"
            if lo is None:
                lo, hi = ((e["start"], e["start"] + L - 1) if side == "left"
                          else (e["end"] - L + 1, e["end"]))
            # 5'/3' relative to the element's own orientation
            name = ("5p" if side == "left" else "3p") if e["strand"] == "+" else \
                   ("3p" if side == "left" else "5p")
            rows.append({**base, "ltr": name, "start": lo, "end": hi, "source": src})
    return pd.DataFrame(rows)


def extract(contig_seq: str, start: int, end: int, strand: str) -> str:
    s = contig_seq[start - 1:end].upper()
    return revcomp(s) if strand == "-" else s


# ---------------------------------------------------------------- genomic windows
def merge(iv: np.ndarray) -> np.ndarray:
    """Merge [lo, hi) intervals (n, 2) -> sorted disjoint."""
    if len(iv) == 0:
        return iv.reshape(0, 2)
    iv = iv[np.argsort(iv[:, 0])]
    out = [list(iv[0])]
    for lo, hi in iv[1:]:
        if lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return np.array(out)


def free_mask(lo: np.ndarray, hi: np.ndarray, blocked: np.ndarray) -> np.ndarray:
    """True where [lo, hi) overlaps none of the merged blocked intervals."""
    if len(blocked) == 0:
        return np.ones(len(lo), bool)
    k = np.searchsorted(blocked[:, 0], hi, side="left")   # blocks starting before hi
    prev_hi = np.where(k > 0, blocked[np.maximum(k - 1, 0), 1], -1)
    return prev_hi <= lo


def sample_windows(contig_len: dict[str, int], blocked: dict[str, np.ndarray], lengths: np.ndarray,
                   n_rep: int, rng: np.random.Generator, max_rounds: int = 200,
                   allow_missing: bool = False) -> pd.DataFrame:
    """n_rep x len(lengths) windows, uniform over the genome, avoiding blocked intervals.

    Windows still unplaced after max_rounds raise, or with allow_missing are dropped
    (attrs["n_missing"] says how many) -- long windows may not fit gene-free gaps.
    """
    names = list(contig_len)
    sizes = np.array([contig_len[c] for c in names], float)
    lengths = np.tile(np.asarray(lengths, int), n_rep)
    rep = np.repeat(np.arange(n_rep), len(lengths) // n_rep)
    out_c = np.full(len(lengths), -1)
    out_lo = np.zeros(len(lengths), int)
    todo = np.arange(len(lengths))
    for _ in range(max_rounds):
        if len(todo) == 0:
            break
        L = lengths[todo]
        ci = rng.choice(len(names), size=len(todo), p=sizes / sizes.sum())
        lo = (rng.random(len(todo)) * np.maximum(sizes[ci] - L, 1)).astype(int)
        ok = lo + L <= sizes[ci]
        for k in np.unique(ci):
            sel = ci == k
            ok[sel] &= free_mask(lo[sel], lo[sel] + L[sel], blocked.get(names[k], np.empty((0, 2))))
        out_c[todo[ok]], out_lo[todo[ok]] = ci[ok], lo[ok]
        todo = todo[~ok]
    if len(todo) and not allow_missing:
        raise RuntimeError(f"{len(todo)} windows could not be placed")
    keep = out_c >= 0
    w = pd.DataFrame({"rep": rep[keep], "contig": [names[i] for i in out_c[keep]], "lo": out_lo[keep],
                      "hi": out_lo[keep] + lengths[keep]})
    w.attrs["n_missing"] = int(len(todo))
    return w


def count_contained(hit_lo: np.ndarray, hit_hi: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Per window, number of hits with lo <= hit_lo and hit_hi <= hi (hits sorted by lo)."""
    if len(hit_lo) == 0:
        return np.zeros(len(lo), int)
    a = np.searchsorted(hit_lo, lo, side="left")
    b = np.searchsorted(hit_lo, hi, side="left")
    out = np.zeros(len(lo), int)
    for i in np.flatnonzero(b > a):   # few windows contain hits; check their ends
        out[i] = int((hit_hi[a[i]:b[i]] <= hi[i]).sum())
    return out


def window_counts(hits: dict[str, Hits], windows: pd.DataFrame) -> pd.DataFrame:
    """Add n_cdei / n_cdeiii / n_units to windows from per-contig genome-wide hits."""
    w = windows.copy()
    for col in ("n_cdei", "n_cdeiii", "n_units"):
        w[col] = 0
    for contig, idx in w.groupby("contig").groups.items():
        h = hits[contig]
        lo, hi = w.loc[idx, "lo"].to_numpy(), w.loc[idx, "hi"].to_numpy()
        for col, arr in (("n_cdei", h.cdei[:, :2]), ("n_cdeiii", h.cdeiii[:, :2]),
                         ("n_units", h.units[["lo", "hi"]].to_numpy(int))):
            arr = arr[np.argsort(arr[:, 0], kind="stable")] if len(arr) else arr.reshape(0, 2)
            w.loc[idx, col] = count_contained(arr[:, 0], arr[:, 1], lo, hi)
    return w


# ---------------------------------------------------------------- statistics
def empirical_p(null: np.ndarray, obs: float) -> float:
    """Two-sided empirical p with +1 correction."""
    null = np.asarray(null, float)
    lo = (np.sum(null <= obs) + 1) / (len(null) + 1)
    hi = (np.sum(null >= obs) + 1) / (len(null) + 1)
    return float(min(1.0, 2 * min(lo, hi)))


def compare(obs: dict[str, float], obs_bp: float, null: pd.DataFrame, null_bp: np.ndarray) -> list[dict]:
    """Observed counts vs replicate sums (null: rows = replicates, cols = n_*)."""
    rows = []
    for col in ("n_cdei", "n_cdeiii", "n_units"):
        rate_null = null[col].to_numpy() / null_bp * 1000
        rate_obs = obs[col] / obs_bp * 1000
        rows.append({"motif": col[2:], "obs": int(obs[col]), "obs_per_kb": round(rate_obs, 4),
                     "null_mean_per_kb": round(float(rate_null.mean()), 4),
                     "null_q025": round(float(np.quantile(rate_null, 0.025)), 4),
                     "null_q975": round(float(np.quantile(rate_null, 0.975)), 4),
                     "ratio": round(rate_obs / rate_null.mean(), 3) if rate_null.mean() > 0 else np.nan,
                     "p": empirical_p(rate_null, rate_obs)})
    return rows
