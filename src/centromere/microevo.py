"""Analysis B: CDEII length in helical steps and the indels that change it.

CDEII wraps once around the Cse4 nucleosome, so a length change of about one
DNA turn (~10 bp) keeps CDEI and CDEIII on the same face of the helix, while a
change of half a turn (~5 bp) rotates one relative to the other (Helsen 2026).

Steps are measured against the population: per centromere, the median CDEII
length and a reference allele (the most frequent CDEII sequence among those of
median length). S288C is reported alongside but is not the baseline, because its
own CEN4 CDEII (77 bp) is an 8-bp deletion relative to almost every other strain.

Indel events come from a global alignment of each CDEII to the reference allele.
Each event is shifted to its leftmost and rightmost equivalent position in the
sequence that carries the extra bases; the width of that range (`slide`) tells
whether the extra bases repeat their neighbour (tandem: slide >= length) and how
much junction microhomology there is (slide - length for tandem events, slide
otherwise).
"""
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner

from centromere.cde import CDEI_LEN, CDEIII_LEN

HELICAL_REPEAT = 10.2
CLASS_ORDER = ["within", "off_step", "one_turn", "two_turns", "large"]
_COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def step_class(delta: int) -> str:
    """Within +-2 bp; one turn 8-12; two turns 18-22; off-step 3-7 and 13-17; large > 22."""
    a = abs(int(delta))
    if a <= 2:
        return "within"
    if 8 <= a <= 12:
        return "one_turn"
    if 18 <= a <= 22:
        return "two_turns"
    if a <= 17:
        return "off_step"
    return "large"


# ------------------------------------------------------------ reference allele

def reference_alleles(cdeii: pd.DataFrame, weight: str | None = None) -> pd.DataFrame:
    """Per centromere: median CDEII length and the modal sequence of that length.

    `cdeii` has columns centromere, cdeii_seq (and optionally a count column `weight`).
    Ties between equally frequent sequences go to the alphabetically first.
    """
    rows = []
    cdeii = cdeii.reset_index(drop=True)
    for cen, g in cdeii.groupby("centromere"):
        w = g[weight] if weight else pd.Series(1, index=g.index)
        lengths = np.repeat(g["cdeii_seq"].str.len().to_numpy(), w.to_numpy())
        median = int(np.round(np.median(lengths)))
        same = g[g["cdeii_seq"].str.len() == median]
        if same.empty:  # even count with no allele at the rounded median
            same = g.iloc[[int(np.argmin(np.abs(g["cdeii_seq"].str.len() - median)))]]
        counts = (same.assign(_w=w.loc[same.index]).groupby("cdeii_seq")["_w"].sum()
                  .reset_index().sort_values(["_w", "cdeii_seq"], ascending=[False, True]))
        rows.append({"centromere": cen, "median_length": median,
                     "ref_seq": counts["cdeii_seq"].iloc[0], "ref_count": int(counts["_w"].iloc[0])})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ indels

def _aligner() -> PairwiseAligner:
    # Cheap extension so a 20-bp duplication stays one gap in AT-rich DNA
    return PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                           open_gap_score=-6, extend_gap_score=-0.5)


def _slide(seq: str, start: int, length: int) -> tuple[int, int]:
    """Leftmost and rightmost start of an extra segment seq[start:start+length]
    that removes to the same shorter sequence."""
    left = start
    while left > 0 and seq[left - 1] == seq[left + length - 1]:
        left -= 1
    right = start
    while right + length < len(seq) and seq[right] == seq[right + length]:
        right += 1
    return left, right


def _edit_distance(a: str, b: str) -> int:
    """Levenshtein distance (a copy shifted by one base inside the unit costs 2, not L)."""
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]


def _smallest_period(seg: str) -> int:
    for u in range(1, len(seg) + 1):
        if len(seg) % u == 0 and seg == seg[:u] * (len(seg) // u):
            return u
    return len(seg)


def describe_extra(seq: str, start: int, length: int) -> dict:
    """Tandem status, unit, orientation and microhomology of seq[start:start+length]."""
    left, right = _slide(seq, start, length)
    slide = right - left
    seg = seq[left:left + length]
    out = {"left": left, "slide": slide, "extra_seq": seg}
    if slide >= length:
        unit = _smallest_period(seg)
        out.update(tandem="exact", orientation="direct", unit_length=unit,
                   copies_added=length // unit, microhomology=slide - length)
        return out

    out.update(unit_length=length, copies_added=1, microhomology=slide)
    # Imperfect copies: <= 20 % edits for >= 8 bp; for 5-7 bp only an exact inverted copy
    max_mm = length // 5 if length >= 8 else 0
    best = None  # (mismatches, orientation)
    for x in range(left, right + 1):
        s = seq[x:x + length]
        for nb in (seq[max(0, x - length):x], seq[x + length:x + 2 * length]):
            if len(nb) < length:
                continue
            for orient, cand in (("direct", s), ("inverted", revcomp(s))):
                mm = _edit_distance(cand, nb)
                if best is None or mm < best[0]:
                    best = (mm, orient)
    if length >= 5 and best is not None and best[0] <= max_mm:
        out.update(tandem="inverted" if best[1] == "inverted" else "imperfect", orientation=best[1])
    else:
        out.update(tandem="none", orientation=None)
    return out


def indel_events(ref: str, alt: str) -> list[dict]:
    """Indels of alt relative to ref. `ref_pos` is 0-based in ref (left-normalised);
    `length` > 0 for insertions, < 0 for deletions; mixed gaps are 'complex'."""
    if ref == alt:
        return []
    aln = _aligner().align(ref, alt)[0]
    tb, qb = aln.aligned
    # Sentinel blocks make leading and trailing gaps look like internal ones
    tb = [(0, 0)] + [tuple(b) for b in tb] + [(len(ref), len(ref))]
    qb = [(0, 0)] + [tuple(b) for b in qb] + [(len(alt), len(alt))]
    events = []
    for i in range(len(tb) - 1):
        t0, t1 = tb[i][1], tb[i + 1][0]
        q0, q1 = qb[i][1], qb[i + 1][0]
        dt, dq = t1 - t0, q1 - q0
        if dt == dq:  # mismatches only (or nothing)
            continue
        if dt and dq:
            events.append({"type": "complex", "length": dq - dt, "ref_pos": t0,
                           "ref_span": dt, "alt_span": dq, "tandem": "none",
                           "extra_seq": alt[q0:q1] if dq > dt else ref[t0:t1]})
            continue
        if dq:  # insertion: extra bases live in alt
            d = describe_extra(alt, q0, dq)
            ref_pos = t0 - (q0 - d["left"])
            events.append({"type": "insertion", "length": dq, "ref_pos": ref_pos, **d})
        else:   # deletion: extra bases live in ref
            d = describe_extra(ref, t0, dt)
            events.append({"type": "deletion", "length": -dt, "ref_pos": d["left"], **d})
    for e in events:
        e.pop("left", None)
        e["rel_pos"] = round(e["ref_pos"] / max(1, len(ref)), 3)
    return events


DIVERGENT_IDENTITY = 0.85


def identity(ref: str, alt: str) -> float:
    """Matches / gap-free columns of the global alignment, so indels do not lower it.

    Alleles below DIVERGENT_IDENTITY are flagged: their alignments to the reference
    often place compensating insertion/deletion pairs that are not real events.
    """
    aln = _aligner().align(ref, alt)[0]
    cols = matches = 0
    for x, y in zip(aln[0], aln[1]):
        if x != "-" and y != "-":
            cols += 1
            matches += x == y
    return matches / cols if cols else 0.0


# --------------------------------------------------------------- per dataset

def call_steps(loci: pd.DataFrame, s288c: dict[str, str], id_cols: list[str],
               weight: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Per-locus step table, per-event table and the reference alleles.

    `loci`: one row per allele with `centromere`, `cdeii_seq` and the `id_cols`.
    `s288c`: centromere -> full S288C CEN sequence.
    """
    refs = reference_alleles(loci, weight)
    ref_by_cen = refs.set_index("centromere")
    steps, events = [], []
    for _, row in loci.iterrows():
        cen, seq = row["centromere"], row["cdeii_seq"]
        ref = ref_by_cen.loc[cen]
        s288c_len = len(s288c[cen]) - CDEI_LEN - CDEIII_LEN
        ev = indel_events(ref["ref_seq"], seq)
        ids = {c: row[c] for c in id_cols}
        if weight:
            ids[weight] = row[weight]
        d_med = len(seq) - ref["median_length"]
        ident = identity(ref["ref_seq"], seq)
        big = [e for e in ev if abs(e["length"]) >= 3]
        steps.append({**ids, "centromere": cen, "cdeii_length": len(seq),
                      "median_length": ref["median_length"], "delta_median": d_med,
                      "class_median": step_class(d_med), "s288c_length": s288c_len,
                      "delta_s288c": len(seq) - s288c_len,
                      "class_s288c": step_class(len(seq) - s288c_len),
                      "identity_to_ref": round(ident, 3), "divergent": ident < DIVERGENT_IDENTITY,
                      "n_indels": len(ev), "n_indels_3bp": len(big),
                      "events": "; ".join(_event_label(e) for e in big)})
        for e in ev:
            events.append({**ids, "centromere": cen, "cdeii_length": len(seq),
                           "delta_median": d_med, "divergent": ident < DIVERGENT_IDENTITY, **e})
    return pd.DataFrame(steps), pd.DataFrame(events), refs


def _event_label(e: dict) -> str:
    kind = {"insertion": "ins", "deletion": "del", "complex": "cx"}[e["type"]]
    tail = ""
    if e.get("tandem") == "exact":
        tail = f" tandem x{e['copies_added']}u{e['unit_length']} mh{e['microhomology']}"
    elif e.get("tandem") in ("imperfect", "inverted"):
        tail = f" {e['tandem']}"
    return f"{kind}{abs(e['length'])}@{e['ref_pos']}{tail}"


# ------------------------------------------------------------ Helsen SD4

def helsen_cdeii(helsen: pd.DataFrame, extra_3p: int = 1,
                 min_len: int = 60, max_len: int = 130) -> pd.DataFrame:
    """Distinct CDEII alleles per centromere in Helsen SD4, with assembly counts.

    Helsen sequences are SGD CEN + `extra_3p` bp, so CDEII = seq[8 : -(26 + extra_3p)].
    Alleles outside min_len..max_len are dropped as annotation failures.
    """
    h = helsen.copy()
    h["cdeii_seq"] = [s[CDEI_LEN:len(s) - CDEIII_LEN - extra_3p] for s in h["sequence"]]
    h = h[h["cdeii_seq"].str.len().between(min_len, max_len)]
    h = h[~h["cdeii_seq"].str.contains("[^ACGT]")]
    out = (h.groupby(["centromere", "cdeii_seq"])
           .agg(n_assemblies=("assembly", "nunique"), example=("assembly", "first"))
           .reset_index())
    out["allele"] = out.groupby("centromere").cumcount().map(lambda i: f"a{i + 1}")
    out["allele"] = out["centromere"] + "_" + out["allele"]
    return out


# --- Analysis F: nucleotide diversity at the CEN vs its flanks ---------------------------

def region_map(s288c_aln: str, cen_start: int, cen_len: int, cdeiii_len: int = CDEIII_LEN,
               window: int = 100) -> pd.DataFrame:
    """Region of every MSA column from the S288C row.

    `cen_start` = position of CDEI's first base in the ungapped S288C window, `cen_len` = S288C
    CEN length (CDEI..CDEIII). Columns that are gaps in S288C get the preceding S288C position
    (leading ones the first). rel = position relative to CDEI start. Flank windows are counted
    outward from the CEN: left_1 = the 100 bp just before CDEI, right_1 = just after CDEIII.
    """
    is_base = np.array([c != "-" for c in s288c_aln])
    pos = np.cumsum(is_base) - 1
    pos[pos < 0] = 0
    rel = pos - cen_start
    region = np.empty(len(rel), dtype=object)
    offset = np.zeros(len(rel), dtype=int)
    cdeii_end = cen_len - cdeiii_len
    region[(rel >= 0) & (rel < CDEI_LEN)] = "CDEI"
    region[(rel >= CDEI_LEN) & (rel < cdeii_end)] = "CDEII"
    region[(rel >= cdeii_end) & (rel < cen_len)] = "CDEIII"
    left, right = rel < 0, rel >= cen_len
    offset[left] = -((-rel[left] - 1) // window + 1)
    offset[right] = (rel[right] - cen_len) // window + 1
    region[left] = [f"left_{-k}" for k in offset[left]]
    region[right] = [f"right_{k}" for k in offset[right]]
    return pd.DataFrame({"col": np.arange(len(rel)), "s288c_pos": pos, "rel": rel,
                         "s288c_base": is_base, "region": region, "offset": offset})


def _encode(seqs: list[str]) -> np.ndarray:
    """0 = missing (N or terminal gap), 1-4 ACGT, 5 = internal gap."""
    lut = np.zeros(256, dtype=np.int8)
    for i, b in enumerate("ACGT", 1):
        lut[ord(b)] = lut[ord(b.lower())] = i
    lut[ord("-")] = 5
    out = np.stack([lut[np.frombuffer(s.encode(), dtype=np.uint8)] for s in seqs])
    for row in out:
        bases = np.flatnonzero((row >= 1) & (row <= 4))
        if len(bases) == 0:
            row[:] = 0
            continue
        row[:bases[0]] = 0
        row[bases[-1] + 1:] = 0
    return out


def pairwise_pi(seqs: list[str], groups: np.ndarray, n_sites: dict | None = None,
                min_site_frac: float = 0.5, min_pair_frac: float = 0.5) -> pd.DataFrame:
    """Per group of columns: pi, pi_sub and indel_rate, averaged over sequence pairs.

    Per pair, comparable sites = columns where both carry a base; columns where either is
    missing (N, terminal gap) or both are gapped are skipped. A substitution is a comparable
    column with different bases; an indel event is a run of consecutive (non-skipped) columns
    where exactly one of the two has an internal gap, counted once at the column where it
    starts. pi = (substitutions + indel events) / comparable sites, pi_sub = substitutions /
    comparable sites, indel_rate = indel events / comparable sites.

    A pair enters a group only with comparable sites >= min_site_frac x n_sites[group]
    (n_sites default: the group's column count); a group with fewer than min_pair_frac of all
    pairs usable is NA.
    """
    labels, n_sites, use, pi_sub, indel = pair_rates(seqs, groups, n_sites, min_site_frac)
    n = len(seqs)
    n_pairs_total = n * (n - 1) // 2
    with np.errstate(invalid="ignore", divide="ignore"):
        n_use = use.sum(axis=0)
        ok = n_use >= min_pair_frac * n_pairs_total
        mean = lambda m: np.where(ok, np.nansum(m, axis=0) / np.maximum(n_use, 1), np.nan)
        return pd.DataFrame({"group": labels, "n_sites": [n_sites[g] for g in labels],
                             "n_pairs": n_use, "pi": mean(pi_sub + indel),
                             "pi_sub": mean(pi_sub), "indel_rate": mean(indel)})


def pair_rates(seqs: list[str], groups: np.ndarray, n_sites: dict | None = None,
               min_site_frac: float = 0.5):
    """Per pair (i < j, row order of `pair_index`) x group: usable mask, pi_sub, indel rate
    (NaN where unusable). Definitions as in `pairwise_pi`. Returns labels, n_sites, use,
    pi_sub, indel."""
    x = _encode(seqs)
    groups = np.asarray(groups)
    labels = pd.unique(groups)
    gi = pd.Index(labels).get_indexer(groups)
    ng = len(labels)
    n = len(seqs)
    if n_sites is None:
        n_sites = dict(zip(labels, np.bincount(gi, minlength=ng)))
    need = np.array([min_site_frac * n_sites[g] for g in labels])
    G = np.zeros((x.shape[1], ng), dtype=np.float32)       # column -> group one-hot
    G[np.arange(x.shape[1]), gi] = 1
    cols = np.arange(x.shape[1])
    acc = {k: [] for k in ("comp", "sub", "ind")}
    for i in range(n - 1):            # sequence i against all j > i at once
        a, b = x[i][None, :], x[i + 1:]
        keep = (a > 0) & (b > 0) & ~((a == 5) & (b == 5))
        comp = keep & (a < 5) & (b < 5)
        sub = comp & (a != b)
        state = np.where(keep & (a == 5), 1, np.where(keep & (b == 5), 2, 0)).astype(np.int8)
        # state at the previous kept column of the same pair (0 before the first one)
        last = np.maximum.accumulate(np.where(keep, cols, -1), axis=1)
        prev = np.concatenate([np.full((len(b), 1), -1), last[:, :-1]], axis=1)
        prev_state = np.where(prev >= 0, np.take_along_axis(state, np.maximum(prev, 0), axis=1), 0)
        opens = keep & (state > 0) & (state != prev_state)
        acc["comp"].append(comp.astype(np.float32) @ G)
        acc["sub"].append(sub.astype(np.float32) @ G)
        acc["ind"].append(opens.astype(np.float32) @ G)
    if n < 2:
        comp = sub = ind = np.zeros((0, ng))
    else:
        comp, sub, ind = (np.concatenate(acc[k]).astype(np.float64) for k in ("comp", "sub", "ind"))
    use = (comp >= need) & (comp > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = lambda m: np.where(use, m / np.where(comp > 0, comp, 1), np.nan)
        return labels, n_sites, use, rate(sub), rate(ind)


def pair_index(n: int) -> tuple[np.ndarray, np.ndarray]:
    """(i, j) of each pair in `pair_rates` row order."""
    i = np.repeat(np.arange(n - 1), np.arange(n - 1, 0, -1))
    j = np.concatenate([np.arange(k + 1, n) for k in range(n - 1)]) if n > 1 else np.array([], int)
    return i, j


def leave_one_out_pi(seqs: list[str], groups: np.ndarray, group, min_site_frac: float = 0.5,
                     min_pair_frac: float = 0.5) -> tuple[float, np.ndarray]:
    """pi of `group` (as `pairwise_pi`) for all sequences, and for each k without sequence k;
    one pass over the pairs instead of n + 1."""
    labels, _, use, pi_sub, indel = pair_rates(seqs, groups, min_site_frac=min_site_frac)
    g = list(labels).index(group)
    u, v = use[:, g], np.where(use[:, g], pi_sub[:, g] + indel[:, g], 0.0)
    n = len(seqs)
    i, j = pair_index(n)
    nu = np.bincount(i, u, n) + np.bincount(j, u, n)      # usable pairs touching k
    sv = np.bincount(i, v, n) + np.bincount(j, v, n)
    full_ok = u.sum() >= min_pair_frac * n * (n - 1) / 2
    full = v.sum() / max(u.sum(), 1) if full_ok else np.nan
    n_wo = u.sum() - nu
    ok = n_wo >= min_pair_frac * (n - 1) * (n - 2) / 2
    wo = np.where(ok, (v.sum() - sv) / np.maximum(n_wo, 1), np.nan)
    return full, wo


def diversity_windows(aln: dict[str, str], s288c: str, cen_start: int, cen_len: int,
                      window: int = 100, min_pair_frac: float = 0.5) -> pd.DataFrame:
    """pi per CDE and per 100-bp flank window from one window MSA (dict name -> aligned seq)."""
    rm = region_map(aln[s288c], cen_start, cen_len, window=window)
    n_sites = rm[rm["s288c_base"]].groupby("region").size().to_dict()
    res = pairwise_pi(list(aln.values()), rm["region"].to_numpy(), n_sites=n_sites,
                      min_pair_frac=min_pair_frac)
    off = rm.groupby("region")["offset"].first()
    res = res.rename(columns={"group": "region"})
    res.insert(1, "offset", res["region"].map(off).astype(int))
    return res.sort_values(["offset", "region"], key=lambda s: s if s.name == "offset" else
                           s.map({"CDEI": 0, "CDEII": 1, "CDEIII": 2}).fillna(0)).reset_index(drop=True)
