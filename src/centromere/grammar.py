"""Analysis F2: CDEII A-tract grammar and predicted DNA shape.

Spec docs/specs/2026-11-15-cdeii-grammar-shape.md.

A-tract = maximal A/T run of >= `min_len` bp without a TpA step (AnTm form), so
AAAATTTT is one tract and TTTTAAAA is two. Phasing = fraction of tract-centre pairs
spaced within +-tol of a multiple of the helical period. Shape comes from DNAshapeR
(pentamer model), run once from R on a FASTA of every sequence; values for the two
terminal bases (MGW, ProT) or steps (Roll, HelT) are NA.
"""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

SHAPES = ["MGW", "Roll", "HelT", "ProT"]
_AT = set("AT")
_COMP = str.maketrans("ACGT", "TGCA")


# ------------------------------------------------------------------ grammar

def atracts(seq: str, min_len: int = 4) -> list[tuple[int, int]]:
    """(start, length) of every A-tract, 0-based."""
    out, start = [], None
    for i, b in enumerate(seq + "G"):
        breaks_ta = start is not None and b == "A" and seq[i - 1] == "T"
        if b in _AT and not breaks_ta:
            if start is None:
                start = i
            continue
        if start is not None and i - start >= min_len:
            out.append((start, i - start))
        start = i if breaks_ta else None
    return out


def phase_score(tracts: list[tuple[int, int]], period: float = 10.5, tol: float = 1.0,
                max_d: float = 42.0) -> float:
    """Fraction of tract-centre pairs (distance <= max_d) within tol of k * period, k >= 1."""
    c = np.array([s + (n - 1) / 2 for s, n in tracts])
    if len(c) < 2:
        return np.nan
    d = np.abs(c[:, None] - c[None, :])[np.triu_indices(len(c), 1)]
    d = d[d <= max_d]
    if len(d) == 0:
        return np.nan
    k = np.maximum(np.round(d / period), 1)
    return float(np.mean(np.abs(d - k * period) <= tol))


def aatt_power(seq: str, periods=(10.0, 10.25, 10.5, 10.75, 11.0)) -> float:
    """Mean normalised periodogram power of the AA/TT dinucleotide indicator at the periods.

    Normalised so that white noise has expectation 1.
    """
    x = np.array([seq[i:i + 2] in ("AA", "TT") for i in range(len(seq) - 1)], float)
    v = x.var()
    if v == 0:
        return np.nan
    x = x - x.mean()
    i = np.arange(len(x))
    p = [abs(np.sum(x * np.exp(-2j * np.pi * i / T))) ** 2 / (len(x) * v) for T in periods]
    return float(np.mean(p))


def grammar_stats(seq: str, min_len: int = 4, period: float = 10.5, tol: float = 1.0,
                  max_d: float = 42.0) -> dict:
    t = atracts(seq, min_len)
    lens = [n for _, n in t]
    return {"length": len(seq), "at": sum(b in _AT for b in seq) / len(seq), "n_tracts": len(t),
            "mean_tract": float(np.mean(lens)) if lens else 0.0, "max_tract": max(lens, default=0),
            "frac_in_tracts": sum(lens) / len(seq), "phase": phase_score(t, period, tol, max_d),
            "power": aatt_power(seq)}


def empirical_p(null: np.ndarray, obs: float) -> float:
    """Two-sided empirical p with +1 correction; NaN values in the null are dropped."""
    null = np.asarray(null, float)
    null = null[~np.isnan(null)]
    if np.isnan(obs) or len(null) == 0:
        return np.nan
    lo = (np.sum(null <= obs) + 1) / (len(null) + 1)
    hi = (np.sum(null >= obs) + 1) / (len(null) + 1)
    return float(min(1.0, 2 * min(lo, hi)))


def percentile(null: np.ndarray, obs: float) -> float:
    """Mid-rank percentile of obs in the null (0 = below all, 1 = above all)."""
    null = np.asarray(null, float)
    null = null[~np.isnan(null)]
    if np.isnan(obs) or len(null) == 0:
        return np.nan
    return float((np.sum(null < obs) + 0.5 * np.sum(null == obs)) / len(null))


# ------------------------------------------------------- AT-matched windows

class GenomeWindows:
    """AT-matched windows of a genome, avoiding blocked intervals and non-ACGT bases."""

    def __init__(self, genome: dict[str, str], blocked: dict[str, list[tuple[int, int]]]):
        self.names = list(genome)
        self.seqs = genome
        self.at_cs, self.bad_cs = {}, {}
        for c, s in genome.items():
            a = np.frombuffer(s.upper().encode(), np.uint8)
            is_at = (a == ord("A")) | (a == ord("T"))
            bad = ~(is_at | (a == ord("C")) | (a == ord("G")))
            for lo, hi in blocked.get(c, []):
                bad[max(lo, 0):min(hi, len(a))] = True
            self.at_cs[c] = np.concatenate([[0], np.cumsum(is_at, dtype=np.int32)])
            self.bad_cs[c] = np.concatenate([[0], np.cumsum(bad, dtype=np.int32)])
        self._cache: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    def _windows(self, length: int):
        if length not in self._cache:
            ci, st, at = [], [], []
            for k, c in enumerate(self.names):
                cs, bc = self.at_cs[c], self.bad_cs[c]
                if len(cs) - 1 < length:
                    continue
                ok = (bc[length:] - bc[:-length]) == 0
                s = np.nonzero(ok)[0]
                ci.append(np.full(len(s), k, np.int16))
                st.append(s.astype(np.int32))
                at.append((cs[length:] - cs[:-length])[s].astype(np.int16))
            self._cache = {length: (np.concatenate(ci), np.concatenate(st), np.concatenate(at))}
        return self._cache[length]

    def sample(self, length: int, at: float, tol: float, n: int, rng: np.random.Generator,
               tol_wide: float | None = None, min_n: int = 0) -> tuple[list[str], float]:
        """Up to n windows with AT within tol of `at` (widened to tol_wide if < min_n),
        random strand. Returns the sequences and the tolerance used."""
        ci, st, cnt = self._windows(length)
        for t in [tol] + ([tol_wide] if tol_wide else []):
            sel = np.nonzero(np.abs(cnt / length - at) <= t + 1e-9)[0]
            if len(sel) >= min_n:
                break
        pick = rng.choice(sel, size=min(n, len(sel)), replace=False)
        out = []
        for i in pick:
            s = self.seqs[self.names[ci[i]]][st[i]:st[i] + length].upper()
            out.append(s.translate(_COMP)[::-1] if rng.random() < 0.5 else s)
        return out, t


# ------------------------------------------------------ substitution null

def substitutions(ref: str, alt: str) -> list[tuple[int, str, str]]:
    if len(ref) != len(alt):
        raise ValueError("substitutions need equal-length sequences")
    return [(i, a, b) for i, (a, b) in enumerate(zip(ref, alt)) if a != b]


def spectrum_null(ref: str, subs: list[tuple[int, str, str]], n: int,
                  rng: np.random.Generator) -> list[str]:
    """n mutants of ref with the same substitutions (b1 -> b2) at random positions that
    carry b1; no position is used twice within a mutant."""
    by_base = {b: [i for i, x in enumerate(ref) if x == b] for b in "ACGT"}
    out = []
    for _ in range(n):
        s, used = list(ref), set()
        for _, b1, b2 in subs:
            free = [i for i in by_base[b1] if i not in used]
            i = free[rng.integers(len(free))]
            used.add(i)
            s[i] = b2
        out.append("".join(s))
    return out


def mutate(ref: str, pos: int, base: str) -> str:
    return ref[:pos] + base + ref[pos + 1:]


def event_null_positions(ref: str, pos: int, context: bool = False) -> list[int]:
    """Exhaustive null for one substitution at `pos`: every other position carrying the same
    base, or (context=True) the same trinucleotide (ends use the available neighbour only)."""
    b = ref[pos]

    def ctx(i):
        return ref[max(i - 1, 0):i + 2] if 0 < i < len(ref) - 1 else None

    if context:
        c = ctx(pos)
        return [] if c is None else [i for i in range(len(ref)) if i != pos and ctx(i) == c]
    return [i for i, x in enumerate(ref) if x == b and i != pos]


# ------------------------------------------------------------------ shape

def write_fasta(seqs: dict[str, str], path: Path) -> None:
    with open(path, "w") as f:
        for k, s in seqs.items():
            f.write(f">{k}\n{s}\n")


def parse_shape_file(path: Path) -> dict[str, np.ndarray]:
    """DNAshapeR output: FASTA-like, comma-separated values over wrapped lines, NA allowed."""
    out, name, buf = {}, None, []

    def flush():
        if name is not None:
            vals = ",".join(buf).replace(",,", ",").strip(",").split(",")
            out[name] = np.array([np.nan if v == "NA" else float(v) for v in vals])

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                flush()
                name, buf = line[1:], []
            elif line:
                buf.append(line)
    flush()
    return out


def run_dnashape(seqs: dict[str, str], workdir: Path, rscript: str = "Rscript",
                 script: str = "scripts/dnashape.R") -> dict[str, dict[str, np.ndarray]]:
    """Shape per feature per sequence id; sequences are deduplicated before the R call."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    uniq = {s: f"u{i}" for i, s in enumerate(dict.fromkeys(seqs.values()))}
    fa = workdir / "shape_input.fasta"
    write_fasta({v: k for k, v in uniq.items()}, fa)
    subprocess.run([rscript, script, str(fa)], check=True, capture_output=True)
    res = {}
    for f in SHAPES:
        vals = parse_shape_file(Path(f"{fa}.{f}"))
        res[f] = {k: vals[uniq[s]] for k, s in seqs.items()}
    return res


def shape_distance(a: np.ndarray, b: np.ndarray) -> float:
    ok = ~(np.isnan(a) | np.isnan(b))
    return float(np.mean(np.abs(a[ok] - b[ok]))) if ok.any() else np.nan


def shape_summary(mgw: np.ndarray, prot: np.ndarray, narrow: float = 4.0) -> dict:
    m = mgw[~np.isnan(mgw)]
    return {"mean_mgw": float(m.mean()), "frac_narrow": float(np.mean(m < narrow)),
            "mean_prot": float(np.nanmean(prot))}


def shape_long(shape: dict[str, dict[str, np.ndarray]], ids: list[str]) -> pd.DataFrame:
    """Long table id, feature, pos (1-based; step i = between bases i and i+1), value."""
    rows = []
    for f in SHAPES:
        for k in ids:
            v = shape[f][k]
            rows.append(pd.DataFrame({"id": k, "feature": f, "pos": np.arange(1, len(v) + 1), "value": v}))
    return pd.concat(rows, ignore_index=True).dropna(subset=["value"])
