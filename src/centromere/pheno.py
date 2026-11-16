"""Centromere alleles vs phenotypes (Loegler 2025 traits) with an EMMAX-style LMM.

Kinship from the Peter 2018 1011 GWAS matrix (plink), leave-one-chromosome-out (LOCO) for a
CEN's own chromosome. Spec: docs/specs/2026-11-16-cen-phenotypes.md.
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

CDEI_LEN, CDEIII_LEN = 8, 26

# ------------------------------------------------------------------ plink I/O

# 2-bit code -> dosage of allele A1 (bim column 5); 01 = missing
_CODE_TO_DOSAGE = np.array([2, -1, 1, 0], dtype=np.int8)


def read_plink(prefix: str | Path) -> tuple[list[str], pd.DataFrame, np.ndarray]:
    """(sample ids, bim table, genotypes samples x markers, int8 A1 dosage, -1 = missing)."""
    prefix = str(prefix)
    fam = pd.read_csv(prefix + ".fam", sep=r"\s+", header=None)
    bim = pd.read_csv(prefix + ".bim", sep=r"\s+", header=None,
                      names=["chrom", "id", "cm", "pos", "a1", "a2"])
    n, m = len(fam), len(bim)
    raw = np.fromfile(prefix + ".bed", dtype=np.uint8)
    if raw[:3].tolist() != [0x6C, 0x1B, 0x01]:
        raise ValueError("not a SNP-major plink .bed file")
    nbytes = (n + 3) // 4
    body = raw[3:].reshape(m, nbytes)
    codes = np.stack([(body >> s) & 3 for s in (0, 2, 4, 6)], axis=2).reshape(m, nbytes * 4)[:, :n]
    return fam[1].astype(str).tolist(), bim, _CODE_TO_DOSAGE[codes].T.copy()


def write_plink(prefix: str | Path, samples: list[str], bim: pd.DataFrame, geno: np.ndarray) -> None:
    """Inverse of read_plink (used by tests)."""
    prefix = str(prefix)
    pd.DataFrame({0: samples, 1: samples, 2: 0, 3: 0, 4: 0, 5: -9}).to_csv(
        prefix + ".fam", sep=" ", header=False, index=False)
    bim[["chrom", "id", "cm", "pos", "a1", "a2"]].to_csv(prefix + ".bim", sep="\t", header=False, index=False)
    n, m = geno.shape
    code = np.select([geno == 2, geno == 1, geno == 0], [0, 2, 3], default=1).astype(np.uint8).T
    pad = (-n) % 4
    code = np.concatenate([code, np.zeros((m, pad), dtype=np.uint8)], axis=1).reshape(m, -1, 4)
    packed = code[:, :, 0] | (code[:, :, 1] << 2) | (code[:, :, 2] << 4) | (code[:, :, 3] << 6)
    with open(prefix + ".bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        fh.write(packed.astype(np.uint8).tobytes())


def cnv_states(geno: np.ndarray, a1: np.ndarray, a2: np.ndarray) -> np.ndarray:
    """Peter 2018 CNV markers -> copy state 0 (absent) / 1 (0.5-1 copy) / 2 (multiple copies).

    The allele labels in the .bim are the states themselves; state 0 is written as plink's
    missing allele, so a missing genotype means 'absent'. Heterozygous calls -> mean of a1, a2.
    """
    a1, a2 = a1.astype(float), a2.astype(float)
    return np.select([geno == 2, geno == 0, geno == 1], [a1, a2, (a1 + a2) / 2], default=0.0)


def local_pcs(geno: np.ndarray, n_pcs: int) -> np.ndarray:
    """Top principal components of standardised local SNPs (used to condition on the haplotype)."""
    x = standardise(geno)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    return u[:, :n_pcs] * s[:n_pcs]


# ------------------------------------------------------------------ kinship

def standardise(geno: np.ndarray) -> np.ndarray:
    """VanRaden standardisation; missing -> 0 (mean imputation); monomorphic markers -> 0."""
    g = geno.astype(np.float64)
    miss = geno < 0
    g[miss] = np.nan
    p = np.nanmean(g, axis=0) / 2
    sd = np.sqrt(2 * p * (1 - p))
    with np.errstate(invalid="ignore", divide="ignore"):
        x = (g - 2 * p) / sd
    x[:, ~(sd > 0)] = 0
    x[np.isnan(x)] = 0
    return x


def loco_kinships(geno: np.ndarray, chroms: np.ndarray) -> dict:
    """{chrom: K without that chromosome, 'all': K from all markers}; K = X X' / m."""
    num, cnt = {}, {}
    for c in np.unique(chroms):
        x = standardise(geno[:, chroms == c])
        num[c], cnt[c] = x @ x.T, x.shape[1]
    total, m = sum(num.values()), sum(cnt.values())
    out = {"all": total / m}
    for c in num:
        if m > cnt[c]:
            out[c] = (total - num[c]) / (m - cnt[c])
    return out


# ------------------------------------------------------------------ LMM (EMMAX)

def inverse_normal(y: np.ndarray) -> np.ndarray:
    r = stats.rankdata(y)
    return stats.norm.ppf((r - 0.5) / len(y))


def _reml(logdelta: float, s: np.ndarray, yt: np.ndarray, xt: np.ndarray) -> float:
    """Negative REML log-likelihood (Kang et al. 2008) at delta = exp(logdelta)."""
    n, p = xt.shape
    w = 1.0 / (s + np.exp(logdelta))
    xtw = xt.T * w
    a = xtw @ xt
    beta = np.linalg.solve(a, xtw @ yt)
    r = yt - xt @ beta
    rss = np.sum(w * r * r)
    _, logdet_a = np.linalg.slogdet(a)
    ll = 0.5 * ((n - p) * np.log((n - p) / (2 * np.pi)) - (n - p) - (n - p) * np.log(rss)
                + np.sum(np.log(w)) - logdet_a)
    return -ll


def eigen(k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    s, u = np.linalg.eigh(k)
    return np.clip(s, 0, None), u


def fit_null(y: np.ndarray, x: np.ndarray, s: np.ndarray, u: np.ndarray) -> dict:
    """REML fit of y = X b + g + e; returns delta, h2 and the rotated data for lmm_test."""
    yt, xt = u.T @ y, u.T @ x
    grid = np.linspace(-10, 10, 81)
    vals = [_reml(g, s, yt, xt) for g in grid]
    i = int(np.argmin(vals))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    best = optimize.minimize_scalar(_reml, bounds=(lo, hi), args=(s, yt, xt), method="bounded")
    logd = best.x if best.fun <= vals[i] else grid[i]
    delta = float(np.exp(logd))
    return {"delta": delta, "h2": 1.0 / (1.0 + delta), "w": 1.0 / (s + delta),
            "yt": yt, "xt": xt, "u": u}


def lmm_test(null: dict, pred: np.ndarray) -> dict:
    """GLS Wald test of one predictor with the null model's delta (EMMAX)."""
    xt = np.column_stack([null["xt"], null["u"].T @ pred])
    w, yt = null["w"], null["yt"]
    n, p = xt.shape
    xtw = xt.T * w
    a = xtw @ xt
    try:
        ainv = np.linalg.inv(a)
    except np.linalg.LinAlgError:
        return {"beta": np.nan, "se": np.nan, "p": np.nan}
    beta = ainv @ (xtw @ yt)
    r = yt - xt @ beta
    sigma2 = np.sum(w * r * r) / (n - p)
    se = float(np.sqrt(sigma2 * ainv[-1, -1]))
    t = beta[-1] / se
    return {"beta": float(beta[-1]), "se": se, "p": float(2 * stats.t.sf(abs(t), n - p))}


def ols_test(y: np.ndarray, x: np.ndarray, pred: np.ndarray) -> float:
    """Plain OLS p for comparison in the calibration test."""
    xx = np.column_stack([x, pred])
    beta, *_ = np.linalg.lstsq(xx, y, rcond=None)
    r = y - xx @ beta
    n, p = xx.shape
    se = np.sqrt(r @ r / (n - p) * np.linalg.inv(xx.T @ xx)[-1, -1])
    return float(2 * stats.t.sf(abs(beta[-1] / se), n - p))


def genomic_lambda(p: pd.Series) -> float:
    p = p.dropna()
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)) if len(p) else np.nan


def bh(p: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p.index)
    ok = p.dropna()
    if len(ok):
        from statsmodels.stats.multitest import multipletests
        out[ok.index] = multipletests(ok.values, method="fdr_bh")[1]
    return out


# ------------------------------------------------------------------ CEN predictors

def _off_step(delta: int) -> bool:
    d = abs(delta)
    return d >= 3 and not any(abs(d - k) <= 2 for k in (10, 20))


def cen_predictors(helsen: pd.DataFrame, extra_3p: int = 1, rare_freq: float = 0.05) -> pd.DataFrame:
    """strain x CEN predictors from Helsen SD4 long table (assembly = strain code).

    non_modal = share of the strain's variants differing from the population major allele;
    cdeii_delta = mean CDEII length - per-CEN median; off_step = any variant off a whole turn;
    rare = any variant with population frequency < rare_freq. Variants with non-ACGT are dropped.
    """
    h = helsen[~helsen["sequence"].str.contains("[^ACGT]")].copy()
    h["cdeii_len"] = h["sequence"].str.len() - CDEI_LEN - CDEIII_LEN - extra_3p
    freq = h.groupby(["centromere", "sequence"]).size().rename("count").reset_index()
    freq["freq"] = freq["count"] / freq.groupby("centromere")["count"].transform("sum")
    modal = freq.sort_values("count", ascending=False).drop_duplicates("centromere")
    h = h.merge(freq[["centromere", "sequence", "freq"]], on=["centromere", "sequence"])
    h = h.merge(modal[["centromere", "sequence"]].rename(columns={"sequence": "modal"}), on="centromere")
    median = h[h["variant"] == "A"].groupby("centromere")["cdeii_len"].median()
    h["delta"] = h["cdeii_len"] - h["centromere"].map(median)
    h["is_off"] = h["delta"].map(_off_step)
    g = h.groupby(["assembly", "centromere"])
    out = pd.DataFrame({
        "non_modal": g.apply(lambda d: float((d["sequence"] != d["modal"]).mean()), include_groups=False),
        "cdeii_delta": g["delta"].mean(),
        "off_step": g["is_off"].any().astype(float),
        "rare": g["freq"].min().lt(rare_freq).astype(float),
        "het": g["sequence"].nunique().gt(1).astype(float),
    }).reset_index().rename(columns={"assembly": "strain"})
    return out


def rare_burden(pred: pd.DataFrame, n_cen: int = 16) -> pd.Series:
    """Strains with all CENs called -> number of CENs carrying a rare allele."""
    g = pred.groupby("strain")
    full = g.size() == n_cen
    return g["rare"].sum()[full]


# ------------------------------------------------------------------ genes

def read_gff_genes(path: str | Path) -> pd.DataFrame:
    """Protein-coding genes: chrom, start, end, orf (systematic id), name (standard or orf)."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene" or "protein_coding_gene" not in f[8]:
                continue
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            rows.append((f[0], int(f[3]), int(f[4]), attr["ID"], attr.get("gene", attr["ID"])))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "orf", "name"])


def cis_genes(genes: pd.DataFrame, cens: pd.DataFrame, window: int) -> pd.DataFrame:
    """Genes overlapping CEN +- window; cens has chrom, start, end, centromere (BED, 0-based start)."""
    out = []
    for c in cens.itertuples():
        g = genes[(genes["chrom"] == c.chrom) & (genes["end"] > c.start - window)
                  & (genes["start"] <= c.end + window)]
        out.append(g.assign(centromere=c.centromere))
    return pd.concat(out, ignore_index=True)


# ------------------------------------------------------------------ phenotypes

def read_phen(path: str | Path) -> pd.Series:
    d = pd.read_csv(path, sep="\t")
    s = pd.to_numeric(d.iloc[:, 1], errors="coerce")
    s.index = d.iloc[:, 0].astype(str)
    return s.dropna()[lambda x: ~x.index.duplicated()]


def list_traits(root: str | Path) -> pd.DataFrame:
    """trait table: trait, klass (growth/expression/protein), orf (None for growth), path."""
    rows = []
    for klass, sub in (("growth", "GrowthTraits"), ("expression", "ExpressionTraits"),
                       ("protein", "ProteinTraits")):
        for p in sorted(Path(root, sub).glob("*.phen")):
            name = p.name[:-5]
            orf = None if klass == "growth" else name.split(".")[0]
            rows.append((name, klass, orf, str(p)))
    return pd.DataFrame(rows, columns=["trait", "klass", "orf", "path"])


ROMAN = {r: i for i, r in enumerate(
    "I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI".split(), 1)}


def chrom_number(name: str) -> int:
    return ROMAN[re.sub("^chr", "", name)]
