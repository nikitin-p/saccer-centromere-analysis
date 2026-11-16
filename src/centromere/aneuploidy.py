"""Analysis G: aneuploidy (Peter 2018) vs centromere alleles on the 748-strain BLAST set."""
import re

import numpy as np
import pandas as pd

_EVENT = re.compile(r"([+-])(\d+)\*(\d+)")


def parse_aneuploidies(value) -> list[tuple[int, int]] | None:
    """'aneu;+1*9;-1*3;' -> [(9, +1), (3, -1)]; 'euploid' -> []; missing -> None."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    value = str(value).strip()
    if value == "euploid":
        return []
    if not value.startswith("aneu"):
        return None
    return [(int(c), int(n) * (1 if s == "+" else -1)) for s, n, c in _EVENT.findall(value)]


def normalize_clade(value) -> str:
    """'1. Wine/European (subclade 3)' -> '1. Wine/European'; trailing spaces removed; NaN -> 'unassigned'."""
    if not isinstance(value, str) or not value.strip():
        return "unassigned"
    return re.sub(r"\s*\(subclade \d+\)\s*$", "", value).strip()


def load_peter2018(path) -> pd.DataFrame:
    """Strain metadata: strain, clade, ploidy, aneuploidies (raw), events (parsed list or None)."""
    meta = pd.read_excel(path, header=3)
    meta = meta.loc[meta["Standardized name"].notna()]
    out = pd.DataFrame({
        "strain": meta["Standardized name"].astype(str).str.strip(),
        "clade": meta["Clades"].map(normalize_clade),
        "ploidy": meta["Ploidy"],
        "aneuploidies": meta["Aneuploidies"],
    })
    out["events"] = out["aneuploidies"].map(parse_aneuploidies)
    return out.reset_index(drop=True)


def chromosome_status(meta: pd.DataFrame, chroms=range(1, 17)) -> pd.DataFrame:
    """strain x chromosome copy change (0 = euploid for that chromosome); strains with unknown status dropped."""
    known = meta.loc[meta["events"].notna()]
    rows = []
    for strain, events in zip(known["strain"], known["events"]):
        change = dict.fromkeys(chroms, 0)
        for c, n in events:
            change[c] = change.get(c, 0) + n
        rows += [{"strain": strain, "chrom": c, "change": change[c]} for c in chroms]
    return pd.DataFrame(rows)


PREDICTORS = ("non_modal", "length_variant", "non_s288c")


def locus_calls(loci: pd.DataFrame, reference: dict[str, str], end_tolerance: int = 2) -> pd.DataFrame:
    """Per strain x centromere from BLAST loci (sequences in query orientation).

    BLAST trims terminal mismatches (S288C CEN1 base 1 is G, most strains differ), so a
    locus is complete when the hit reaches within `end_tolerance` bp of both query ends;
    trimmed ends count as substitutions (length-neutral). Only complete loci are called:
      non_s288c       differs from S288C (inside the hit, or a trimmed end)
      non_modal       differs from the most common complete sequence of that CEN
      length_variant  full length (hit + trimmed ends) differs from the modal full length
    A strain with several complete loci is a variant if any of them is.
    """
    loci = loci.copy()
    ref_u = {c: s.upper() for c, s in reference.items()}
    loci["ref_len"] = loci["centromere"].map(lambda c: len(ref_u[c]))
    trim5, trim3 = loci["qstart"] - 1, loci["ref_len"] - loci["qend"]
    loci["complete"] = (trim5 <= end_tolerance) & (trim3 <= end_tolerance)
    loci["seq_u"] = loci["sequence"].str.upper()
    inside_ref = [ref_u[c][s - 1:e] for c, s, e in zip(loci["centromere"], loci["qstart"], loci["qend"])]
    loci["is_ref"] = (loci["seq_u"] == inside_ref) & (trim5 == 0) & (trim3 == 0)
    loci["full_len"] = loci["sequence_length"] + trim5 + trim3
    comp = loci.loc[loci["complete"]]
    modal = pd.DataFrame({
        "modal_len": comp.groupby("centromere")["full_len"].agg(lambda s: s.mode().min()),
        "modal_seq": comp.groupby("centromere")["seq_u"].agg(lambda s: s.value_counts().index[0]),
    })
    loci = loci.join(modal, on="centromere")

    def summarise(g: pd.DataFrame) -> pd.Series:
        c = g.loc[g["complete"]]
        if c.empty:
            return pd.Series({"n_loci": len(g), "n_complete": 0, **dict.fromkeys(PREDICTORS, np.nan),
                              "delta_len": np.nan, "n_alleles": 0})
        delta = c["full_len"] - c["modal_len"]
        return pd.Series({
            "n_loci": len(g), "n_complete": len(c),
            "non_modal": float((c["seq_u"] != c["modal_seq"]).any()),
            "length_variant": float((delta != 0).any()),
            "non_s288c": float((~c["is_ref"]).any()),
            "delta_len": int(delta.iloc[delta.abs().argmax()]),
            "n_alleles": c["seq_u"].nunique(),
        })

    out = loci.groupby(["strain", "centromere"]).apply(summarise, include_groups=False).reset_index()
    out["chrom"] = out["centromere"].str.removeprefix("CEN").astype(int)
    return out


def _bh(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg q-values, NaN kept."""
    ok = p.dropna().sort_values()
    m = len(ok)
    q = (ok * m / np.arange(1, m + 1))[::-1].cummin()[::-1].clip(upper=1)
    return q.reindex(p.index)


def stratified_permutation(y: np.ndarray, x: np.ndarray, strata: np.ndarray,
                           n_perm: int = 10000, seed: int = 0) -> float:
    """Two-sided p for sum(x*y), permuting x within strata (the conditional-logit score statistic)."""
    rng = np.random.default_rng(seed)
    obs = float((x * y).sum())
    expected = 0.0
    null = np.zeros(n_perm)
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s)
        xs, ys = x[idx], y[idx]
        expected += xs.mean() * ys.sum()
        if xs.min() == xs.max() or ys.min() == ys.max():
            null += (xs * ys).sum()
            continue
        perm = rng.permuted(np.tile(xs, (n_perm, 1)), axis=1)
        null += perm @ ys
    return float((np.sum(np.abs(null - expected) >= abs(obs - expected) - 1e-9) + 1) / (n_perm + 1))


def test_chromosome(df: pd.DataFrame, predictor: str, n_perm: int = 10000) -> dict:
    """df: one row per strain with y (aneuploid for this chromosome), predictor (0/1), clade."""
    from scipy.stats import fisher_exact
    from statsmodels.discrete.conditional_models import ConditionalLogit
    from statsmodels.stats.contingency_tables import StratifiedTable

    y = df["y"].to_numpy(int)
    x = df[predictor].to_numpy(int)
    strata = df["clade"].to_numpy()
    tab = pd.crosstab(pd.Categorical(y, [0, 1]), pd.Categorical(x, [0, 1]), dropna=False).to_numpy()
    res = {
        "n": len(df), "n_aneuploid": int(y.sum()), "n_variant": int(x.sum()),
        "aneu_variant": int(tab[1, 1]), "aneu_ref": int(tab[1, 0]),
        "eu_variant": int(tab[0, 1]), "eu_ref": int(tab[0, 0]),
        "fisher_or": np.nan, "fisher_p": np.nan, "mh_or": np.nan,
        "clogit_or": np.nan, "clogit_lo": np.nan, "clogit_hi": np.nan, "clogit_p": np.nan,
        "perm_p": np.nan, "separation": False,
    }
    if x.min() == x.max() or y.min() == y.max():
        return res
    res["fisher_or"], res["fisher_p"] = fisher_exact(tab)

    # strata where both y and x vary carry information; the rest drop out of the conditional likelihood
    informative = [s for s in np.unique(strata)
                   if np.ptp(y[strata == s]) > 0 and np.ptp(x[strata == s]) > 0]
    tables = [pd.crosstab(pd.Categorical(x[strata == s], [0, 1]),
                          pd.Categorical(y[strata == s], [0, 1]), dropna=False).to_numpy()
              for s in informative]
    if tables:
        res["mh_or"] = StratifiedTable(tables).oddsratio_pooled
        res["perm_p"] = stratified_permutation(y, x, strata, n_perm)
        keep = np.isin(strata, informative)
        pooled = np.sum(tables, axis=0)
        if pooled.min() == 0:  # quasi-separation: conditional MLE diverges, keep only perm_p / MH
            res["separation"] = True
            return res
        try:
            fit = ConditionalLogit(y[keep], x[keep, None].astype(float), groups=strata[keep]).fit(disp=0)
            beta, (lo, hi) = fit.params[0], fit.conf_int()[0]
            res.update(clogit_or=np.exp(beta), clogit_lo=np.exp(lo), clogit_hi=np.exp(hi),
                       clogit_p=fit.pvalues[0])
        except (np.linalg.LinAlgError, ValueError):
            pass
    return res


def association(calls: pd.DataFrame, status: pd.DataFrame, meta: pd.DataFrame,
                min_aneuploid: int = 10, n_perm: int = 10000) -> pd.DataFrame:
    """Per chromosome k and predictor: aneuploid for k (gain or loss) ~ allele class at CENk | clade."""
    data = (calls.merge(status, on=["strain", "chrom"])
                 .merge(meta[["strain", "clade"]], on="strain"))
    data = data.loc[data["n_complete"] > 0].copy()
    data["y"] = (data["change"] != 0).astype(int)
    data["gain"] = (data["change"] > 0).astype(int)
    rows = []
    for k, d in data.groupby("chrom"):
        for predictor in PREDICTORS:
            r = {"chrom": k, "centromere": f"CEN{k}", "predictor": predictor,
                 "n_gain": int(d["gain"].sum()), "tested": int(d["y"].sum()) >= min_aneuploid}
            r.update(test_chromosome(d, predictor, n_perm))
            rows.append(r)
    out = pd.DataFrame(rows)
    out["perm_q"] = np.nan
    for predictor, idx in out.loc[out["tested"]].groupby("predictor").groups.items():
        out.loc[idx, "perm_q"] = _bh(out.loc[idx, "perm_p"])
    return out
