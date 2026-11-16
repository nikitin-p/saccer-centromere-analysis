#!/usr/bin/env python3
"""CEN alleles (Helsen SD4, Peter 2018 genomes) vs Loegler 2025 phenotypes, EMMAX LMM with LOCO kinship.

Outputs in <output.dir>/phenotypes/:
  cen_predictors.tsv     strain x CEN predictors (+ rare_burden per strain)
  trait_samples.tsv      trait, class, family, n, null h2 (full K)
  assoc_primary.tsv      hypothesis-driven tests (segregation-stress growth, cis genes, kinetochore genes)
  assoc_secondary.tsv    all growth traits x all predictors
  positive_controls.tsv  CUP1 CNV vs copper growth; ORF CNV vs own expression

  PYTHONPATH=src .venv/bin/python scripts/15_phenotypes.py --config config.peter748.yaml
"""
import argparse
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from centromere.aneuploidy import chromosome_status, load_peter2018
from centromere.io import load_config, save_tsv
from centromere.pheno import (bh, chrom_number, cis_genes, cen_predictors, cnv_states, eigen, fit_null,
                              genomic_lambda, inverse_normal, list_traits, lmm_test, local_pcs,
                              loco_kinships, rare_burden, read_gff_genes, read_phen, read_plink)
from centromere.validate import load_helsen_sd4

CEN_PREDICTORS = ("non_modal", "cdeii_delta", "off_step")


class Model:
    """Caches eigendecompositions per (kinship key, sample set) and null fits per trait."""

    def __init__(self, kin: dict, kin_ids: list[str]):
        self.kin, self.idx = kin, {s: i for i, s in enumerate(kin_ids)}
        self._eig = {}

    def eig(self, key, strains: tuple):
        ck = (key, hash(strains))
        if ck not in self._eig:
            ii = [self.idx[s] for s in strains]
            self._eig[ck] = eigen(self.kin[key][np.ix_(ii, ii)])
        return self._eig[ck]


def minor_count(x: np.ndarray, predictor: str) -> int:
    nz = int(np.sum(x != 0))
    return nz if predictor == "cdeii_delta" else min(nz, len(x) - nz)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    config = load_config(args.config)
    cfg = config["phenotypes"]
    out = Path(config["output"]["dir"]) / "phenotypes"
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # --- genotypes and kinship
    kin_ids, bim, geno = read_plink(cfg["plink"])
    snp = (bim["chrom"] >= 1) & (bim["chrom"] <= 16)
    kin = loco_kinships(geno[:, snp.values], bim.loc[snp, "chrom"].values)
    is_cnv = (bim["chrom"] == 17).values
    cnv = pd.DataFrame(cnv_states(geno[:, is_cnv], bim.loc[is_cnv, "a1"].values, bim.loc[is_cnv, "a2"].values),
                       index=kin_ids, columns=bim.loc[is_cnv, "id"].str.replace("^rs", "", regex=True))
    model = Model(kin, kin_ids)
    print(f"kinship: {len(kin_ids)} strains, {int(snp.sum())} SNPs, {cnv.shape[1]} CNV markers "
          f"({time.time() - t0:.0f} s)")

    # --- CEN predictors
    helsen = load_helsen_sd4(cfg["helsen_sd4"])
    helsen = helsen[helsen["assembly"].isin(kin_ids)]
    pred = cen_predictors(helsen, rare_freq=cfg["rare_freq"])
    burden = rare_burden(pred)
    save_tsv(pred.merge(burden.rename("rare_burden"), left_on="strain", right_index=True, how="left"),
             out / "cen_predictors.tsv")
    complete = set(burden.index)
    wide = {p: pred.pivot(index="strain", columns="centromere", values=p) for p in CEN_PREDICTORS}
    print(f"CEN predictors: {pred['strain'].nunique()} strains, {len(complete)} with all 16 CENs")

    # --- covariates
    meta = load_peter2018(cfg["metadata"]).drop_duplicates("strain").set_index("strain")
    status = chromosome_status(meta.reset_index())
    change = status.pivot(index="strain", columns="chrom", values="change")
    ploidy = pd.to_numeric(meta["ploidy"], errors="coerce")
    ploidy = ploidy.fillna(ploidy.median())
    n_aneu = (change != 0).sum(axis=1)
    print(f"covariates: {len(meta)} strains in Peter 2018, {meta['events'].isna().sum()} with unknown "
          "aneuploidy status (counted as euploid)")

    def covariates(strains, chrom):
        c = np.column_stack([np.ones(len(strains)), ploidy.reindex(strains).fillna(ploidy.median())])
        extra = (n_aneu if chrom == "all" else change[chrom]).reindex(strains).fillna(0).values
        return np.column_stack([c, extra])

    # --- traits
    traits = list_traits(cfg["traits_dir"])
    genes = read_gff_genes(cfg["gff"])
    cens = pd.read_csv(cfg["cen_bed"], sep="\t", header=None,
                       names=["chrom", "start", "end", "centromere", "score", "strand"])
    cis = cis_genes(genes, cens, cfg["cis_window"])
    kt_orfs = genes.loc[genes["name"].isin(cfg["kinetochore_genes"]), "orf"].tolist()
    missing_kt = set(cfg["kinetochore_genes"]) - set(genes.loc[genes["orf"].isin(kt_orfs), "name"])
    if missing_kt:
        print("kinetochore genes not in GFF:", missing_kt)
    stress = re.compile(cfg["stress_regex"])
    cen_chrom = {c.centromere: chrom_number(c.chrom) for c in cens.itertuples()}
    cens16 = sorted(cen_chrom, key=lambda c: cen_chrom[c])

    # tests: trait -> list of (family, CEN or 'burden', predictor)
    plan: dict[str, list] = {}
    for t in traits.itertuples():
        tests = []
        if t.klass == "growth":
            full = [(c, p) for c in cens16 for p in CEN_PREDICTORS] + [("burden", "rare_burden")]
            tests += [("secondary", c, p) for c, p in full]
            if stress.search(t.trait):
                tests += [("primary_growth", c, p) for c, p in full]
        else:
            for c in cis.loc[cis["orf"] == t.orf, "centromere"]:
                tests += [("primary_cis", c, p) for p in CEN_PREDICTORS]
            if t.orf in kt_orfs:
                tests += [("primary_kinetochore", c, "non_modal") for c in cens16]
                tests += [("primary_kinetochore", "burden", "rare_burden")]
        if tests:
            plan[t.trait] = tests
    print(f"traits with tests: {len(plan)} ({sum(len(v) for v in plan.values())} tests)")

    # --- run
    rows, samples = [], []
    tinfo = traits.set_index("trait")
    for k, (trait, tests) in enumerate(plan.items(), 1):
        y_raw = read_phen(tinfo.at[trait, "path"])
        strains = tuple(sorted(set(y_raw.index) & complete & set(model.idx)))
        if len(strains) < cfg["min_n"]:
            continue
        y = inverse_normal(y_raw.reindex(list(strains)).values)
        nulls = {}
        for fam, cen, p in tests:
            key = "all" if cen == "burden" else cen_chrom[cen]
            if key not in nulls:
                s, u = model.eig(key, strains)
                nulls[key] = fit_null(y, covariates(list(strains), key), s, u)
            x = (burden.reindex(list(strains)) if cen == "burden"
                 else wide[p][cen].reindex(list(strains))).values.astype(float)
            n_minor = minor_count(x, p)
            res = {"beta": np.nan, "se": np.nan, "p": np.nan}
            if n_minor >= cfg["min_minor"] or p == "rare_burden":
                res = lmm_test(nulls[key], x)
            rows.append({"family": fam, "trait": trait, "klass": tinfo.at[trait, "klass"],
                         "centromere": cen, "predictor": p, "n": len(strains), "n_minor": n_minor,
                         "h2": nulls[key]["h2"], **res})
        samples.append({"trait": trait, "klass": tinfo.at[trait, "klass"], "n": len(strains),
                        "families": ",".join(sorted({f for f, _, _ in tests}))})
        if k % 50 == 0:
            print(f"  {k}/{len(plan)} traits ({time.time() - t0:.0f} s)")
    res = pd.DataFrame(rows)
    save_tsv(pd.DataFrame(samples), out / "trait_samples.tsv")

    primary = res[res["family"].str.startswith("primary")].copy()
    primary["q"] = bh(primary["p"])
    secondary = res[res["family"] == "secondary"].copy()
    secondary["q"] = bh(secondary["p"])
    # --- conditional on the local haplotype: top tests refitted with PCs of SNPs near the CEN
    local_idx = {}
    for c in cens.itertuples():
        ch = cen_chrom[c.centromere]
        near = ((bim["chrom"] == ch) & (bim["pos"] > c.start - cfg["local_window"])
                & (bim["pos"] <= c.end + cfg["local_window"])).values
        local_idx[c.centromere] = np.where(near)[0]
    kin_pos = model.idx

    def conditional(row):
        if row["centromere"] == "burden" or not np.isfinite(row["p"]):
            return pd.Series({"n_local_snps": np.nan, "beta_local": np.nan, "p_local": np.nan})
        y_raw = read_phen(tinfo.at[row["trait"], "path"])
        st = tuple(sorted(set(y_raw.index) & complete & set(kin_pos)))
        key = cen_chrom[row["centromere"]]
        ii = [kin_pos[s] for s in st]
        snps = local_idx[row["centromere"]]
        pcs = local_pcs(geno[np.ix_(ii, snps)], cfg["local_pcs"])
        s, u = model.eig(key, st)
        nl = fit_null(inverse_normal(y_raw.reindex(list(st)).values),
                      np.column_stack([covariates(list(st), key), pcs]), s, u)
        r = lmm_test(nl, wide[row["predictor"]][row["centromere"]].reindex(list(st)).values.astype(float))
        return pd.Series({"n_local_snps": len(snps), "beta_local": r["beta"], "p_local": r["p"]})

    for d in (primary, secondary):
        top = d["p"].nsmallest(cfg["n_conditional"]).index
        d.loc[top, ["n_local_snps", "beta_local", "p_local"]] = d.loc[top].apply(conditional, axis=1).values

    for name, d in (("primary", primary), ("secondary", secondary)):
        lam = {f: genomic_lambda(g["p"]) for f, g in d.groupby("family")}
        d["lambda_family"] = d["family"].map(lam)
        save_tsv(d.sort_values("p"), out / f"assoc_{name}.tsv")
        print(f"{name}: {d['p'].notna().sum()} tests, lambda {lam}, "
              f"min p {d['p'].min():.2g}, q<0.05: {(d['q'] < 0.05).sum()}")

    # --- positive controls (full kinship; CNV markers are not in K)
    ctrl = []
    def control(label, trait_path, marker):
        yr = read_phen(trait_path)
        st = tuple(sorted(set(yr.index) & set(model.idx)))
        x = cnv.loc[list(st), marker].astype(float)
        ok = x >= 0
        st = tuple(np.array(st)[ok.values])
        s, u = model.eig("all", st)
        nl = fit_null(inverse_normal(yr.reindex(list(st)).values), covariates(list(st), "all"), s, u)
        r = lmm_test(nl, x[ok].values)
        ctrl.append({"control": label, "trait": Path(trait_path).name[:-5], "marker": marker,
                     "n": len(st), "h2": nl["h2"], **r})
    for t in traits[traits["trait"].str.contains(cfg["copper_regex"])].itertuples():
        control("CUP1 CNV vs copper growth", t.path, cfg["cup1_marker"])
    expr = traits[traits["klass"] == "expression"].drop_duplicates("orf").set_index("orf")
    cands = [m for m in cnv.columns if m in expr.index]
    var = cnv[cands].replace(-1, np.nan).var().sort_values(ascending=False)
    for m in var.index[:cfg["n_cis_cnv"]]:
        control("ORF CNV vs own expression", expr.at[m, "path"], m)
    save_tsv(pd.DataFrame(ctrl), out / "positive_controls.tsv")
    print(pd.DataFrame(ctrl)[["control", "trait", "n", "beta", "p"]].to_string(index=False))
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
