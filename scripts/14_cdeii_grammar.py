#!/usr/bin/env python3
"""Analysis F2: CDEII A-tract grammar and predicted DNA shape.

Spec docs/specs/2026-11-15-cdeii-grammar-shape.md. Needs scripts/06_microevolution.py outputs
(helsen_cdeii_steps.tsv) and DNAshapeR (scripts/dnashape.R). Outputs in
`<output.dir>/cdeii_grammar/`:
  atracts.tsv             per CDEII (S288C, Helsen modal, ScRAP): A-tract stats + null percentiles
  grammar_null.tsv        per set x CEN x statistic x null: obs vs null; centromere "all" = Wilcoxon
  shape_profiles.tsv      MGW/Roll/HelT/ProT along CDEII, Helsen modal alleles and S288C
  shape_null.tsv          per set x CEN shape summaries vs the two nulls (+ pooled Wilcoxon)
  shape_variants.tsv      per substitution-only allele: shape distance to modal, null percentile
  shape_variants_summary.tsv  per set x CEN x feature: median percentile, one-sided Wilcoxon,
                              calibration KS (null mutant 0 as pseudo-observed)
  shape_events.tsv        per unique Helsen substitution event (CEN, pos, b1>b2): shape change vs
                          every other same-base / same-trinucleotide position (exhaustive null)
  shape_events_summary.tsv    per null x CEN x feature: median percentile, Wilcoxon, Spearman vs frequency
"""
import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kstest, spearmanr, wilcoxon

from centromere.cde import CDEI_LEN, CDEIII_LEN, read_reference
from centromere.grammar import (SHAPES, GenomeWindows, empirical_p, event_null_positions,
                                grammar_stats, mutate, percentile, run_dnashape, shape_distance, shape_long, shape_summary,
                                spectrum_null, substitutions)
from centromere.io import load_config, save_tsv
from centromere.locate import read_fasta
from centromere.ltr_cde import dinuc_shuffle
from centromere.microevo import helsen_cdeii, reference_alleles
from centromere.validate import HELSEN_EXTRA_3P, load_helsen_sd4

STATS = ["n_tracts", "mean_tract", "max_tract", "frac_in_tracts", "phase", "power"]
SUMMARIES = ["mean_mgw", "frac_narrow", "mean_prot"]
NULLS = ["shuffle", "genomic"]


def load_inputs(config, out):
    s288c = {c: s[CDEI_LEN:len(s) - CDEIII_LEN] for c, s in read_reference(config["reference"]).items()}
    scrap = pd.read_csv(out / "preprocessed.tsv", sep="\t").dropna(subset=["cdeii_seq"])
    helsen = load_helsen_sd4(config["validation"]["helsen_sd4"])
    alleles = helsen_cdeii(helsen, HELSEN_EXTRA_3P)
    steps = pd.read_csv(out / "helsen_cdeii_steps.tsv", sep="\t")[["allele", "centromere", "divergent"]]
    alleles = alleles.merge(steps, on=["allele", "centromere"], how="left", validate="1:1")
    if alleles["divergent"].isna().any():
        raise RuntimeError("helsen_cdeii_steps.tsv does not match helsen_cdeii(); rerun 06")
    alleles = alleles[~alleles["divergent"].astype(bool)].reset_index(drop=True)
    modal = reference_alleles(alleles, weight="n_assemblies").set_index("centromere")["ref_seq"].to_dict()
    return s288c, scrap, alleles, modal


def genome_windows(gc, config):
    genome = read_fasta(Path(config["kinetochore"]["s288c_genome"]).expanduser())
    bed = pd.read_csv(config["transposons"]["s288c_cen"], sep="\t", header=None)
    ex = gc["cen_exclusion"]
    blocked = defaultdict(list)
    for _, r in bed.iterrows():
        blocked[r[0]].append((r[1] - ex, r[2] + ex))
    genome = {c: s for c, s in genome.items() if c != "chrM"}
    return GenomeWindows(genome, blocked)


def grammar_part(seqsets, gw, gc, rng):
    kw = dict(min_len=gc["min_tract"], period=gc["period"], tol=gc["phase_tol"], max_d=gc["phase_max_d"])
    uniq = list(dict.fromkeys(seqsets["seq"]))
    obs = {s: grammar_stats(s, **kw) for s in uniq}
    null_seqs = {s: {"shuffle": [dinuc_shuffle(s, rng) for _ in range(gc["n_null"])]} for s in uniq}
    tol_used = {}
    for L in sorted({len(s) for s in uniq}):
        for s in [u for u in uniq if len(u) == L]:
            w, t = gw.sample(L, obs[s]["at"], gc["at_tol"], gc["n_null"], rng,
                             tol_wide=gc["at_tol_wide"], min_n=gc["min_windows"])
            null_seqs[s]["genomic"], tol_used[s] = w, t
    null_stats = {s: {n: pd.DataFrame([grammar_stats(x, **kw) for x in null_seqs[s][n]]) for n in NULLS}
                  for s in uniq}

    rows = []
    for _, r in seqsets.iterrows():
        row = r.drop("seq").to_dict() | obs[r["seq"]]
        row.update(n_genomic=len(null_seqs[r["seq"]]["genomic"]), at_tol=tol_used[r["seq"]])
        for n in NULLS:
            for st in ("phase", "power", "max_tract"):
                row[f"pct_{st}_{n}"] = percentile(null_stats[r["seq"]][n][st].to_numpy(), obs[r["seq"]][st])
        rows.append(row)
    atr = pd.DataFrame(rows)

    rows = []
    for sname in ("modal", "S288C"):
        sub = seqsets[seqsets["set"] == sname]
        for n in NULLS:
            for st in STATS:
                diffs = []
                for _, r in sub.iterrows():
                    nv = null_stats[r["seq"]][n][st].to_numpy(float)
                    o = obs[r["seq"]][st]
                    med = float(np.nanmedian(nv))
                    rows.append({"set": sname, "centromere": r["centromere"], "null": n, "statistic": st,
                                 "obs": o, "null_median": med, "null_q05": float(np.nanquantile(nv, 0.05)),
                                 "null_q95": float(np.nanquantile(nv, 0.95)), "p": empirical_p(nv, o)})
                    diffs.append(o - med)
                d = np.array(diffs, float)
                d = d[~np.isnan(d)]
                p = wilcoxon(d).pvalue if np.any(d != 0) else np.nan
                rows.append({"set": sname, "centromere": "all", "null": n, "statistic": st,
                             "obs": float(np.median(d)), "null_median": 0.0, "n_higher": int(np.sum(d > 0)),
                             "n_cen": len(d), "p": p})
    return atr, pd.DataFrame(rows), null_seqs


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    gc = config["grammar"]
    out = Path(config["output"]["dir"])
    gdir = out / "cdeii_grammar"
    gdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(gc["seed"])

    s288c, scrap, alleles, modal = load_inputs(config, out)
    cens = sorted(modal, key=lambda c: int(c[3:]))
    seqsets = pd.concat([
        pd.DataFrame({"set": "S288C", "centromere": cens, "id": "S288C", "seq": [s288c[c] for c in cens]}),
        pd.DataFrame({"set": "modal", "centromere": cens, "id": "modal", "seq": [modal[c] for c in cens]}),
        pd.DataFrame({"set": "ScRAP", "centromere": scrap["centromere"], "id": scrap["strain"],
                      "seq": scrap["cdeii_seq"]})], ignore_index=True)

    # ---- F2.1 grammar vs composition
    gw = genome_windows(gc, config)
    atr, gnull, null_seqs = grammar_part(seqsets, gw, gc, rng)
    save_tsv(atr, gdir / "atracts.tsv")
    save_tsv(gnull, gdir / "grammar_null.tsv")
    wide = atr[atr["at_tol"] > gc["at_tol"]]
    print(f"F2.1: {len(atr)} CDEIIs; AT tolerance widened for {len(wide)}; "
          f"min genomic windows {atr['n_genomic'].min()}")
    pooled = gnull[gnull["centromere"] == "all"]
    print("Across the 16 CENs (median obs - null median, CENs higher, Wilcoxon p):")
    print(pooled[["set", "null", "statistic", "obs", "n_higher", "n_cen", "p"]].to_string(index=False))

    # ---- F2.3 substitution-only alleles
    var_rows, shape_in = [], {}
    events = defaultdict(lambda: {"n_alleles": 0, "n_assemblies": 0})
    nn = gc["n_null"]
    for sname, df, idcol in (("Helsen", alleles, "allele"), ("ScRAP", scrap, "strain")):
        for _, r in df.iterrows():
            ref, alt = modal[r["centromere"]], r["cdeii_seq"]
            if len(ref) != len(alt) or alt == ref:
                continue
            subs = substitutions(ref, alt)
            if len(subs) > gc["max_subs"]:
                continue
            key = f"{sname}|{r[idcol]}|{r['centromere']}"
            shape_in[f"var|{key}"] = alt
            for i, m in enumerate(spectrum_null(ref, subs, nn + 1, rng)):
                shape_in[f"vnull|{key}|{i}"] = m
            var_rows.append({"set": sname, "centromere": r["centromere"], "id": r[idcol],
                             "n_assemblies": r.get("n_assemblies", 1), "n_subs": len(subs), "key": key})
            if sname == "Helsen":
                for p, b1, b2 in subs:
                    e = events[(r["centromere"], p, b1, b2)]
                    e["n_alleles"] += 1
                    e["n_assemblies"] += int(r["n_assemblies"])
    print(f"\nF2.3: substitution-only alleles (1-{gc['max_subs']} subs vs modal): "
          + ", ".join(f"{k} {v}" for k, v in pd.DataFrame(var_rows)["set"].value_counts().items())
          + f"; unique Helsen substitution events {len(events)}")
    for c in cens:  # every single mutant of each modal allele, for the per-event test
        for p in range(len(modal[c])):
            for b in "ACGT":
                if b != modal[c][p]:
                    shape_in[f"single|{c}|{p}|{b}"] = mutate(modal[c], p, b)

    # ---- shape for everything in one DNAshapeR call
    for c in cens:
        shape_in[f"modal|{c}"] = modal[c]
        shape_in[f"S288C|{c}"] = s288c[c]
        for n in NULLS:
            for i, s in enumerate(null_seqs[modal[c]][n]):
                shape_in[f"null|modal|{c}|{n}|{i}"] = s
            for i, s in enumerate(null_seqs[s288c[c]][n]):
                shape_in[f"null|S288C|{c}|{n}|{i}"] = s
    print(f"DNAshapeR on {len(set(shape_in.values()))} unique sequences ...", flush=True)
    sh = run_dnashape(shape_in, gdir / "shape", gc.get("rscript", "Rscript"))

    # F2.2 profiles and summaries vs nulls
    prof = shape_long(sh, [f"{s}|{c}" for s in ("modal", "S288C") for c in cens])
    prof[["set", "centromere"]] = prof["id"].str.split("|", expand=True)
    save_tsv(prof.drop(columns="id"), gdir / "shape_profiles.tsv")
    rows = []
    for sname in ("modal", "S288C"):
        for n in NULLS:
            diffs = defaultdict(list)
            for c in cens:
                o = shape_summary(sh["MGW"][f"{sname}|{c}"], sh["ProT"][f"{sname}|{c}"])
                src = s288c[c] if sname == "S288C" else modal[c]
                keys = [f"null|{sname}|{c}|{n}|{i}" for i in range(len(null_seqs[src][n]))]
                nul = pd.DataFrame([shape_summary(sh["MGW"][k], sh["ProT"][k]) for k in keys])
                for st in SUMMARIES:
                    med = float(nul[st].median())
                    rows.append({"set": sname, "centromere": c, "null": n, "statistic": st, "obs": o[st],
                                 "null_median": med, "p": empirical_p(nul[st].to_numpy(), o[st])})
                    diffs[st].append(o[st] - med)
            for st in SUMMARIES:
                d = np.array(diffs[st])
                rows.append({"set": sname, "centromere": "all", "null": n, "statistic": st,
                             "obs": float(np.median(d)), "null_median": 0.0, "n_higher": int(np.sum(d > 0)),
                             "n_cen": len(d), "p": wilcoxon(d).pvalue})
    snull = pd.DataFrame(rows)
    save_tsv(snull, gdir / "shape_null.tsv")
    print("\nF2.2 shape summaries, 16 modal CDEIIs vs nulls:")
    print(snull[(snull["centromere"] == "all") & (snull["set"] == "modal")]
          [["null", "statistic", "obs", "n_higher", "p"]].to_string(index=False))

    # F2.3 distances and percentiles
    for v in var_rows:
        c = v["centromere"]
        for f in SHAPES:
            ref = sh[f][f"modal|{c}"]
            d_obs = shape_distance(ref, sh[f][f"var|{v['key']}"])
            d_null = np.array([shape_distance(ref, sh[f][f"vnull|{v['key']}|{i}"]) for i in range(nn + 1)])
            v[f"d_{f}"] = d_obs
            v[f"null_median_{f}"] = float(np.median(d_null[1:]))
            v[f"pct_{f}"] = percentile(d_null[1:], d_obs)
            v[f"calib_{f}"] = percentile(d_null[1:], d_null[0])
    var = pd.DataFrame(var_rows).drop(columns="key")
    save_tsv(var, gdir / "shape_variants.tsv")
    rows = []
    for sname, g in var.groupby("set"):
        for c, gg in [("all", g)] + list(g.groupby("centromere")):
            for f in SHAPES:
                x = gg[f"pct_{f}"].dropna().to_numpy() - 0.5
                row = {"set": sname, "centromere": c, "feature": f, "n": len(x),
                       "median_pct": float(np.median(x + 0.5)) if len(x) else np.nan,
                       "frac_below_half": float(np.mean(x < 0)) if len(x) else np.nan,
                       "p_less": wilcoxon(x, alternative="less").pvalue if len(x) >= 5 and np.any(x) else np.nan}
                if c == "all":
                    row["calib_ks_p"] = kstest(gg[f"calib_{f}"].dropna(), "uniform").pvalue
                rows.append(row)
    vs = pd.DataFrame(rows)
    save_tsv(vs, gdir / "shape_variants_summary.tsv")
    # F2.3b unique substitution events, exhaustive single-mutant nulls
    erows = []
    for (c, p, b1, b2), e in events.items():
        ref = {f: sh[f][f"modal|{c}"] for f in SHAPES}
        d = {f: (lambda q, f=f: shape_distance(ref[f], sh[f][f"single|{c}|{q}|{b2}"])) for f in SHAPES}
        row = {"centromere": c, "pos": p + 1, "b1": b1, "b2": b2, "context": modal[c][max(p - 1, 0):p + 2]} | e
        for nname, ctx in (("base", False), ("trinuc", True)):
            qs = event_null_positions(modal[c], p, ctx)
            row[f"n_null_{nname}"] = len(qs)
            for f in SHAPES:
                row[f"pct_{f}_{nname}"] = percentile(np.array([d[f](q) for q in qs]), d[f](p)) \
                    if len(qs) >= gc["min_event_null"] else np.nan
        for f in SHAPES:
            row[f"d_{f}"] = d[f](p)
        erows.append(row)
    ev = pd.DataFrame(erows).sort_values(["centromere", "pos", "b2"])
    save_tsv(ev, gdir / "shape_events.tsv")
    rows = []
    for nname in ("base", "trinuc"):
        for c, g in [("all", ev)] + list(ev.groupby("centromere")):
            for f in SHAPES:
                col = f"pct_{f}_{nname}"
                gg = g.dropna(subset=[col])
                x = gg[col].to_numpy() - 0.5
                row = {"null": nname, "centromere": c, "feature": f, "n_events": len(x),
                       "median_pct": float(np.median(x + 0.5)) if len(x) else np.nan,
                       "p_less": wilcoxon(x, alternative="less").pvalue if len(x) >= 5 and np.any(x) else np.nan}
                if c == "all":
                    rho, pr = spearmanr(gg[col], gg["n_assemblies"])
                    row.update(rho_freq=rho, p_rho=pr)
                rows.append(row)
    es = pd.DataFrame(rows)
    save_tsv(es, gdir / "shape_events_summary.tsv")
    ea = es[es["centromere"] == "all"]
    print("\nF2.3b unique substitution events vs every other position (same base / same trinucleotide):")
    print(ea.to_string(index=False))
    for nname in ("base", "trinuc"):
        x = es[(es["centromere"] != "all") & (es["null"] == nname)]
        print(nname, "CENs with median < 0.5:",
              x.groupby("feature")["median_pct"].apply(lambda s: f"{int((s < 0.5).sum())}/{s.notna().sum()}").to_dict())

    print("\nF2.3 natural substitutions vs spectrum-matched random positions "
          "(percentile < 0.5 = less shape change than random):")
    print(vs[vs["centromere"] == "all"].to_string(index=False))


if __name__ == "__main__":
    main()
