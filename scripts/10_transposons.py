#!/usr/bin/env python3
"""Analysis D: Ty elements and solo LTRs around centromeres (10 assemblies + S288C).

Outputs in `<output.dir>/` (spec docs/specs/2026-11-13-transposons-near-cen.md):
  te_elements.tsv      genome-wide Ty/LTR calls per genome
  trna_loci.tsv        tRNA genes per genome
  te_near_cen.tsv      elements within +-20 kb of each CEN (trna_near_cen.tsv: same for tRNAs)
  te_cen_summary.tsv   genome x CEN: counts in +-1/5/20 kb, nearest element
  te_null.tsv          observed vs tRNA-matched and uniform null, per genome and pooled
  te_sites.tsv         orthologous insertion sites near CENs x strain (present/absent/...)
  te_flank_calls.tsv   raw junction-flank calls behind te_sites.tsv
  te_s288c_check.tsv   SGD TEs near S288C CENs and whether they were recovered
  te_loegler_s5.tsv    Loegler 2025 Table S5 hotspots within +-20 kb of a CEN
"""
import argparse
import hashlib
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from centromere.genes import make_db
from centromere.io import load_config, save_tsv, tool
from centromere.transposons import (HIT_COLS, UnionFind, call_elements, call_trnas, candidate_positions,
                                    empirical_p, flank_call, interval_distance, matched_pools, near_cen,
                                    nearest_trna, parse_hits, replicate_sums, te_clusters, unique_best,
                                    window_counts)

ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
         "XI", "XII", "XIII", "XIV", "XV", "XVI"]


def genomes(config: dict) -> list[tuple[str, Path]]:
    samples = pd.read_csv(config["input"]["sample_map"], sep="\t")
    gdir = Path(config["cde"]["genome_dir"])
    out = [(r["strain"], gdir / f"{r['accession']}.fasta.gz") for _, r in samples.iterrows()]
    return out + [("S288C", Path(config["kinetochore"]["s288c_genome"]).expanduser())]


def cen_table(config: dict) -> pd.DataFrame:
    f = pd.read_csv(Path(config["output"]["dir"]) / "filtered.tsv", sep="\t")
    cens = f[["strain", "centromere", "sseqid", "cen_start", "cen_end", "strand"]].rename(
        columns={"strain": "genome", "sseqid": "contig"})
    bed = pd.read_csv(config["transposons"]["s288c_cen"], sep="\t", header=None,
                      names=["contig", "s0", "cen_end", "centromere", "score", "strand"])
    s288c = bed.assign(genome="S288C", cen_start=bed["s0"] + 1)[cens.columns]
    return pd.concat([cens, s288c], ignore_index=True)


def blastn(query: Path, db: Path, blast_bin: Path, out: Path, task: str, evalue: float) -> pd.DataFrame:
    if not out.exists():
        cmd = [tool(blast_bin, "blastn"), "-task", task, "-query", str(query), "-db", str(db),
               "-evalue", str(evalue), "-dust", "no", "-max_target_seqs", "5000",
               "-num_threads", "4", "-outfmt", "6 " + " ".join(HIT_COLS), "-out", str(out)]
        subprocess.run(cmd, check=True)
    df = pd.read_csv(out, sep="\t", header=None, names=HIT_COLS)
    return parse_hits(df)


def s288c_check(elements: pd.DataFrame, cens: pd.DataFrame, sgd: pd.DataFrame) -> pd.DataFrame:
    """SGD Ty/LTRs within 20 kb of an S288C CEN: recovered if >= 50 % overlapped by our calls."""
    s_el = elements[elements["genome"] == "S288C"]
    rows = []
    for _, c in cens[cens["genome"] == "S288C"].iterrows():
        t = sgd[sgd["chrom"] == c["contig"]]
        d = interval_distance(t["start"], t["end"], c["cen_start"], c["cen_end"])
        for _, x in t[d <= 20000].iterrows():
            e = s_el[s_el["contig"] == x["chrom"]]
            ov = (np.minimum(e["end"], x["end"]) - np.maximum(e["start"], x["start"]) + 1).clip(lower=0)
            L = x["end"] - x["start"] + 1
            rows.append({"centromere": c["centromere"], "name": x["name"], "type": x["type"],
                         "family": x["family"], "start": x["start"], "end": x["end"], "length": L,
                         "recovered": bool((ov >= 0.5 * L).any()),
                         "our_class": ",".join(e.loc[ov > 0, "class"])})
    return pd.DataFrame(rows)


def genome_recall(elements: pd.DataFrame, sgd: pd.DataFrame) -> pd.DataFrame:
    s_el = elements[elements["genome"] == "S288C"]
    rec = []
    for _, x in sgd.iterrows():
        e = s_el[s_el["contig"] == x["chrom"]]
        ov = (np.minimum(e["end"], x["end"]) - np.maximum(e["start"], x["start"]) + 1).clip(lower=0)
        rec.append(bool((ov >= 0.5 * (x["end"] - x["start"] + 1)).any()))
    sgd = sgd.assign(recovered=rec, length=sgd["end"] - sgd["start"] + 1)
    sgd["size"] = np.where(sgd["type"] == "LTR_retrotransposon", "Ty element",
                           np.where(sgd["length"] >= 0.9 * sgd.groupby("family")["length"].transform(
                               lambda x: x.round(-1).mode().iloc[0]), "full-length LTR", "LTR fragment"))
    return sgd.groupby("size")["recovered"].agg(["sum", "count"])


def null_test(elements, trnas, cens, contig_len, tcfg, rng) -> pd.DataFrame:
    windows = tcfg["windows"]
    rows, pooled = [], {"tRNA-matched": [], "uniform": []}
    obs_pool = np.zeros(len(windows))
    for g, cg in cens.groupby("genome"):
        el, tr = elements[elements["genome"] == g], trnas[trnas["genome"] == g]
        obs = np.zeros(len(windows))
        cen_d, cen_len = [], int(np.median(cg["cen_end"] - cg["cen_start"] + 1))
        for _, c in cg.iterrows():
            obs += window_counts(el, c["contig"], [c["cen_start"]], [c["cen_end"]], windows)[0]
            cen_d.append(nearest_trna(tr, c["contig"], [c["cen_start"]], [c["cen_end"]])[0])
        cand = candidate_positions(contig_len[g], cg, tcfg["null_candidates"], cen_len,
                                   max(windows), rng)
        cand["d_trna"] = np.nan
        counts = np.zeros((len(cand), len(windows)), dtype=int)
        for contig, idx in cand.groupby("contig").groups.items():
            lo, hi = cand.loc[idx, "lo"].to_numpy(), cand.loc[idx, "hi"].to_numpy()
            cand.loc[idx, "d_trna"] = nearest_trna(tr, contig, lo, hi)
            counts[np.asarray(idx)] = window_counts(el, contig, lo, hi, windows)
        per_cen = int(np.ceil(tcfg["null_per_genome"] / len(cg)))
        pools = matched_pools(cand, np.array(cen_d), tcfg["null_tolerance"], per_cen, rng)
        uniform = [rng.choice(len(cand), size=per_cen, replace=False) for _ in cen_d]
        for name, pl in (("tRNA-matched", pools), ("uniform", uniform)):
            assert all(len(p) >= 20 for p in pl), f"{g}: null pool too small ({[len(p) for p in pl]})"
            reps = replicate_sums([counts[p] for p in pl], tcfg["null_replicates"], rng)
            pooled[name].append(reps)
            for j, w in enumerate(windows):
                rows.append({"genome": g, "null": name, "window_kb": w // 1000, "observed": obs[j],
                             "null_mean": reps[:, j].mean(), "null_q05": np.quantile(reps[:, j], 0.05),
                             "null_q95": np.quantile(reps[:, j], 0.95),
                             "ratio": obs[j] / reps[:, j].mean() if reps[:, j].mean() else np.nan,
                             "p": empirical_p(reps[:, j], obs[j]),
                             "pool_min": min(len(p) for p in pl), "independent": True})
        obs_pool += obs
    per = pd.DataFrame(rows)
    for name, reps in pooled.items():
        # pooled rows are NOT independent: insertions shared by descent are counted once per genome
        tot = np.sum(reps, axis=0)
        for j, w in enumerate(windows):
            pg = per[(per["null"] == name) & (per["window_kb"] == w // 1000)]
            rows.append({"genome": "all", "null": name, "window_kb": w // 1000, "observed": obs_pool[j],
                         "genomes_above_null": int((pg["observed"] > pg["null_mean"]).sum()),
                         "null_mean": tot[:, j].mean(), "null_q05": np.quantile(tot[:, j], 0.05),
                         "null_q95": np.quantile(tot[:, j], 0.95),
                         "ratio": obs_pool[j] / tot[:, j].mean() if tot[:, j].mean() else np.nan,
                         "p": empirical_p(tot[:, j], obs_pool[j]), "pool_min": np.nan,
                         "independent": False})
    return pd.DataFrame(rows)


def flank_seqs(clusters: pd.DataFrame, contigs: dict[str, str], bp: int) -> list[tuple[str, str]]:
    out = []
    for _, c in clusters.iterrows():
        s = contigs[c["contig"]]
        if c["start"] - bp < 1 or c["end"] + bp > len(s):
            continue
        out.append((f"{c['cluster_id']}|L", s[c["start"] - 1 - bp:c["start"] - 1]))
        out.append((f"{c['cluster_id']}|R", s[c["end"]:c["end"] + bp]))
    return out


def polymorphic_sites(elements, near, contigs, dbs, blast_bin, work, tcfg):
    """Junction-flank orthology for TE clusters within 20 kb of a CEN."""
    clusters = []
    for g, el in elements.groupby("genome"):
        cl = te_clusters(el, tcfg["flank_bp"]).assign(genome=g)
        cl["cluster_id"] = [f"{g}:c{i}" for i in range(len(cl))]
        clusters.append(cl)
    clusters = pd.concat(clusters, ignore_index=True)
    near_ids = set(near["element_id"])
    has_near = clusters["element_ids"].map(lambda s: any(x in near_ids for x in s.split(",")))
    eid_near = near.sort_values("distance").drop_duplicates("element_id").set_index("element_id")
    eid2cen = eid_near["centromere"]
    src = clusters[has_near].copy()
    src["centromere"] = src["element_ids"].map(
        lambda s: ",".join(sorted({eid2cen[x] for x in s.split(",") if x in eid2cen})))

    # One query with every source genome's junction flanks, one BLAST per target genome (n runs,
    # not n^2); the hits of a genome's own flanks in its own assembly are not used. The cache
    # name carries the query's md5, so a run on another genome set never reuses old hits.
    lines = [f">{name}\n{seq}\n" for g, sg in src.groupby("genome")
             for name, seq in flank_seqs(sg, contigs[g], tcfg["flank_bp"])]
    tag = hashlib.md5("".join(lines).encode()).hexdigest()[:10]
    q = work / f"flanks_all_{tag}.fa"
    q.write_text("".join(lines))
    by_genome = {g: sg for g, sg in src.groupby("genome")}
    tcls = {t: c.rename(columns={"cluster_id": "element_id"}) for t, c in clusters.groupby("genome")}
    empty = clusters.iloc[:0].rename(columns={"cluster_id": "element_id"})
    calls = []
    for t, db in dbs.items():
        hits = blastn(q, db, blast_bin, work / f"flanks_{tag}_vs_{t}.tsv", "megablast", 1e-20)
        best = unique_best(hits, 90, 0.8)
        best = dict(zip(best["qseqid"], (r for _, r in best.iterrows())))
        tcl = tcls.get(t, empty)
        for g, sg in by_genome.items():
            if g == t:
                continue
            for c in sg[["cluster_id", "centromere"]].itertuples(index=False):
                call, gap, tid = flank_call(best.get(f"{c.cluster_id}|L"), best.get(f"{c.cluster_id}|R"), tcl)
                calls.append((c.cluster_id, g, c.centromere, t, call, gap, tid))
    calls = pd.DataFrame(calls, columns=["source", "source_genome", "centromere", "target_genome",
                                         "call", "gap", "target_cluster"])

    uf = UnionFind()
    for cid in src["cluster_id"]:
        uf.find(cid)
    for r in calls[calls["call"] == "present"].itertuples(index=False):
        for tid in r.target_cluster.split(","):
            uf.union(r.source, tid)
    genome_of = dict(zip(clusters["cluster_id"], clusters["genome"]))
    label = dict(zip(clusters["cluster_id"], clusters["families"] + "/" + clusters["classes"]))
    src["site"] = src["cluster_id"].map(uf.find)
    calls["site"] = calls["source"].map(uf.find)
    site_members: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for m in list(uf.parent):
        site_members[uf.find(m)][genome_of[m]].append(m)
    # majority non-present call per site and target genome (ties: first seen, as value_counts)
    np_calls = calls[calls["call"] != "present"].reset_index(drop=True)
    cnt = (np_calls.assign(_i=np_calls.index).groupby(["site", "target_genome", "call"], sort=False)
           .agg(n=("_i", "size"), first=("_i", "min")).reset_index()
           .sort_values(["n", "first"], ascending=[False, True])
           .drop_duplicates(["site", "target_genome"]))
    majority = dict(zip(zip(cnt["site"], cnt["target_genome"]), cnt["call"]))
    strains = sorted(elements["genome"].unique())
    rows = []
    for site, members in src.groupby("site"):
        ids = [x for m in members["element_ids"] for x in m.split(",") if x in eid_near.index]
        best = eid_near.loc[ids].sort_values("distance").iloc[0]
        row = {"site": site, "centromere": ",".join(sorted(set(",".join(members["centromere"]).split(",")))),
               "distance": int(best["distance"]), "side": best["side"]}
        mem = site_members[site]
        for s in strains:
            row[s] = ("present:" + ";".join(label[m] for m in sorted(mem[s]))) if mem.get(s) \
                else majority.get((site, s), "unresolved")
        states = [row[s].split(":")[0] for s in strains]
        row["n_present"] = states.count("present")
        row["n_absent"] = states.count("absent")
        row["polymorphic"] = row["n_present"] > 0 and row["n_absent"] > 0
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["centromere", "n_present"]), calls


def loegler_s5(path: str, cens: pd.DataFrame, elements: pd.DataFrame) -> pd.DataFrame:
    s5 = pd.read_excel(path, sheet_name="TableS5", header=2)
    s5["contig"] = "chr" + s5["Chromosome"].str.replace("chromosome", "").astype(int).map(lambda i: ROMAN[i - 1])
    rows = []
    for _, c in cens[cens["genome"] == "S288C"].iterrows():
        h = s5[s5["contig"] == c["contig"]]
        d = interval_distance(h["Start"], h["End"], c["cen_start"], c["cen_end"])
        for _, x in h[d <= 20000].iterrows():
            n = window_counts(elements[elements["genome"] == "S288C"], c["contig"],
                              [x["Start"]], [x["End"]], [0])[0, 0]
            rows.append({"centromere": c["centromere"], "region": x["Region"], "start": x["Start"],
                         "end": x["End"], "type": x["Type"], "distance": int(d[h.index.get_loc(x.name)]),
                         "TEs": x["TEs"], "LTR": x["LTR"], "tRNA": x["tRNA"],
                         "our_s288c_elements_in_hotspot": int(n)})
    return pd.DataFrame(rows, columns=["centromere", "region", "start", "end", "type", "distance",
                                       "TEs", "LTR", "tRNA", "our_s288c_elements_in_hotspot"])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    tcfg = config["transposons"]
    out = Path(config["output"]["dir"])
    work = out / "transposons"
    work.mkdir(parents=True, exist_ok=True)
    blast_bin = config["kinetochore"].get("blast_bin")   # null -> PATH
    rng = np.random.default_rng(tcfg["seed"])
    cens = cen_table(config)

    elements, trnas, contigs, contig_len, dbs = [], [], {}, {}, {}
    for g, path in genomes(config):
        db = make_db(path, out / "kinetochore_blastdb", blast_bin)
        dbs[g] = db
        contigs[g] = {r.id: str(r.seq) for r in SeqIO.parse(db, "fasta")}
        contig_len[g] = {k: len(v) for k, v in contigs[g].items()}
        hits = blastn(Path(tcfg["library"]), db, blast_bin, work / f"{g}.te.tsv", "blastn", 1e-10)
        el, dropped = call_elements(hits, tcfg["min_pident"], tcfg["min_qcov"])
        el.insert(0, "genome", g)
        el["element_id"] = [f"{g}:e{i}" for i in range(len(el))]
        elements.append(el)
        tr = call_trnas(blastn(Path(tcfg["trna"]), db, blast_bin, work / f"{g}.trna.tsv", "blastn", 1e-10))
        tr.insert(0, "genome", g)
        trnas.append(tr)
        print(f"{g}: {len(el)} elements {el['class'].value_counts().to_dict()}, "
              f"{dropped} fragments dropped, {len(tr)} tRNAs")
    elements, trnas = pd.concat(elements, ignore_index=True), pd.concat(trnas, ignore_index=True)
    save_tsv(elements, out / "te_elements.tsv")
    save_tsv(trnas, out / "trna_loci.tsv")

    sgd = pd.read_csv(tcfg["sgd_te"], sep="\t")
    chk = s288c_check(elements, cens, sgd)
    save_tsv(chk, out / "te_s288c_check.tsv")
    print(f"\nCheck 1 - SGD TEs within 20 kb of S288C CENs recovered: "
          f"{chk['recovered'].sum()}/{len(chk)}")
    if (~chk["recovered"]).any():
        print(chk[~chk["recovered"]].to_string(index=False))
    print("S288C genome-wide recall of SGD annotation:\n" + genome_recall(elements, sgd).to_string())

    near = near_cen(elements, cens, max(tcfg["windows"]))
    save_tsv(near, out / "te_near_cen.tsv")
    tr_ids = trnas.assign(element_id=[f"t{i}" for i in range(len(trnas))], family="tRNA", **{"class": "tRNA"})
    save_tsv(near_cen(tr_ids, cens, max(tcfg["windows"])), out / "trna_near_cen.tsv")
    summ = []
    for _, c in cens.iterrows():
        el = elements[elements["genome"] == c["genome"]]
        n = window_counts(el, c["contig"], [c["cen_start"]], [c["cen_end"]], tcfg["windows"])[0]
        e = el[el["contig"] == c["contig"]]
        d = interval_distance(e["start"], e["end"], c["cen_start"], c["cen_end"])
        row = {"genome": c["genome"], "centromere": c["centromere"],
               "nearest_te": int(d.min()) if len(d) else None,
               "nearest_trna": nearest_trna(trnas[trnas["genome"] == c["genome"]], c["contig"],
                                            [c["cen_start"]], [c["cen_end"]])[0]}
        row |= {f"n_{w // 1000}kb": int(k) for w, k in zip(tcfg["windows"], n)}
        summ.append(row)
    summ = pd.DataFrame(summ)
    save_tsv(summ, out / "te_cen_summary.tsv")

    nul = null_test(elements, trnas, cens, contig_len, tcfg, rng)
    save_tsv(nul, out / "te_null.tsv")
    print("\nNull (pooled):\n" + nul[nul["genome"] == "all"].round(3).to_string(index=False))

    sites, calls = polymorphic_sites(elements, near, contigs, dbs, blast_bin, work, tcfg)
    save_tsv(sites, out / "te_sites.tsv")
    save_tsv(calls, out / "te_flank_calls.tsv")
    print(f"\nSites near CENs: {len(sites)}, polymorphic: {int(sites['polymorphic'].sum())}; "
          f"flank calls {calls['call'].value_counts().to_dict()}")

    s5 = loegler_s5(tcfg["loegler"], cens, elements)
    save_tsv(s5, out / "te_loegler_s5.tsv")
    print(f"Loegler S5 hotspots within 20 kb of a CEN: {len(s5)}")


if __name__ == "__main__":
    main()
