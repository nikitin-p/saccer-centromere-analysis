#!/usr/bin/env python3
"""L3 liftover + analysis F: nucleotide diversity at CDEI/CDEII/CDEIII vs the +-2 kb flanks.

Spec docs/specs/2026-11-13-flanks-diversity-trees.md, steps 0.2-0.5 and 1. Needs the PAFs
from scripts/run_minimap2.sh. Outputs in `<output.dir>/`:
  l3_loci.tsv                    S288C CEN lifted per assembly vs the L1 (BLAST+CDE) locus
  flanks/<CEN>_windows.fasta     CEN +- flank_bp per assembly + S288C, CDEI on the left
  flanks/<CEN>_windows.aln.fasta MAFFT of the above
  flanks/windows.tsv             per assembly x CEN: lifted window, flank kept on each side
  diversity_windows.tsv          pi, pi_sub, indel_rate per CEN x region (CDEs, 100-bp windows)
  diversity_summary.tsv          per CEN x element: percentile among its flank windows;
                                 pooled rows (centromere = "all"): medians + Wilcoxon
  diversity_crosscheck.tsv       CEN pi from the window MSA vs from aligned.tsv (step 1.6)
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from centromere.align import run_mafft
from centromere.io import load_config, save_tsv
from centromere.locate import assemblies, lift_interval, lift_window, parse_paf, read_fasta, s288c_cens
from centromere.microevo import diversity_windows, leave_one_out_pi, pairwise_pi, region_map, revcomp

ELEMENTS = ["CDEI", "CDEII", "CDEIII"]
CROSSCHECK_MARGIN = 500   # left-flank columns kept for the crosscheck (see cen_columns)


def l3_loci(asm, cens, blocks, l1, bridge):
    """Exact liftover first; a CEN in an alignment gap between collinear blocks is 'bridged'
    (the L3 interval is then the whole unaligned query segment and agree = contains L1)."""
    rows = []
    for _, a in asm.iterrows():
        for _, c in cens.iterrows():
            b = blocks[a["accession"]]
            r = lift_interval(b, c["chrom"], c["start"], c["end"])
            status = "lifted"
            if isinstance(r, str):
                r2 = lift_interval(b, c["chrom"], c["start"], c["end"], bridge)
                status, r = (r, r) if isinstance(r2, str) else ("bridged", r2)
            row = {"strain": a["strain"], "centromere": c["centromere"]}
            k = l1.get((a["strain"], c["centromere"]))
            if isinstance(r, str):
                row.update(status=r)
            else:
                contig, qs, qe, strand = r
                row.update(status=status, contig=contig, l3_start=qs + 1, l3_end=qe, l3_strand=strand)
            if k is not None:
                row.update(l1_contig=k["sseqid"], l1_start=k["cen_start"], l1_end=k["cen_end"])
                if row["status"] in ("lifted", "bridged"):
                    same = row["contig"] == k["sseqid"]
                    ov = min(row["l3_end"], k["cen_end"]) - max(row["l3_start"], k["cen_start"]) + 1
                    contains = row["l3_start"] <= k["cen_start"] and row["l3_end"] >= k["cen_end"]
                    row.update(d_start=row["l3_start"] - k["cen_start"], d_end=row["l3_end"] - k["cen_end"],
                               agree=bool(same and (ov > 0 if row["status"] == "lifted" else contains)))
            rows.append(row)
    return pd.DataFrame(rows)


def windows(asm, cens, blocks, ref, flank, outdir, mafft_args, bridge):
    rows, alns = [], {}
    for _, c in cens.iterrows():
        cen = c["centromere"]
        s288c = ref[c["chrom"]][c["start"] - flank:c["end"] + flank]
        if c["strand"] == "-":
            s288c = revcomp(s288c)
        assert s288c[flank:flank + c["end"] - c["start"]] == c["cen_seq"], f"{cen}: BED != reference FASTA"
        seqs = {"S288C": s288c}
        for _, a in asm.iterrows():
            w = lift_window(blocks[a["accession"]], c["chrom"], c["start"], c["end"], flank, bridge=bridge)
            row = {"centromere": cen, "strain": a["strain"], **w}
            if w["status"] in ("ok", "shortened"):
                s = a["contigs"][w["contig"]][w["qstart"]:w["qend"]]
                if w["strand"] == "-":
                    s = revcomp(s)
                left, right = w["left_kept"], w["right_kept"]
                if c["strand"] == "-":
                    s, left, right = revcomp(s), right, left
                row.update(left_kept=left, right_kept=right, length=len(s))
                seqs[a["strain"]] = s
            rows.append(row)
        fa = outdir / f"{cen}_windows.fasta"
        fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        aln_path = outdir / f"{cen}_windows.aln.fasta"
        if aln_path.exists():
            aln = read_fasta(aln_path)
        else:
            aln = dict(zip(seqs, (x.upper() for x in run_mafft(list(seqs.values()), mafft_args))))
            aln_path.write_text("".join(f">{k}\n{v}\n" for k, v in aln.items()))
        if set(aln) != set(seqs):   # cached MSA from a different strain set
            aln_path.unlink()
            return windows(asm, cens, blocks, ref, flank, outdir, mafft_args, bridge)
        alns[cen] = aln
    return pd.DataFrame(rows), alns


def percentile(x, ref):
    ref = np.asarray(ref, dtype=float)
    ref = ref[~np.isnan(ref)]
    if np.isnan(x) or len(ref) == 0:
        return np.nan
    return 100 * (np.mean(ref < x) + 0.5 * np.mean(ref == x))


def summarise(div):
    rows = []
    for cen, d in div.groupby("centromere", sort=False):
        fl = d[~d["region"].isin(ELEMENTS)]
        for el in ELEMENTS:
            e = d[d["region"] == el].iloc[0]
            row = {"centromere": cen, "element": el, "n_sites": e["n_sites"], "n_flank_windows": int(fl["pi"].notna().sum())}
            for m in ("pi", "pi_sub", "indel_rate"):
                row[m] = e[m]
                row[f"flank_mean_{m}"] = fl[m].mean()
                row[f"percentile_{m}"] = percentile(e[m], fl[m])
            rows.append(row)
    per = pd.DataFrame(rows)
    pooled = []
    for el, g in per.groupby("element", sort=False):
        row = {"centromere": "all", "element": el, "n_cen": int(g["pi"].notna().sum())}
        for m in ("pi", "pi_sub", "indel_rate"):
            row[f"median_{m}"] = g[m].median()
            row[f"median_flank_{m}"] = g[f"flank_mean_{m}"].median()
            row[f"n_below_flank_{m}"] = int((g[m] < g[f"flank_mean_{m}"]).sum())
        ok = g[["pi", "flank_mean_pi"]].dropna()
        t = wilcoxon(ok["pi"], ok["flank_mean_pi"], zero_method="wilcox")   # the one test (pi)
        row.update(wilcoxon_stat=t.statistic, wilcoxon_p=t.pvalue)
        pooled.append(row)
    return pd.concat([per, pd.DataFrame(pooled)], ignore_index=True)


def cen_columns(aln, cols, margin=CROSSCHECK_MARGIN):
    """Window MSA cut to the CEN columns plus `margin` columns on the left, with their group
    labels. The CEN group's pi depends on flank columns only through the indel state at the
    last comparable column before the CEN (right-flank columns never matter), so the cut is
    exact for every pair with a comparable column within `margin`; for the rest a gap run
    entering the CEN counts as opening there, as it does in the CEN-only MSA. Avoids the full
    pass over 4-45 kb wide windows (1086 strains: hours per CEN on the cluster)."""
    idx = np.flatnonzero(cols)
    lo, hi = max(idx[0] - margin, 0), idx[-1] + 1
    return [s[lo:hi] for s in aln.values()], np.where(cols[lo:hi], "CEN", "flank")


def crosscheck(div_alns, cens, aligned):
    """CEN pi from the window MSA vs the CEN-only MSA in aligned.tsv, same sequences. Where they
    differ by > 10 %, `explained_by` = the strain whose removal brings both within 10 %
    (the two MSAs place that strain's indels differently)."""
    def both(w, groups, a, n_cen):
        # site threshold from the S288C CEN length (as diversity_windows), not the column
        # count: one strain's insertion (CEN5: 244 bp) widens the MSA and would make every pair NA
        m = [a["ref_aligned_seq"].iloc[0]] + list(a["aligned_seq"])
        pw = pairwise_pi(w, groups, n_sites={"CEN": n_cen, "flank": 1}).set_index("group").loc["CEN"]
        pa = pairwise_pi(m, np.array(["CEN"] * len(m[0])), n_sites={"CEN": n_cen}).iloc[0]
        return pw, pa

    rows = []
    for _, c in cens.iterrows():
        cen = c["centromere"]
        aln = div_alns[cen]
        rm = region_map(aln["S288C"], aln["S288C"].replace("-", "").find(c["cen_seq"]), len(c["cen_seq"]))
        cols = rm["region"].isin(ELEMENTS).to_numpy()
        wseqs, groups = cen_columns(aln, cols)
        a = aligned[(aligned["centromere"] == cen) & aligned["strain"].isin(aln)]
        n_cen = len(c["cen_seq"])
        w, m = both(wseqs, groups, a, n_cen)
        rel = (w["pi"] - m["pi"]) / m["pi"]
        expl, rel_wo = "", np.nan
        if abs(rel) > 0.10:
            # pi without each strain, from one pass over the pairs per MSA (was 2n full passes)
            _, w_wo = leave_one_out_pi(wseqs, groups, "CEN", n_sites={"CEN": n_cen, "flank": 1})
            _, m_wo = leave_one_out_pi([a["ref_aligned_seq"].iloc[0]] + list(a["aligned_seq"]),
                                       np.array(["CEN"] * len(a["ref_aligned_seq"].iloc[0])), "CEN",
                                       n_sites={"CEN": n_cen})
            widx = {k: n for n, k in enumerate(aln)}
            alt = {k: (w_wo[widx[k]] - m_wo[1 + n]) / m_wo[1 + n] for n, k in enumerate(a["strain"])}
            k = min(alt, key=lambda k: abs(alt[k]))
            if abs(alt[k]) <= 0.10:
                expl, rel_wo = k, alt[k]
        rows.append({"centromere": cen, "n_seqs": len(aln), "pi_window": w["pi"], "pi_aligned": m["pi"],
                     "pi_sub_window": w["pi_sub"], "pi_sub_aligned": m["pi_sub"], "rel_diff": rel,
                     "explained_by": expl, "rel_diff_without": rel_wo})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    cfg = load_config(ap.parse_args().config)
    out = Path(cfg["output"]["dir"])
    flank = cfg["diversity"]["flank_bp"]
    fdir = out / "flanks"
    fdir.mkdir(parents=True, exist_ok=True)

    asm = assemblies(cfg)
    missing = [str(p) for p in asm["paf"] if not Path(p).exists()]
    if missing:
        raise SystemExit(f"run scripts/run_minimap2.sh first; missing {missing}")
    cens = s288c_cens(cfg)
    ref = read_fasta(Path(cfg["locate"]["paf_dir"]) / "sacCer3.nuclear.fa")
    blocks = {a: parse_paf(p, cfg["locate"]["min_mapq"]) for a, p in zip(asm["accession"], asm["paf"])}
    asm["contigs"] = [read_fasta(p) for p in asm["fasta"]]

    f = pd.read_csv(out / "filtered.tsv", sep="\t")
    l1 = {(r["strain"], r["centromere"]): r for _, r in f.iterrows()}
    bridge = cfg["locate"].get("bridge_bp", 0)
    loci = l3_loci(asm, cens, blocks, l1, bridge)
    save_tsv(loci, out / "l3_loci.tsv")
    n_l1 = loci["l1_contig"].notna().sum()
    n_agree = int(loci["agree"].fillna(False).sum())
    print(f"L1 vs L3: {n_agree}/{n_l1} loci agree (same contig, overlapping)")
    lifted = loci[loci["status"] == "lifted"]
    print(f"  status: {loci['status'].value_counts().to_dict()}; exact lifts with |d_start| <= 2: "
          f"{(lifted['d_start'].abs() <= 2).sum()}, |d_end| <= 2: {(lifted['d_end'].abs() <= 2).sum()}")
    br = loci[loci["status"] == "bridged"]
    if len(br):
        print("  bridged (CEN in a gap between collinear minimap2 blocks; L3 = whole gap):")
        print(br[["strain", "centromere", "d_start", "d_end", "agree"]].to_string(index=False))
    bad = loci[~loci["agree"].fillna(False).astype(bool)]
    if len(bad):
        print(bad[["strain", "centromere", "status", "contig", "l1_contig", "d_start", "d_end"]].to_string(index=False))

    win, alns = windows(asm, cens, blocks, ref, flank, fdir, cfg["align"]["mafft_args"], bridge)
    save_tsv(win, fdir / "windows.tsv")
    print(f"windows: {win['status'].value_counts().to_dict()}")

    div = []
    for _, c in cens.iterrows():
        d = diversity_windows(alns[c["centromere"]], "S288C", flank, len(c["cen_seq"]),
                              window=cfg["diversity"]["window"], min_pair_frac=cfg["diversity"]["min_pair_frac"])
        div.append(d.assign(centromere=c["centromere"], n_seqs=len(alns[c["centromere"]])))
    div = pd.concat(div, ignore_index=True)
    div = div[["centromere", "region", "offset", "n_seqs", "n_sites", "n_pairs", "pi", "pi_sub", "indel_rate"]]
    save_tsv(div, out / "diversity_windows.tsv")

    summ = summarise(div)
    save_tsv(summ, out / "diversity_summary.tsv")
    pooled = summ[summ["centromere"] == "all"]
    print(pooled[["element", "n_cen", "median_pi", "median_flank_pi", "n_below_flank_pi", "wilcoxon_p",
                  "median_pi_sub", "median_flank_pi_sub", "median_indel_rate", "median_flank_indel_rate"]]
          .to_string(index=False))

    aligned = pd.read_csv(out / "aligned.tsv", sep="\t")
    cc = crosscheck(alns, cens, aligned)
    save_tsv(cc, out / "diversity_crosscheck.tsv")
    print(f"cross-check (CEN pi, window MSA vs aligned.tsv): {(cc['rel_diff'].abs() <= 0.10).sum()}/{len(cc)} within 10 %, "
          f"{(cc['explained_by'] != '').sum()} more explained by one strain's alignment")
    print(cc[cc["rel_diff"].abs() > 0.10].to_string(index=False))


if __name__ == "__main__":
    main()
