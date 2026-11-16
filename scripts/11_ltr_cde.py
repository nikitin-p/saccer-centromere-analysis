#!/usr/bin/env python3
"""Analysis E: CDE-like motifs (CDEI, CDEIII, full units) in Ty LTRs vs controls.

Needs analysis D outputs (te_elements.tsv, trna_loci.tsv, transposons/<g>.te.tsv).
Outputs in `<output.dir>/` (spec docs/specs/2026-11-13-ltr-cde-motifs.md):
  ltr_seqs.fasta           every LTR, oriented like its element (header: ltr_id seq_id)
  ltr_cde_scores.tsv       one row per LTR: best CDEI/CDEIII scores, hit and unit counts
  ltr_cde_control.tsv      observed vs shuffled / genomic / S288C-intergenic rates
  ltr_cde_poscontrol.tsv   scanner on the 176 CEN loci (160 ScRAP + 16 S288C)
  ltr_cde_candidates.tsv   LTRs containing a full CDEI-CDEII-CDEIII unit
  ltr_cde_score_dist.tsv   best CDEI/CDEIII score per sequence for LTR/shuffled/genomic (figure)
  ltr_cde_genome_units.tsv every full unit genome-wide, with the CEN it falls in (if any)
  ltr_cde_thresholds.tsv   the thresholds used
"""
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from statsmodels.stats.multitest import multipletests

from centromere.cde import CDEI_LEN, CDEModel, read_reference
from centromere.genes import make_db
from centromere.io import load_config, save_tsv
from centromere.ltr_cde import (CDEScanner, compare, dinuc_shuffle, extract, ltr_intervals, merge,
                                revcomp, sample_windows, window_counts)
from centromere.transposons import HIT_COLS, interval_distance, parse_hits

COUNTS = ["n_cdei", "n_cdeiii", "n_units"]
# SGD feature types that are not "intergenic" (gene-free control on S288C)
GENIC = {"gene", "tRNA_gene", "snoRNA_gene", "ncRNA_gene", "rRNA_gene", "snRNA_gene",
         "telomerase_RNA_gene", "transposable_element_gene", "blocked_reading_frame",
         "pseudogene", "LTR_retrotransposon", "long_terminal_repeat", "centromere",
         "telomere", "X_element", "Y_prime_element"}


def genomes(config: dict) -> list[tuple[str, Path]]:
    samples = pd.read_csv(config["input"]["sample_map"], sep="\t")
    gdir = Path(config["cde"]["genome_dir"])
    out = [(r["strain"], gdir / f"{r['accession']}.fasta.gz") for _, r in samples.iterrows()]
    return out + [("S288C", Path(config["kinetochore"]["s288c_genome"]).expanduser())]


def cen_table(config: dict) -> pd.DataFrame:
    f = pd.read_csv(Path(config["output"]["dir"]) / "filtered.tsv", sep="\t")
    cens = f[["strain", "centromere", "sseqid", "cen_start", "cen_end", "strand", "cdeii_seq",
              "cde_ok", "cdeii_length"]].rename(columns={"strain": "genome", "sseqid": "contig"})
    bed = pd.read_csv(config["transposons"]["s288c_cen"], sep="\t", header=None,
                      names=["contig", "s0", "cen_end", "centromere", "score", "strand"])
    ref = read_reference(config["reference"])
    s288c = bed.assign(genome="S288C", cen_start=bed["s0"] + 1,
                       cdeii_seq=bed["centromere"].map(lambda c: ref[c][CDEI_LEN:-26]), cde_ok=True)
    s288c["cdeii_length"] = s288c["cdeii_seq"].str.len()
    return pd.concat([cens, s288c[cens.columns]], ignore_index=True)


def ltr_lengths(library: str) -> dict[str, int]:
    lens: dict[str, list[int]] = {}
    for rec in SeqIO.parse(library, "fasta"):
        kind, fam, _ = rec.id.split("|")
        if kind == "LTR":
            lens.setdefault(fam, []).append(len(rec.seq))
    return {f: int(np.median(v)) for f, v in lens.items()}


def genic_intervals(gff: str) -> dict[str, np.ndarray]:
    iv: dict[str, list] = {}
    with gzip.open(gff, "rt") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) == 9 and p[2] in GENIC:
                iv.setdefault(p[0], []).append((int(p[3]) - 1, int(p[4])))
    return {c: merge(np.array(v)) for c, v in iv.items()}


def positive_control(scanner: CDEScanner, cens: pd.DataFrame, contigs: dict, flank: int = 150) -> pd.DataFrame:
    rows = []
    for _, c in cens.iterrows():
        seq = contigs[c["genome"]][c["contig"]]
        lo = max(0, c["cen_start"] - 1 - flank)
        win = seq[lo:c["cen_end"] + flank].upper()
        off = c["cen_start"] - 1 - lo
        if c["strand"] == "-":
            win, off = revcomp(win), (c["cen_end"] + flank if c["cen_end"] + flank <= len(seq)
                                      else len(seq)) - c["cen_end"]
        u = scanner.hits(win).units
        u = u[(u["strand"] == 1) & ((u["lo"] - off).abs() <= 2)]
        same = [win[a + CDEI_LEN:b - 26] == c["cdeii_seq"] for a, b in zip(u["lo"], u["hi"])]
        in_def = bool(c["cde_ok"]) and scanner.cdeii_min <= c["cdeii_length"] <= scanner.cdeii_max
        rows.append({"genome": c["genome"], "centromere": c["centromere"], "found": len(u) > 0,
                     "same_cdeii_as_cde_py": any(same), "meets_definition": in_def,
                     "cdeii_length": c["cdeii_length"], "cde_ok": c["cde_ok"]})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    lcfg, tcfg = config["ltr_cde"], config["transposons"]
    out = Path(config["output"]["dir"])
    blast_bin = config["kinetochore"].get("blast_bin")   # null -> PATH
    rng = np.random.default_rng(lcfg["seed"])
    model = CDEModel(read_reference(config["reference"]))
    scanner = CDEScanner(model, lcfg["cdeii_min_len"], lcfg["cdeii_max_len"], lcfg["cdeii_min_at"])
    print(f"Thresholds: CDEI >= {model.cdei_min:.2f}, CDEIII >= {model.cdeiii_min:.2f}, CDEII "
          f"{scanner.cdeii_min}-{scanner.cdeii_max} bp, AT >= {scanner.cdeii_min_at}")

    elements = pd.read_csv(out / "te_elements.tsv", sep="\t")
    trnas = pd.read_csv(out / "trna_loci.tsv", sep="\t")
    cens = cen_table(config)
    flen = ltr_lengths(tcfg["library"])
    contigs = {}
    ltrs = []
    for g, path in genomes(config):
        fa = make_db(path, out / "kinetochore_blastdb", blast_bin)
        contigs[g] = {r.id: str(r.seq) for r in SeqIO.parse(fa, "fasta")}
        raw = pd.read_csv(out / "transposons" / f"{g}.te.tsv", sep="\t", header=None, names=HIT_COLS)
        ltrs.append(ltr_intervals(elements[elements["genome"] == g], parse_hits(raw), flen,
                                  tcfg["min_pident"], tcfg["min_qcov"]))
    ltrs = pd.concat(ltrs, ignore_index=True)
    ltrs["ltr_id"] = ltrs["element_id"] + ":" + ltrs["ltr"]
    ltrs["seq"] = [extract(contigs[r.genome][r.contig], r.start, r.end, r.strand) for r in ltrs.itertuples()]
    uniq = {s: f"s{i}" for i, s in enumerate(dict.fromkeys(ltrs["seq"]))}
    ltrs["seq_id"] = ltrs["seq"].map(uniq)
    with open(out / "ltr_seqs.fasta", "w") as fh:
        for r in ltrs.itertuples():
            fh.write(f">{r.ltr_id} {r.seq_id}\n{r.seq}\n")
    print(f"LTRs: {len(ltrs)} ({ltrs['ltr'].value_counts().to_dict()}; full-element LTR source "
          f"{ltrs.loc[ltrs['ltr'] != 'solo', 'source'].value_counts().to_dict()}), "
          f"{len(uniq)} unique sequences")

    # ---- positive control
    pc = positive_control(scanner, cens, contigs)
    save_tsv(pc, out / "ltr_cde_poscontrol.tsv")
    d = pc[pc["meets_definition"]]
    print(f"\nPositive control: unit found at {pc['found'].sum()}/{len(pc)} CEN loci; "
          f"{d['found'].sum()}/{len(d)} of loci meeting the definitions; "
          f"same CDEII as cde.py: {pc['same_cdeii_as_cde_py'].sum()}/{len(pc)}")
    if (~pc["found"]).any():
        print(pc[~pc["found"]].to_string(index=False))
    if not d["found"].all():
        raise SystemExit("Scanner misses CEN loci that meet the definitions - stop (spec step 4).")

    # ---- scan unique LTR sequences
    seqs = pd.DataFrame({"seq": list(uniq), "seq_id": list(uniq.values())})
    summ = pd.DataFrame([scanner.summary(s) for s in seqs["seq"]])
    seqs = pd.concat([seqs, summ], axis=1)
    fam = ltrs.groupby("seq_id")["family"].agg(lambda x: x.mode().iloc[0])
    seqs["family"] = seqs["seq_id"].map(fam)
    scores = ltrs.drop(columns="seq").merge(seqs.drop(columns=["seq", "family"]), on="seq_id")
    scores["n_strains_with_seq"] = scores["seq_id"].map(ltrs.groupby("seq_id")["genome"].nunique())
    save_tsv(scores, out / "ltr_cde_scores.tsv")

    # ---- controls
    rows = []
    # (a) dinucleotide shuffles of each unique sequence
    shuf = np.zeros((len(seqs), lcfg["n_shuffles"], 3), int)
    shuf_best = []
    for i, s in enumerate(seqs["seq"]):
        for k in range(lcfg["n_shuffles"]):
            sh = dinuc_shuffle(s, rng)
            shuf[i, k] = scanner.counts(sh)
            if k == 0:
                b = scanner.summary(sh)
                shuf_best.append({"set": "shuffled", "best_cdei": b["best_cdei"], "best_cdeiii": b["best_cdeiii"]})
    for fam_name, sel in [("all", np.ones(len(seqs), bool))] + \
            [(f, (seqs["family"] == f).to_numpy()) for f in sorted(seqs["family"].unique())]:
        obs = seqs.loc[sel, COUNTS].sum().to_dict()
        bp = seqs.loc[sel, "length"].sum()
        null = pd.DataFrame(shuf[sel].sum(axis=0), columns=COUNTS)
        for r in compare(obs, bp, null, np.full(len(null), bp)):
            rows.append({"control": "shuffled", "set": f"unique LTRs ({fam_name})", "n_seq": int(sel.sum()),
                         "bp": int(bp), **r})

    # (b) genomic and (c) S288C intergenic windows, matched to each genome's LTR lengths
    genic = genic_intervals(lcfg["sgd_gff"])
    gen_reps, gen_obs, gen_best, gen_units = [], [], [], []
    for g in contigs:
        L = ltrs[ltrs["genome"] == g]
        if L.empty:
            continue
        hits = {c: scanner.hits(s) for c, s in contigs[g].items()}
        for c, h in hits.items():
            for u in h.units.itertuples():
                cc = cens[(cens["genome"] == g) & (cens["contig"] == c)]
                hit_cen = cc[(cc["cen_start"] - 1 < u.hi) & (cc["cen_end"] > u.lo)]
                gen_units.append({"genome": g, "contig": c, "start": u.lo + 1, "end": u.hi,
                                  "strand": "+" if u.strand == 1 else "-", "cdei_score": round(u.cdei_score, 2),
                                  "cdeiii_score": round(u.cdeiii_score, 2), "cdeii_len": u.cdeii_len,
                                  "centromere": hit_cen["centromere"].iloc[0] if len(hit_cen) else None})
        blocked = {}
        for c in contigs[g]:
            iv = [elements.loc[(elements["genome"] == g) & (elements["contig"] == c), ["start", "end"]].to_numpy()
                  - np.array([1, 0]),
                  trnas.loc[(trnas["genome"] == g) & (trnas["contig"] == c), ["start", "end"]].to_numpy()
                  - np.array([1, 0])]
            cc = cens[(cens["genome"] == g) & (cens["contig"] == c)]
            iv.append(np.c_[cc["cen_start"] - 1 - lcfg["cen_exclusion"], cc["cen_end"] + lcfg["cen_exclusion"]])
            blocked[c] = merge(np.vstack([x.reshape(-1, 2) for x in iv]).astype(int))
        lengths = L["end"].to_numpy() - L["start"].to_numpy() + 1
        clen = {c: len(s) for c, s in contigs[g].items()}
        obs = scores.loc[scores["genome"] == g, COUNTS].sum().to_dict()
        controls = [("genomic", blocked)]
        if g == "S288C":
            controls.append(("intergenic", {c: merge(np.vstack([blocked[c], genic.get(c, np.empty((0, 2)))]))
                                            for c in blocked}))
        for name, blk in controls:
            win = sample_windows(clen, blk, lengths, lcfg["n_resamples"], rng, allow_missing=(name == "intergenic"))
            if win.attrs["n_missing"]:
                print(f"  {g} {name}: {win.attrs['n_missing']} of {len(lengths) * lcfg['n_resamples']} "
                      f"windows did not fit (dropped; rates use placed bp)")
            w = window_counts(hits, win)
            w["bp"] = w["hi"] - w["lo"]
            null = w.groupby("rep")[COUNTS + ["bp"]].sum()
            for r in compare(obs, lengths.sum(), null[COUNTS], null["bp"].to_numpy()):
                rows.append({"control": name, "set": f"{g} LTRs", "n_seq": len(L), "bp": int(lengths.sum()), **r})
            null = null[COUNTS]
            if name == "genomic":
                gen_reps.append(null.to_numpy())
                gen_obs.append(obs)
                first = w[w["rep"] == 0]
                gen_best += [{"set": "genomic", **{k: v for k, v in scanner.summary(
                    extract(contigs[g][r.contig], r.lo + 1, r.hi, "+")).items() if k in ("best_cdei", "best_cdeiii")}}
                    for r in first.itertuples()]
        print(f"{g}: {len(L)} LTRs, genome-wide units {sum(len(h.units) for h in hits.values())}")
    gu = pd.DataFrame(gen_units)
    save_tsv(gu, out / "ltr_cde_genome_units.tsv")
    print(f"Genome-wide units: {len(gu)}, outside a CEN locus: {gu['centromere'].isna().sum()}")
    save_tsv(pd.DataFrame([{"cdei_min": round(model.cdei_min, 2), "cdeiii_min": round(model.cdeiii_min, 2),
                            "cdeii_min_len": scanner.cdeii_min, "cdeii_max_len": scanner.cdeii_max,
                            "cdeii_min_at": scanner.cdeii_min_at, "n_shuffles": lcfg["n_shuffles"],
                            "n_resamples": lcfg["n_resamples"]}]), out / "ltr_cde_thresholds.tsv")
    tot_bp = ltrs.eval("end - start + 1").sum()
    obs = pd.DataFrame(gen_obs).sum().to_dict()
    null = pd.DataFrame(np.sum(gen_reps, axis=0), columns=COUNTS)
    for r in compare(obs, tot_bp, null, np.full(len(null), tot_bp)):
        rows.append({"control": "genomic", "set": "all genomes LTRs", "n_seq": len(ltrs), "bp": int(tot_bp), **r})
    ctl = pd.DataFrame(rows)
    ctl["q"] = np.nan
    for name, idx in ctl.groupby("control").groups.items():
        ctl.loc[idx, "q"] = multipletests(ctl.loc[idx, "p"], method="fdr_bh")[1]
    save_tsv(ctl, out / "ltr_cde_control.tsv")
    print("\nControls (pooled rows):")
    print(ctl[ctl["set"].isin(["unique LTRs (all)", "all genomes LTRs", "S288C LTRs"])].to_string(index=False))

    dist = pd.concat([seqs[["best_cdei", "best_cdeiii"]].assign(set="LTR"), pd.DataFrame(shuf_best),
                      pd.DataFrame(gen_best)], ignore_index=True)
    save_tsv(dist, out / "ltr_cde_score_dist.tsv")

    # ---- candidates: LTRs with a full unit
    cand = scores[scores["n_units"] > 0].copy()
    if not cand.empty:
        plus = cand["strand"] == "+"
        cand["unit_contig_start"] = np.where(plus, cand["start"] + cand["unit_lo"], cand["end"] - cand["unit_hi"] + 1)
        cand["unit_contig_end"] = np.where(plus, cand["start"] + cand["unit_hi"] - 1, cand["end"] - cand["unit_lo"])
        near = []
        for r in cand.itertuples():
            cc = cens[(cens["genome"] == r.genome) & (cens["contig"] == r.contig)]
            if cc.empty:
                near.append((None, np.nan))
                continue
            dd = interval_distance(r.start, r.end, cc["cen_start"].to_numpy(), cc["cen_end"].to_numpy())
            near.append((cc["centromere"].iloc[int(np.argmin(dd))], int(dd.min())))
        cand["nearest_cen"], cand["nearest_cen_bp"] = zip(*near)
        cand["strains_with_seq"] = cand["seq_id"].map(ltrs.groupby("seq_id")["genome"].agg(lambda x: ",".join(sorted(set(x)))))
    save_tsv(cand, out / "ltr_cde_candidates.tsv")
    print(f"\nCandidates: {len(cand)} LTRs, {cand['seq_id'].nunique() if len(cand) else 0} unique sequences")
    if len(cand):
        print(cand[["ltr_id", "family", "class", "seq_id", "unit_cdei_score", "unit_cdeiii_score",
                    "unit_cdeii_len", "unit_cdeii_at", "nearest_cen", "nearest_cen_bp"]].to_string(index=False))


if __name__ == "__main__":
    main()
