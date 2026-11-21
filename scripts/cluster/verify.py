#!/usr/bin/env python
"""Sanity checks after a cluster run: is every step done, did any step log an error, did every
genome make it through, how do the locators agree. Prints FAIL / WARN / ok lines; exit code 1 if
any FAIL.

    python scripts/cluster/verify.py              # scrap142 loegler1086
    python scripts/cluster/verify.py loegler1086
"""
import re
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[2]
STEPS = ("blast_cen pipeline flanks_blast l2_flanks minimap2 l3_diversity hmm reconcile cen_like "
         "synthetic cdeii_steps kinetochore kinetochore_domains kinetochore_assoc transposons "
         "ltr_cde diversity trees grammar report").split()
# key outputs that must exist and be non-empty (relative to output.dir)
OUTPUTS = ["centromeres.tsv", "loci_reconciled.tsv", "aligned.tsv", "l2_loci.tsv", "l3_loci.tsv",
           "hmm/hmm_loci.tsv", "cdeii_steps.tsv", "kinetochore_variants.tsv", "te_near_cen.tsv",
           "ltr_cde_scores.tsv", "diversity_summary.tsv", "cen_like_elements.tsv",
           "report/report.html"]
ERR = re.compile(r"Traceback|Error|Killed|MemoryError|CANCELLED|DUE TO TIME LIMIT|oom", re.I)
fails = 0


def say(level, msg):
    global fails
    fails += level == "FAIL"
    print(f"  {level:4} {msg}")


def check(name):
    cfg = yaml.safe_load(open(REPO / f"config.{name}.yaml"))
    out = REPO / cfg["output"]["dir"]
    print(f"== {name}  ({out.relative_to(REPO)})")

    missing = [s for s in STEPS if not (out / ".done" / s).exists()]
    say("FAIL" if missing else "ok", f"steps done {len(STEPS) - len(missing)}/{len(STEPS)}"
        + (f"; missing: {' '.join(missing)}" if missing else ""))

    for log in sorted((out / "logs").glob("*.log")):
        hits = [ln.strip() for ln in log.read_text(errors="replace").splitlines() if ERR.search(ln)]
        hits = [h for h in hits if not re.search(r"\b0 errors?\b|error_rate|stderr", h, re.I)]
        if hits:
            say("WARN", f"{log.name}: {len(hits)} suspicious lines, first: {hits[0][:120]}")

    for f in OUTPUTS:
        p = out / f
        if not p.exists() or p.stat().st_size == 0:
            say("FAIL", f"missing/empty {f}")

    genomes = {p.name.removesuffix(".fasta.gz")
               for p in (REPO / cfg["cde"]["genome_dir"]).glob("*.fasta.gz")}
    smap = REPO / cfg["input"].get("sample_map", "")
    if smap.is_file():
        n = len(pd.read_csv(smap, sep="\t"))
        say("ok" if n == len(genomes) else "WARN", f"genomes on disk {len(genomes)}, sample table {n}")
    blast = {p.name.removesuffix("_db_results.tsv") for p in (out / "blast_cen").glob("*_db_results.tsv")}
    lost = genomes - blast
    say("FAIL" if lost else "ok", f"BLAST results for {len(blast)}/{len(genomes)} genomes"
        + (f"; e.g. no result: {sorted(lost)[:5]}" if lost else ""))

    rec_p = out / "loci_reconciled.tsv"
    if not rec_p.exists():
        return
    rec = pd.read_csv(rec_p, sep="\t")
    ns = rec["sample"].nunique()
    say("ok" if ns == len(genomes) else "FAIL", f"reconciled: {ns} samples x "
        f"{rec['centromere'].nunique()} CEN = {len(rec)} loci (expect {len(genomes) * 16})")
    cats = rec["category"].value_counts()
    print("       categories: " + ", ".join(f"{k} {v} ({v / len(rec):.1%})" for k, v in cats.items()))
    full = rec["ran"].fillna("").eq("L1,L2,L3,L4").mean()
    say("ok" if full > 0.99 else "WARN", f"all four locators ran on {full:.1%} of loci")
    bad = rec[~rec["category"].isin(["agree"])]
    if len(bad):
        print("       non-agree by CEN: " + ", ".join(
            f"{c} {n}" for c, n in bad["centromere"].value_counts().items()))
        worst = bad["sample"].value_counts().head(5)
        print("       samples with most non-agree loci: " + ", ".join(f"{s} {n}" for s, n in worst.items()))
    if cats.get("missing", 0) / len(rec) > 0.01:
        say("WARN", "more than 1% of loci missing")

    aln = out / "aligned.tsv"
    if aln.exists():
        a = pd.read_csv(aln, sep="\t", usecols=["sample", "centromere"])
        say("ok", f"aligned.tsv: {len(a)} loci, {a['sample'].nunique()} samples")
    rep = out / "report" / "report.html"
    if rep.exists():
        say("ok", f"report.html {rep.stat().st_size / 1e6:.1f} MB")


for name in sys.argv[1:] or ["scrap142", "loegler1086"]:
    check(name)
print("FAILED" if fails else "all checks passed")
sys.exit(1 if fails else 0)
