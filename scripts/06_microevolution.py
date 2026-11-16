#!/usr/bin/env python3
"""Microevolution analyses (spec docs/specs/2026-10-20-centromere-microevolution.md).

B. CDEII length in helical steps. Outputs in config `output.dir`:
  cdeii_reference_alleles.tsv  per CEN: median CDEII length and the reference allele
  cdeii_steps.tsv              per locus (10 assemblies + S288C): delta, step class, indels
  cdeii_events.tsv             per indel: type, length, position, tandem unit, microhomology
  cdeii_haplotypes.tsv         non-'within' loci vs the Peter 2018 assembly calls in Helsen SD4
  helsen_cdeii_steps.tsv       same as cdeii_steps for every distinct Helsen SD4 CDEII allele
  helsen_cdeii_events.tsv      same as cdeii_events for those alleles
"""
import argparse
from pathlib import Path

import pandas as pd

from centromere.cde import CDEI_LEN, CDEIII_LEN, read_reference
from centromere.io import load_config, save_tsv
from centromere.microevo import CLASS_ORDER, call_steps, helsen_cdeii, step_class
from centromere.validate import HELSEN_EXTRA_3P, load_helsen_sd4

MIN_EVENT = 3  # bp; shorter indels are homopolymer-run noise for the step question


def _phasing(assembly: str) -> str:
    """ScRAP naming: `_collapsed` = both haplotypes merged; `_1a`/`_3a` = one phased haplotype."""
    a = str(assembly)
    if a.endswith("_collapsed"):
        return "collapsed"
    return "phased_haplotype" if a.rsplit("_", 1)[-1][:1].isdigit() else "homozygous"


def haplotype_evidence(steps: pd.DataFrame, helsen: pd.DataFrame) -> pd.DataFrame:
    """For non-'within' loci: the strain's CDEII lengths in Helsen's Peter 2018 calls."""
    rows = []
    var = steps[(steps["class_median"] != "within") & (steps["strain"] != "S288C")]
    for _, r in var.iterrows():
        calls = helsen[(helsen["assembly"] == r["strain"]) & (helsen["centromere"] == r["centromere"])]
        lengths = sorted(len(s) - CDEI_LEN - CDEIII_LEN - HELSEN_EXTRA_3P for s in calls["sequence"])
        ours, med = r["cdeii_length"], r["median_length"]
        if not lengths:
            verdict = "no_other_assembly"
        elif ours in lengths and any(abs(x - med) <= 2 for x in lengths):
            verdict = "both_alleles_in_other_assembly"
        elif ours in lengths:
            verdict = "confirmed_by_other_assembly"
        elif any(abs(x - med) <= 2 for x in lengths):
            # Collapsed ScRAP assembly: candidate heterozygote (each assembly kept one haplotype)
            verdict = "other_assembly_at_median"
        else:
            verdict = "other_assembly_differs"
        rows.append({"sample": r["sample"], "strain": r["strain"], "assembly": r["assembly"],
                     "centromere": r["centromere"], "cdeii_length": ours, "median_length": med,
                     "class_median": r["class_median"],
                     "phasing": _phasing(r["assembly"]),
                     "other_assembly_cdeii": ",".join(map(str, lengths)), "verdict": verdict})
    return pd.DataFrame(rows)


def summarise_events(events: pd.DataFrame, weight: str | None = None) -> pd.DataFrame:
    """Indels >= MIN_EVENT bp by |length| step class and tandem status (counts, optionally weighted)."""
    e = events[(events["length"].abs() >= MIN_EVENT) & (events["type"] != "complex")
               & ~events["divergent"]].copy()
    e["len_class"] = e["length"].map(step_class)
    e["w"] = e[weight] if weight else 1
    t = e.pivot_table(index="len_class", columns="tandem", values="w", aggfunc="sum", fill_value=0)
    return t.reindex([c for c in CLASS_ORDER if c in t.index])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    out = Path(config["output"]["dir"])
    s288c = read_reference(config["reference"])

    # ---- B on our assemblies (+ S288C as an extra row)
    loci = pd.read_csv(out / "preprocessed.tsv", sep="\t")
    loci = loci[["sample", "assembly", "strain", "centromere", "cdeii_seq"]]
    ref_rows = pd.DataFrame([{"sample": "S288C", "assembly": "S288C", "strain": "S288C",
                              "centromere": c, "cdeii_seq": s[CDEI_LEN:len(s) - CDEIII_LEN]}
                             for c, s in s288c.items()])
    loci = pd.concat([loci, ref_rows], ignore_index=True)
    steps, events, refs = call_steps(loci, s288c, ["sample", "assembly", "strain"])
    save_tsv(refs, out / "cdeii_reference_alleles.tsv")
    save_tsv(steps, out / "cdeii_steps.tsv")
    save_tsv(events, out / "cdeii_events.tsv")

    print("B. ScRAP + S288C, step class vs the population median:")
    print(steps["class_median"].value_counts().reindex(CLASS_ORDER).dropna().astype(int).to_string())
    cols = ["strain", "centromere", "cdeii_length", "delta_median", "class_median",
            "delta_s288c", "events"]
    print(steps[steps["class_median"] != "within"][cols].to_string(index=False))

    sd4 = config.get("validation", {}).get("helsen_sd4")
    if not sd4 or not Path(sd4).exists():
        print("Helsen SD4 not found; skipping haplotype evidence and the SD4 allele scan")
        return
    helsen = load_helsen_sd4(sd4)

    hap = haplotype_evidence(steps, helsen)
    save_tsv(hap, out / "cdeii_haplotypes.tsv")
    print("\nB.4 non-'within' loci vs the Peter 2018 assembly (Helsen SD4):")
    print(hap[["strain", "centromere", "cdeii_length", "phasing", "other_assembly_cdeii",
               "verdict"]].to_string(index=False))

    # ---- B.3 on every distinct Helsen SD4 allele
    alleles = helsen_cdeii(helsen, HELSEN_EXTRA_3P)
    h_steps, h_events, _ = call_steps(alleles, s288c, ["allele", "example"], weight="n_assemblies")
    save_tsv(h_steps, out / "helsen_cdeii_steps.tsv")
    save_tsv(h_events, out / "helsen_cdeii_events.tsv")
    n_asm = helsen["assembly"].nunique()
    print(f"\nHelsen SD4: {len(alleles)} distinct CDEII alleles from {n_asm} assemblies")
    print("Alleles / assemblies per step class (vs the SD4 median):")
    print(h_steps[~h_steps["divergent"]].groupby("class_median")
          .agg(alleles=("allele", "size"), assemblies=("n_assemblies", "sum"))
          .reindex(CLASS_ORDER).dropna().astype(int).to_string())
    div = h_steps["divergent"]
    print(f"Divergent alleles (identity < 0.85, excluded below): {int(div.sum())} "
          f"({int(h_steps.loc[div, 'n_assemblies'].sum())} assemblies)")
    print(f"Indels >= {MIN_EVENT} bp by length class x tandem status, distinct alleles:")
    print(summarise_events(h_events).to_string())
    print("... weighted by assemblies:")
    print(summarise_events(h_events, "n_assemblies").to_string())


if __name__ == "__main__":
    main()
