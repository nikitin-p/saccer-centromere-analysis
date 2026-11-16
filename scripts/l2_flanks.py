#!/usr/bin/env python3
"""L2 locator (locator spec phase 1): CEN loci from the 2-kb S288C flanks, compared with L1 and L3.

Needs `blast.sh MODE=flanks` output (`locate.flank_hits`), `<output.dir>/filtered.tsv` (L1) and,
if present, `<output.dir>/l3_loci.tsv` (L3, from scripts/12_flanks_diversity.py).
Writes `<output.dir>/l2_loci.tsv`: one row per assembly x CEN with the L2 locus (1-based,
inclusive), the reason when there is none, and the differences to L1 / L3 (d_* = L2 - other,
start and end; agree_* = same contig and overlapping).

  PYTHONPATH=src .venv/bin/python scripts/l2_flanks.py --config config.scrap.yaml
"""
import argparse
from pathlib import Path

import pandas as pd

from centromere.io import load_config, save_tsv
from centromere.locate import locate_by_flanks, read_flank_hits


def compare(l2: pd.DataFrame, other: pd.DataFrame, tag: str) -> pd.DataFrame:
    """Add {tag}_contig/start/end, d_start_{tag}, d_end_{tag}, agree_{tag}."""
    o = other.rename(columns={"contig": f"{tag}_contig", "start": f"{tag}_start", "end": f"{tag}_end"})
    m = l2.merge(o[["sample", "centromere", f"{tag}_contig", f"{tag}_start", f"{tag}_end"]],
                 on=["sample", "centromere"], how="left")
    m[f"d_start_{tag}"] = m["l2_start"] - m[f"{tag}_start"]
    m[f"d_end_{tag}"] = m["l2_end"] - m[f"{tag}_end"]
    ov = m[["l2_end", f"{tag}_end"]].min(axis=1) - m[["l2_start", f"{tag}_start"]].max(axis=1) + 1
    m[f"agree_{tag}"] = ((m["contig"] == m[f"{tag}_contig"]) & (ov > 0)).astype(object)
    m.loc[m[f"{tag}_contig"].isna() | (m["status"] != "ok"), f"agree_{tag}"] = pd.NA
    return m


def report(m: pd.DataFrame, tag: str) -> None:
    has = m[f"agree_{tag}"].notna()
    n = int(m.loc[has, f"agree_{tag}"].astype(bool).sum())
    exact = int(((m["d_start_" + tag] == 0) & (m["d_end_" + tag] == 0)).sum())
    print(f"L2 vs {tag.upper()}: {n}/{int(has.sum())} agree (same contig, overlapping); "
          f"{exact} with identical boundaries")
    off = m[has & ((m["d_start_" + tag] != 0) | (m["d_end_" + tag] != 0))]
    if len(off):
        print(off[["strain", "centromere", "l2_length", f"d_start_{tag}", f"d_end_{tag}",
                   "ext_left", "ext_right", f"agree_{tag}"]].to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    cfg = load_config(ap.parse_args().config)
    loc, out = cfg["locate"], Path(cfg["output"]["dir"])

    hits = read_flank_hits(loc["flank_hits"])
    l2 = locate_by_flanks(hits, loc.get("min_flank_pident", 90), loc.get("min_flank_cov", 0.5),
                          loc.get("max_locus_bp", 20000))
    smap = pd.read_csv(cfg["input"]["sample_map"], sep="\t").set_index("accession")["strain"]
    l2.insert(1, "strain", l2["sample"].map(smap))
    print(f"L2 status: {l2['status'].value_counts().to_dict()}; "
          f">1 passing flank hit: {int(((l2['n_left'] > 1) | (l2['n_right'] > 1)).sum())}")

    f = pd.read_csv(out / "filtered.tsv", sep="\t")
    l1 = f.rename(columns={"sseqid": "contig", "cen_start": "start", "cen_end": "end"})
    m = compare(l2, l1, "l1")
    report(m, "l1")
    l3_path = out / "l3_loci.tsv"
    if l3_path.exists():
        l3 = pd.read_csv(l3_path, sep="\t")
        l3 = l3[l3["status"].isin(["lifted", "bridged"])].assign(sample=lambda d: d["strain"].map(
            {v: k for k, v in smap.items()}))
        l3 = l3.rename(columns={"l3_start": "start", "l3_end": "end"})
        m = compare(m, l3, "l3")
        report(m, "l3")
    save_tsv(m, out / "l2_loci.tsv")
    print(f"wrote {out / 'l2_loci.tsv'} ({len(m)} rows)")


if __name__ == "__main__":
    main()
