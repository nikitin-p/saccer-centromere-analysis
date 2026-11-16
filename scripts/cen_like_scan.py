#!/usr/bin/env python3
"""Locator spec phase 4: CEN-like scan with the final L4 profiles.

Every assembly's nhmmer hits from the last L4 round (`<hmm_dir>/round_<k>/hits/*.tbl`, E <= 1)
and a fresh nhmmer run on the S288C genome are filtered by `hmm.scan_cen_like`: E <=
`locate.cen_like_evalue`, best profile at the spot by >= `locate.hmm_margin` bits, CDEI and
CDEIII pass `cde_ok`, and no overlap with any known CEN locus of that genome (the reconciled
main table plus every L1-L4 locus; for S288C the 16 SGD CENs). These are candidates
(neocentromeres, dicentric arms, CEN fragments), not calls.

Outputs in `<output.dir>/`:
  cen_like_elements.tsv   candidates passing every check (sample = accession or S288C)
  cen_like_rejected.tsv   spots passing E and margin but failing the CDE gate (context)
  cen_like_nearest.tsv    per genome: the strongest hit outside every known locus, no filter

  PYTHONPATH=src .venv/bin/python scripts/cen_like_scan.py --config config.scrap.yaml
"""
import argparse
from pathlib import Path

import pandas as pd

from centromere.cde import CDEModel, read_reference
from centromere.hmm import add_margin, model_lengths, overlaps, parse_tblout, run_nhmmer, scan_cen_like
from centromere.io import load_config, save_tsv, tool
from centromere.locate import assemblies, locator_loci, read_fasta


def last_round(hmm_dir: Path) -> Path:
    rounds = sorted(hmm_dir.glob("round_*"), key=lambda p: int(p.name.split("_")[1]))
    if not rounds:
        raise SystemExit(f"no L4 rounds in {hmm_dir}; run scripts/hmm_iterate.py first")
    return rounds[-1]


def nearest(hits, known):
    """Strongest hit (lowest E) not overlapping a known locus, with its margin."""
    h = add_margin(hits)
    h = h[[not overlaps(c, lo, hi, known) for c, lo, hi in zip(h["contig"], h["lo"], h["hi"])]]
    return h.sort_values(["evalue", "bits"], ascending=[True, False]).head(1).drop(columns=["lo", "hi"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    loc, out = cfg["locate"], Path(cfg["output"]["dir"])
    rdir = last_round(Path(loc["hmm_dir"]))
    allh = rdir / "all.hmm"
    lengths = model_lengths(allh)
    model = CDEModel(read_reference(cfg["reference"]))
    p = {"evalue": loc.get("cen_like_evalue", 1e-5), "margin": loc.get("hmm_margin", 10),
         "flank": cfg["cde"]["flank"], "radius": cfg["cde"]["search_radius"]}
    print(f"profiles: {allh} ({len(lengths)} models); E <= {p['evalue']}, margin >= {p['margin']}")

    main_tab = pd.read_csv(out / "preprocessed.tsv", sep="\t")
    loci, _ = locator_loci(cfg, main_tab)
    known = pd.concat(list(loci.values()), ignore_index=True)[["sample", "contig", "start", "end"]]

    rows, near = [], []
    for _, a in assemblies(cfg).iterrows():
        hits = parse_tblout(rdir / "hits" / f"{a['accession']}.tbl")
        k = known[known["sample"] == a["accession"]]
        s = scan_cen_like(hits, read_fasta(a["fasta"]), lengths, model, k, **p)
        rows.append(s.assign(sample=a["accession"], strain=a["strain"]))
        near.append(nearest(hits, k).assign(sample=a["accession"], strain=a["strain"]))

    # Spec check 9: S288C itself, known = the 16 SGD CENs (BED is 0-based)
    genome = Path(loc["paf_dir"]) / "sacCer3.nuclear.fa"
    if genome.exists():
        nhmmer = tool(loc.get("hmmer_bin"), "nhmmer")   # null -> PATH
        tbl = out / "cen_like_s288c.tbl"
        hits = run_nhmmer(allh, genome, tbl, nhmmer, loc.get("threads", 4), loc.get("hmm_report_evalue", 1.0))
        bed = pd.read_csv(cfg["transposons"]["s288c_cen"], sep="\t", header=None,
                          names=["contig", "start", "end", "centromere", "score", "strand"])
        s288c_known = bed.assign(start=bed["start"] + 1)[["contig", "start", "end"]]
        own = hits[hits["evalue"] <= loc.get("hmm_evalue", 1e-10)]
        print(f"S288C: {own['profile'].nunique()}/16 profiles hit at E <= {loc.get('hmm_evalue', 1e-10)}")
        s = scan_cen_like(hits, read_fasta(genome), lengths, model, s288c_known, **p)
        rows.append(s.assign(sample="S288C", strain="S288C"))
        near.append(nearest(hits, s288c_known).assign(sample="S288C", strain="S288C"))
    else:
        print(f"S288C genome {genome} missing (written by run_minimap2.sh); check 9 skipped")

    df = pd.concat(rows, ignore_index=True)
    front = ["sample", "strain", "profile", "contig", "start", "end", "strand", "length"]
    df = df.reindex(columns=front + [c for c in df.columns if c not in front])
    ok = df["cde_ok"].astype(bool) if len(df) else pd.Series(dtype=bool)
    save_tsv(df[ok], out / "cen_like_elements.tsv")
    save_tsv(df[~ok], out / "cen_like_rejected.tsv")
    nr = pd.concat(near, ignore_index=True)
    save_tsv(nr[["sample", "strain"] + [c for c in nr.columns if c not in ("sample", "strain")]],
             out / "cen_like_nearest.tsv")
    print("strongest hit outside known loci, per genome: E "
          f"{nr['evalue'].min():.2g}-{nr['evalue'].max():.2g}, {nr['bits'].min()}-{nr['bits'].max()} bits, "
          f"profiles {nr['profile'].value_counts().to_dict()}")
    print(f"CEN-like candidates: {int(ok.sum())} pass all checks, {int((~ok).sum())} fail the CDE gate")
    if len(df):
        cols = ["strain", "profile", "contig", "start", "end", "length", "bits", "evalue", "margin",
                "hmm_cov", "cdei_score", "cdeiii_score", "cde_ok"]
        print(df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
