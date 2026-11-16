#!/usr/bin/env python3
"""L4 locator (locator spec phase 3): iterative, expanding per-CEN nucleotide HMMs.

Round 0 profiles are seeded from `<output.dir>/aligned.tsv` + S288C (`--seed aligned`) or from
the 16 S288C CENs alone (`--seed reference`, shows the profiles growing). Each round searches
every training assembly with all 16 profiles, accepts new sequences (rules in
`centromere.hmm`), realigns and rebuilds; it stops when a round accepts nothing or after
`locate.max_rounds`. A CEN whose `runner_up` rejections rise between rounds is frozen (drift).

`--holdout S1,S2`: those strains are left out of the seed and of the iteration, then searched
with the final profiles (generalisation check, spec check 7).

Outputs in `--out` (default `locate.hmm_dir`):
  round_<k>/CEN*.afa, CEN*.hmm, all.hmm, hits/<accession>.tbl
  hmm_rounds.tsv   round x CEN: profile size, accepted, rejected by reason, S288C self-score
  hmm_profile.tsv  every profile member: round added, strain, locus, sequence
  hmm_loci.tsv     final round's loci, one row per assembly x CEN, with L1/L2/L3 agreement
  hmm_holdout.tsv  (with --holdout) the held-out assemblies' loci from the final profiles

  PYTHONPATH=src .venv/bin/python scripts/hmm_iterate.py --config config.scrap.yaml
"""
import argparse
import gzip
import shutil
import tempfile
from pathlib import Path

import pandas as pd

from centromere.cde import CDEModel, read_reference
from centromere.hmm import (build_profile, drifting, locate_in_genome, model_lengths,
                            new_members, round_summary, run_nhmmer, write_fasta)
from centromere.io import load_config, save_tsv, tool
from centromere.locate import assemblies, compare_loci, l1_l2_l3_loci, read_fasta


def build_round(k, profile, cens, out, tools, mafft_args):
    rdir = out / f"round_{k}"
    (rdir / "hits").mkdir(parents=True, exist_ok=True)
    paths = [build_profile(c, sorted(profile[c]), rdir, tools["hmmbuild"], mafft_args) for c in cens]
    allh = rdir / "all.hmm"
    with open(allh, "w") as fh:
        for p in paths:
            fh.write(p.read_text())
    return rdir, allh


def search(asm, allh, rdir, genomes, model, ref_len, lengths, anchors, p, tools):
    rows = []
    for _, a in asm.iterrows():
        hits = run_nhmmer(allh, genomes[a["accession"]], rdir / "hits" / f"{a['accession']}.tbl",
                          tools["nhmmer"], p["cpu"], p["report_evalue"])
        anc = anchors[anchors["sample"] == a["accession"]] if anchors is not None else None
        loci = locate_in_genome(hits, read_fasta(genomes[a["accession"]]), lengths, model, ref_len,
                                anc, p["evalue"], p["margin"], p["len_tol"], p["flank"], p["radius"])
        rows.append(loci.assign(sample=a["accession"], strain=a["strain"]))
    df = pd.concat(rows, ignore_index=True)
    front = ["sample", "strain", "centromere", "status"]
    return df[front + [c for c in df.columns if c not in front]]


def s288c_self_bits(allh, ref_fa, rdir, tools, p):
    h = run_nhmmer(allh, ref_fa, rdir / "s288c.tbl", tools["nhmmer"], p["cpu"], p["report_evalue"])
    own = h[h["profile"] == h["contig"]]
    return own.groupby("profile")["bits"].max().to_dict()


def with_agreement(loci, others):
    ok = loci["status"].isin(["ok", "insertion"])
    m = loci.rename(columns={"hmm_start": "l4_start", "hmm_end": "l4_end"})
    for tag, o in others.items():
        m = compare_loci(m, o, tag, "l4", ok)
    return m.rename(columns={"l4_start": "hmm_start", "l4_end": "hmm_end"})


def summary(m, others, label):
    st = m["status"].value_counts().to_dict()
    print(f"{label}: {len(m)} assembly x CEN, status {st}")
    for tag in others:
        has = m[f"agree_{tag}"].notna()
        n = int(m.loc[has, f"agree_{tag}"].astype(bool).sum())
        exact = int(((m[f"d_start_{tag}"] == 0) & (m[f"d_end_{tag}"] == 0)).sum())
        print(f"  vs {tag.upper()}: {n}/{int(has.sum())} agree; {exact} identical boundaries")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--seed", choices=["aligned", "reference"], default="aligned")
    ap.add_argument("--holdout", default="", help="comma-separated strains left out of training")
    ap.add_argument("--out", help="output directory (default locate.hmm_dir)")
    args = ap.parse_args()
    cfg = load_config(args.config)
    loc, outdir = cfg["locate"], Path(cfg["output"]["dir"])
    out = Path(args.out or loc["hmm_dir"])
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    tools = {t: tool(loc.get("hmmer_bin"), t) for t in ("hmmbuild", "nhmmer")}   # null -> PATH
    p = {"cpu": loc.get("threads", 4), "evalue": loc.get("hmm_evalue", 1e-10),
         "report_evalue": loc.get("hmm_report_evalue", 1.0), "margin": loc.get("hmm_margin", 10),
         "len_tol": loc.get("hmm_len_tol", 40), "flank": cfg["cde"]["flank"],
         "radius": cfg["cde"]["search_radius"]}
    max_rounds = loc.get("max_rounds", 5)
    mafft_args = cfg["align"]["mafft_args"]

    ref = read_reference(cfg["reference"])
    cens = list(ref)
    ref_len = {c: len(s) for c, s in ref.items()}
    model = CDEModel(ref)
    others = l1_l2_l3_loci(cfg)
    anchors = pd.concat([others[t] for t in ("l2", "l3") if t in others]) if \
        {"l2", "l3"} & set(others) else None

    asm = assemblies(cfg)
    hold = [s for s in args.holdout.split(",") if s]
    unknown = set(hold) - set(asm["strain"])
    if unknown:
        raise SystemExit(f"unknown holdout strains: {sorted(unknown)}")
    train, test = asm[~asm["strain"].isin(hold)], asm[asm["strain"].isin(hold)]

    # Round-0 profile: sequence -> where it came from
    members = [{"round": 0, "centromere": c, "strain": "S288C", "sequence": s} for c, s in ref.items()]
    if args.seed == "aligned":
        a = pd.read_csv(outdir / "aligned.tsv", sep="\t")
        a = a[a["strain"].isin(train["strain"])]
        members += [{"round": 0, "centromere": r.centromere, "strain": r.strain, "sequence": r.sequence,
                     "contig": r.sseqid, "start": r.cen_start, "end": r.cen_end, "source": "L1"}
                    for r in a.itertuples()]
    members = pd.DataFrame(members).drop_duplicates(["centromere", "sequence"])
    profile = {c: set(members.loc[members["centromere"] == c, "sequence"]) for c in cens}

    tmp = Path(tempfile.mkdtemp(prefix="hmm_genomes_"))
    try:
        genomes = {}
        for _, r in asm.iterrows():   # nhmmer needs a rewindable (uncompressed) target
            genomes[r["accession"]] = tmp / f"{r['accession']}.fa"
            with gzip.open(r["fasta"], "rb") as src, open(genomes[r["accession"]], "wb") as dst:
                shutil.copyfileobj(src, dst)
        ref_fa = tmp / "s288c_cen.fa"
        write_fasta(ref_fa, list(ref.items()))

        rounds, frozen = [], set()
        for k in range(max_rounds + 1):
            rdir, allh = build_round(k, profile, cens, out, tools, mafft_args)
            lengths = model_lengths(allh)
            loci = search(train, allh, rdir, genomes, model, ref_len, lengths, anchors, p, tools)
            new = new_members(loci, profile, frozen)
            rounds.append(round_summary(k, loci, profile, new, s288c_self_bits(allh, ref_fa, rdir, tools, p), cens))
            frozen |= drifting(pd.concat(rounds), k)
            print(f"round {k}: profile {sum(map(len, profile.values()))} sequences, "
                  f"{int((loci['status'] == 'ok').sum())}/{len(loci)} ok, accepted {len(new)}"
                  + (f", frozen {sorted(frozen)}" if frozen else ""))
            if k == max_rounds or new.empty:
                break
            for r in new.itertuples():
                profile[r.centromere].add(r.sequence)
            members = pd.concat([members, new.assign(round=k + 1, source="L4").rename(
                columns={"hmm_start": "start", "hmm_end": "end"})[
                ["round", "centromere", "strain", "sequence", "contig", "start", "end", "source"]]],
                ignore_index=True)

        save_tsv(pd.concat(rounds, ignore_index=True), out / "hmm_rounds.tsv")
        save_tsv(members, out / "hmm_profile.tsv")
        m = with_agreement(loci, others)
        save_tsv(m, out / "hmm_loci.tsv")
        summary(m, others, f"final round {k}")
        acc = members[members["source"] == "L4"]
        if len(acc) and anchors is not None:   # spec check 8: every accepted locus anchored by L2/L3
            a = with_agreement(acc.rename(columns={"start": "hmm_start", "end": "hmm_end"})
                               .assign(status="ok", sample=lambda d: d["strain"].map(
                                   asm.set_index("strain")["accession"])),
                               {t: others[t] for t in ("l2", "l3") if t in others})
            cols = [c for c in ("agree_l2", "agree_l3") if c in a]
            n = int(a[cols].fillna(False).astype(bool).any(axis=1).sum())
            print(f"accepted by L4: {len(a)}; anchored by L2 or L3: {n}/{len(a)}")

        if len(test):
            h = search(test, allh, rdir, genomes, model, ref_len, lengths, anchors, p, tools)
            h = with_agreement(h, others)
            save_tsv(h, out / "hmm_holdout.tsv")
            summary(h, others, f"holdout {','.join(hold)}")
            rec = h.assign(hit=(h["status"] == "ok") & h["agree_l1"].fillna(False).astype(bool)) \
                .groupby("strain")["hit"].sum()
            print("  recall (ok + agrees with L1) per held-out strain: " + rec.astype(int).astype(str).add("/16").to_dict().__repr__())
    finally:
        shutil.rmtree(tmp)


if __name__ == "__main__":
    main()
