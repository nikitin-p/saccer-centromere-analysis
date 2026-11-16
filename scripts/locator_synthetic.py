#!/usr/bin/env python3
"""Locator spec phase 4, checks 3-5: plant a damaged centromere in one real assembly, run all
four locators on the copy and reconcile them.

Cases (CEN in its S288C orientation; the rest of the genome untouched):
  none          unmodified copy (control: every CEN `agree`)
  cdeii_random  CDEII replaced by random 90 %-AT sequence of the same length (check 3)
  cde_mutated   as cdeii_random, plus 2 substitutions in CDEI and 3 in CDEIII (check 4)
  insertion     a full Ty1 (LTR-INT-LTR from data/ty_reference.fasta, ~5.9 kb) inserted in the
                middle of CDEII (check 5; longer than preprocess.merge_max_distance)

Per case: L1 = BLAST of the S288C CENs + the 01_preprocess steps + 02_filter; L2 = flank
megablast + `locate_by_flanks`; L3 = minimap2 (locate.preset) + liftover (bridged as in
12_flanks_diversity); L4 = the final L4 profiles (`<hmm_dir>/round_<k>/all.hmm`, not rebuilt).
Then `reconcile` and `rescue_rows` exactly as 01_preprocess.

Outputs in `locate.synthetic_dir`:
  <case>/<accession>.fasta, blast/minimap2/nhmmer intermediates
  locator_synthetic.tsv          case x locator for the target CEN: found, locus, status
  locator_synthetic_summary.tsv  case: category, support, final locus, CDEs vs the planted ones,
                                 categories of the 15 untouched CENs

  PYTHONPATH=src .venv/bin/python scripts/locator_synthetic.py --config config.scrap.yaml \
      [--strain ADS] [--cen CEN5] [--cases none,cdeii_random,cde_mutated,insertion]
"""
import argparse
import random
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from centromere.cde import CDEI_LEN, CDEIII_LEN, CDEModel, annotate_cde, read_reference
from centromere.filter import by_gc, by_length
from centromere.hmm import locate_in_genome, model_lengths, run_nhmmer
from centromere.io import COLUMNS, load_config, save_tsv, tool, tool_file
from centromere.locate import (FLANK_COLS, L4_FOUND, MAIN_TABLE, assemblies, lift_interval,
                               locate_by_flanks, parse_paf, read_fasta, reconcile, rescue_rows,
                               s288c_cens)
from centromere.microevo import revcomp
from centromere.preprocess import assign_loci, find_split_hits, merge_split_hits, remove_gaps

CASES = ["none", "cdeii_random", "cde_mutated", "insertion"]


def mutate(seq: str, n: int, rng: random.Random) -> str:
    s = list(seq)
    for i in rng.sample(range(len(s)), n):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def ty1_element() -> str:
    ty = read_fasta(Path(__file__).parent.parent / "data" / "ty_reference.fasta")
    ltr = next(s for k, s in ty.items() if k.startswith("LTR|Ty1|"))
    return ltr + ty["INT|Ty1|YDRCTy1-1"] + ltr


def planted_cen(case: str, cen: str, rng: random.Random) -> str:
    """The CEN (CEN orientation) after the modification."""
    c1, c2, c3 = cen[:CDEI_LEN], cen[CDEI_LEN:-CDEIII_LEN], cen[-CDEIII_LEN:]
    if case == "none":
        return cen
    if case == "insertion":
        h = len(c2) // 2
        return c1 + c2[:h] + ty1_element() + c2[h:] + c3
    at = "".join(rng.choice("AT") if rng.random() < 0.9 else rng.choice("GC") for _ in c2)
    if case == "cde_mutated":
        c1, c3 = mutate(c1, 2, rng), mutate(c3, 3, rng)
    return c1 + at + c3


def run(cmd, **kw):
    subprocess.run([str(c) for c in cmd], check=True, capture_output=True, **kw)


def blast(cfg, name):
    return tool(cfg["kinetochore"].get("blast_bin"), name)   # null -> PATH


def run_l1(fa, acc, strain, cdir, cfg, reference):
    db = cdir / "db"
    run([blast(cfg, "makeblastdb"), "-in", fa, "-dbtype", "nucl", "-out", db])
    res = subprocess.run([blast(cfg, "blastn"), "-query", cfg["reference"], "-db", db, "-evalue", "1e-5",
                          "-word_size", "7", "-outfmt", "6 qseqid qstart qend sseqid sstart send sseq evalue"],
                         check=True, capture_output=True, text=True).stdout
    rows = [[acc] + line.split("\t") for line in res.strip().splitlines()]
    df = pd.DataFrame(rows, columns=COLUMNS)
    for c in ["qstart", "qend", "sstart", "send"]:
        df[c] = df[c].astype(int)
    df["evalue"] = df["evalue"].astype(float)
    df = assign_loci(remove_gaps(df.assign(assembly=strain, strain=strain)))
    df["blast_pieces"] = df.groupby(["sample", "centromere", "sseqid"])["sseqid"].transform("size")
    split, ok = find_split_hits(df)
    if len(split):
        ok = pd.concat([ok, merge_split_hits(split, cfg["preprocess"]["merge_max_distance"])], ignore_index=True)
    df = annotate_cde(ok, cdir, reference, cfg["cde"]["flank"], cfg["cde"]["search_radius"])
    return by_gc(by_length(df, cfg["filter"]["min_length"]), cfg["filter"]["min_gc"]), db


def run_l2(db, acc, cfg):
    loc = cfg["locate"]
    res = subprocess.run([blast(cfg, "blastn"), "-task", "megablast", "-query", loc["flanks"], "-db", db,
                          "-evalue", "1e-20", "-outfmt",
                          "6 qseqid qstart qend qlen sseqid sstart send pident length evalue bitscore"],
                         check=True, capture_output=True, text=True).stdout
    h = pd.DataFrame([[acc] + l.split("\t") for l in res.strip().splitlines()], columns=FLANK_COLS)
    for c in ["qstart", "qend", "qlen", "sstart", "send", "length"]:
        h[c] = h[c].astype(int)
    for c in ["pident", "evalue", "bitscore"]:
        h[c] = h[c].astype(float)
    h[["centromere", "side"]] = h["qseqid"].str.split("|", expand=True)
    return locate_by_flanks(h, loc["min_flank_pident"], loc["min_flank_cov"], loc["max_locus_bp"])


def run_l3(fa, acc, cdir, cfg, cens):
    loc = cfg["locate"]
    paf = cdir / "aln.paf"
    with open(paf, "w") as fh:
        subprocess.run([tool_file(loc.get("minimap2"), "minimap2"), "-x", loc.get("preset", "asm5"), "-c",
                        "--cs", "-t", str(loc.get("threads", 4)),
                        str(Path(loc["paf_dir"]) / "sacCer3.nuclear.fa"), str(fa)],
                       check=True, stdout=fh, stderr=subprocess.DEVNULL)
    blocks = parse_paf(paf, loc.get("min_mapq", 20))
    rows = []
    for _, c in cens.iterrows():
        r, status = lift_interval(blocks, c["chrom"], c["start"], c["end"]), "lifted"
        if isinstance(r, str):
            r2 = lift_interval(blocks, c["chrom"], c["start"], c["end"], loc.get("bridge_bp", 1000))
            status, r = (r, r) if isinstance(r2, str) else ("bridged", r2)
        row = {"sample": acc, "centromere": c["centromere"], "status": status}
        if not isinstance(r, str):
            ctg, qs, qe, strand = r
            row.update(contig=ctg, start=qs + 1, end=qe, exact=status == "lifted",
                       strand="+" if (strand == "+") == (c["strand"] == "+") else "-")
        rows.append(row)
    return pd.DataFrame(rows)


def run_l4(fa, acc, cdir, cfg, model, ref_len, anchors):
    loc = cfg["locate"]
    rdir = sorted(Path(loc["hmm_dir"]).glob("round_*"), key=lambda p: int(p.name.split("_")[1]))[-1]
    allh = rdir / "all.hmm"
    hits = run_nhmmer(allh, fa, cdir / "l4.tbl", tool(loc.get("hmmer_bin"), "nhmmer"),
                      loc.get("threads", 4), loc.get("hmm_report_evalue", 1.0))
    return locate_in_genome(hits, read_fasta(fa), model_lengths(allh), model, ref_len, anchors,
                            loc.get("hmm_evalue", 1e-10), loc.get("hmm_margin", 10),
                            loc.get("hmm_len_tol", 40), cfg["cde"]["flank"],
                            cfg["cde"]["search_radius"]).assign(sample=acc)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--strain", default="ADS")
    ap.add_argument("--cen", default="CEN5")
    ap.add_argument("--cases", default=",".join(CASES))
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    cfg = load_config(args.config)
    loc = cfg["locate"]
    out = Path(loc.get("synthetic_dir", Path(cfg["output"]["dir"]) / "locator_synthetic"))
    out.mkdir(parents=True, exist_ok=True)
    reference = read_reference(cfg["reference"])
    ref_len = {c: len(s) for c, s in reference.items()}
    model = CDEModel(reference)
    cens = s288c_cens(cfg)

    asm = assemblies(cfg).set_index("strain")
    if args.strain not in asm.index:
        print(f"strain {args.strain} not in {cfg['input']['sample_map']}; planting in {asm.index[0]} instead")
        args.strain = asm.index[0]
    a = asm.loc[args.strain]
    acc = a["accession"]
    main_tab = pd.read_csv(Path(cfg["output"]["dir"]) / "preprocessed.tsv", sep="\t")
    t = main_tab[(main_tab["sample"] == acc) & (main_tab["centromere"] == args.cen)].iloc[0]
    genome = read_fasta(a["fasta"])
    ctg = genome[t["sseqid"]]
    lo, hi, plus = int(t["cen_start"]), int(t["cen_end"]), t["strand"] == "+"
    orig = ctg[lo - 1:hi] if plus else revcomp(ctg[lo - 1:hi])
    assert orig == t["sequence"], "main-table locus does not match the assembly"

    detail, summary = [], []
    for case in args.cases.split(","):
        rng = random.Random(args.seed)
        cdir = out / case
        if cdir.exists():
            shutil.rmtree(cdir)
        cdir.mkdir(parents=True)
        new = planted_cen(case, orig, rng)
        g = dict(genome)
        g[t["sseqid"]] = ctg[:lo - 1] + (new if plus else revcomp(new)) + ctg[hi:]
        fa = cdir / f"{acc}.fasta"
        with open(fa, "w") as fh:
            for k, s in g.items():
                fh.write(f">{k}\n{s}\n")
        print(f"[{case}] planted {args.strain} {args.cen}: {len(orig)} -> {len(new)} bp")

        l1, db = run_l1(fa, acc, args.strain, cdir, cfg, reference)
        l2 = run_l2(db, acc, cfg)
        l3 = run_l3(fa, acc, cdir, cfg, cens)
        l2ok = l2[l2["status"] == "ok"].rename(columns={"l2_start": "start", "l2_end": "end"})
        l3ok = l3[l3["status"].isin(["lifted", "bridged"])]
        cols = ["sample", "centromere", "contig", "start", "end"]
        l4 = run_l4(fa, acc, cdir, cfg, model, ref_len, pd.concat([l2ok[cols], l3ok[cols]]))
        std = ["sample", "centromere", "contig", "start", "end", "strand", "exact"]
        loci = {"l1": l1.rename(columns={"sseqid": "contig", "cen_start": "start", "cen_end": "end"})
                .assign(exact=True)[std],
                "l2": l2ok.assign(exact=True)[std], "l3": l3ok[std],
                "l4": l4[l4["status"].isin(L4_FOUND)].rename(columns={"hmm_start": "start", "hmm_end": "end"})
                .assign(exact=True)[std]}
        rec = reconcile(loci, {k: {acc} for k in loci}, ref_len, loc.get("insertion_bp", 30))
        resc = rescue_rows(rec, cdir, reference, {acc: args.strain}, cfg["cde"]["flank"], cfg["cde"]["search_radius"])

        status = {"l1": "found" if (l1["centromere"] == args.cen).any() else "missed",
                  "l2": l2.set_index("centromere").loc[args.cen, "status"],
                  "l3": l3.set_index("centromere").loc[args.cen, "status"],
                  "l4": l4.set_index("centromere").loc[args.cen, "status"]}
        for tag, d in loci.items():
            x = d[d["centromere"] == args.cen]
            detail.append({"case": case, "locator": tag.upper(), "status": status[tag], "found": len(x) > 0,
                           **({"contig": x["contig"].iat[0], "start": int(x["start"].min()),
                               "end": int(x["end"].max()), "length": int(x["end"].max() - x["start"].min() + 1)}
                              if len(x) else {})})
        r = rec.set_index("centromere").loc[args.cen]
        if r["category"] in MAIN_TABLE:
            fin = resc.set_index("centromere").loc[args.cen]
        else:
            fin = l1.set_index("centromere").loc[args.cen] if (l1["centromere"] == args.cen).any() else None
        others = rec[rec["centromere"] != args.cen]["category"].value_counts().to_dict()
        row = {"case": case, "planted_length": len(new), "category": r["category"], "support": r["support"],
               "source": r.get("source"), "insertion_bp": r.get("insertion_bp"),
               "others": ", ".join(f"{k} {v}" for k, v in others.items())}
        if fin is not None:
            row.update(final_length=int(fin["sequence_length"]), final_cdeii=int(fin["cdeii_length"]),
                       cde_ok=bool(fin["cde_ok"]), sequence_is_planted=fin["sequence"] == new,
                       cdei_is_planted=fin["cdei_seq"] == new[:CDEI_LEN],
                       cdeiii_is_planted=fin["cdeiii_seq"] == new[-CDEIII_LEN:])
        summary.append(row)
        print(f"  {status} -> {r['category']} [{r['support']}]")

    det, summ = pd.DataFrame(detail), pd.DataFrame(summary)
    save_tsv(det, out / "locator_synthetic.tsv")
    save_tsv(summ, out / "locator_synthetic_summary.tsv")
    print(det.to_string(index=False))
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
