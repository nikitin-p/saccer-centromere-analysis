#!/usr/bin/env python3
"""Analysis C.3: kinetochore protein alignments, haplotypes and variable positions by domain.

Reads the outputs of 07_kinetochore.py. Only `intact` alleles enter the alignments: a
frameshifted allele translates out of frame up to its stop and would add spurious columns.
The S288C reference protein is the coordinate system.

Outputs in config `output.dir`:
  kinetochore_aln/<GENE>.faa       MAFFT protein alignment (S288C reference + intact alleles)
  kinetochore_positions.tsv        variable alignment columns with reference position + domain
  kinetochore_domain_summary.tsv   per gene x domain: variable positions, per 100 aa
  kinetochore_haplotypes.tsv       per gene x haplotype: strains, status, variants, domains hit
"""
import argparse
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from centromere.align import run_mafft
from centromere.genes import domain_summary, domains_at, gene_name, msa_positions
from centromere.io import load_config, save_tsv

REF = "S288C_ref"


def label_pos(label: str) -> int | None:
    """First reference position in a variant label (`C262R`, `del50-58`, `ins129_ND`, ...)."""
    m = re.search(r"(\d+)", label)
    return int(m.group(1)) if m else None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    kin = config["kinetochore"]
    out = Path(config["output"]["dir"])
    mafft_args = config["align"]["mafft_args"]
    ref = {gene_name(r.id): str(r.seq) for r in SeqIO.parse(kin["proteins"], "fasta")}
    domains = pd.read_csv(kin["domains"], sep="\t")
    genes = pd.read_csv(out / "kinetochore_genes.tsv", sep="\t", keep_default_na=False)
    genes = genes[(genes["status"] != "absent") & (genes["strain"] != "S288C")]
    proteins = {(r.id.split("|")[0], r.id.split("|")[1]): str(r.seq)
                for r in SeqIO.parse(out / "kinetochore_genes.fasta", "fasta")}

    aln_dir = out / "kinetochore_aln"
    aln_dir.mkdir(parents=True, exist_ok=True)
    positions, summaries, haplos = [], [], []
    for gene, g in genes.groupby("gene"):
        dom = domains[domains["gene"] == gene]
        intact = g[g["status"] == "intact"]
        names = [REF] + intact["strain"].tolist()
        seqs = [ref[gene]] + [proteins[(s, gene)] for s in intact["strain"]]
        aln = dict(zip(names, run_mafft(seqs, mafft_args)))
        (aln_dir / f"{gene}.faa").write_text("".join(f">{k}\n{v}\n" for k, v in aln.items()))

        p = msa_positions(aln, REF)
        p.insert(0, "gene", gene)
        p["domain"] = [domains_at(x, dom) for x in p["ref_pos"]]
        positions.append(p)
        s = domain_summary(p, dom, len(ref[gene]))
        s.insert(0, "gene", gene)
        s.insert(2, "n_intact", len(intact))
        summaries.append(s)

        for hap, h in g.groupby("haplotype", sort=False):
            labels = [v for v in h["variants"].iloc[0].split(",") if v]
            hit = [domains_at(label_pos(v), dom) for v in labels]
            haplos.append({"gene": gene, "haplotype": hap, "n_strains": len(h),
                           "strains": ",".join(h["strain"]), "status": ",".join(sorted(set(h["status"]))),
                           "n_variants": len(labels), "variants": ",".join(labels),
                           "domains_hit": ";".join(sorted({d for x in hit for d in x.split(";")} - {"-"}))})

    positions = pd.concat(positions, ignore_index=True)
    summaries = pd.concat(summaries, ignore_index=True)
    haplos = pd.DataFrame(haplos).sort_values(["gene", "n_strains"], ascending=[True, False])
    save_tsv(positions, out / "kinetochore_positions.tsv")
    save_tsv(summaries, out / "kinetochore_domain_summary.tsv")
    save_tsv(haplos, out / "kinetochore_haplotypes.tsv")

    pd.set_option("display.width", 200)
    print("S288C-private positions (all intact strains carry the non-reference allele):")
    print(positions[positions["ref_private"]][["gene", "ref_pos", "ref", "alleles", "domain"]].to_string(index=False))
    print()
    named = summaries[summaries["gene"].isin(["CBF1", "CEP3", "CSE4"])]
    print(named.to_string(index=False))
    print()
    per_gene = positions[~positions["insertion"]].groupby("gene").agg(
        variable=("ref_pos", "nunique"), substituted=("n_sub", lambda x: int((x > 0).sum())),
        ref_private=("ref_private", "sum"))
    per_gene["length"] = [len(ref[x]) for x in per_gene.index]
    print(per_gene.to_string())


if __name__ == "__main__":
    main()
