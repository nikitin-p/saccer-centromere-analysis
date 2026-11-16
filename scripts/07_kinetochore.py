#!/usr/bin/env python3
"""Analysis C: kinetochore genes in each assembly (+ S288C as a control).

Outputs in config `output.dir`:
  kinetochore_hsps.tsv        raw tblastn HSPs
  kinetochore_genes.tsv       per genome x gene: locus, status, length, identity, haplotype
  kinetochore_variants.tsv    per genome x gene x variant vs the S288C protein
  kinetochore_genes.fasta     extracted proteins;  kinetochore_cds.fasta  their CDS
"""
import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from centromere.genes import (assign_haplotypes, extract_cds, frameshift_sites, gene_name, make_db,
                              pick_locus, protein_variants, run_tblastn, to_plus, variant_label)
from centromere.io import load_config, save_tsv


def genomes(config: dict) -> list[tuple[str, Path]]:
    kin = config["kinetochore"]
    samples = pd.read_csv(config["input"]["sample_map"], sep="\t")
    gdir = Path(config["cde"]["genome_dir"])
    out = [(r["strain"], gdir / f"{r['accession']}.fasta.gz") for _, r in samples.iterrows()]
    return out + [("S288C", Path(kin["s288c_genome"]).expanduser())]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    kin = config["kinetochore"]
    out = Path(config["output"]["dir"])
    blast_bin = kin.get("blast_bin")   # null -> PATH
    db_dir = out / "kinetochore_blastdb"
    ref = {gene_name(r.id): str(r.seq) for r in SeqIO.parse(kin["proteins"], "fasta")}

    all_hsps, rows, variants, faa, fna = [], [], [], [], []
    for strain, path in genomes(config):
        db = make_db(path, db_dir, blast_bin)
        hsps = run_tblastn(Path(kin["proteins"]), db, blast_bin, kin["evalue"])
        hsps.insert(0, "strain", strain)
        all_hsps.append(hsps)
        contigs = {r.id: str(r.seq) for r in SeqIO.parse(db, "fasta")}
        for gene in ref:
            h = hsps[hsps["qseqid"].map(gene_name) == gene]
            row = {"strain": strain, "gene": gene, "ref_length": len(ref[gene])}
            if h.empty:
                rows.append(row | {"status": "absent"})
                continue
            locus, n_loci = pick_locus(h)
            contig, strand = locus["sseqid"].iloc[0], locus["strand"].iloc[0]
            seq = contigs[contig] if strand == "+" else str(Seq(contigs[contig]).reverse_complement())
            cds = extract_cds(seq, locus)
            start, end = to_plus(cds["o_start"], cds["o_end"], strand, len(seq))
            var = protein_variants(ref[gene], cds["protein"])
            n_sub = sum(v["type"] == "substitution" for v in var)
            identity = 1 - n_sub / len(ref[gene]) if not any(v["type"] != "substitution" for v in var) else None
            rows.append(row | {
                "status": cds["status"], "contig": contig, "strand": strand, "start": start, "end": end,
                "protein_length": len(cds["protein"]), "start_shift_codons": cds["start_shift"],
                "hsps": len(locus), "frames": cds["frames"],
                "frameshifts": ";".join(f"q{f['query_pos']}(hp{f['homopolymer']})"
                                        for f in frameshift_sites(seq, locus)),
                "n_loci": n_loci, "n_substitutions": n_sub,
                "n_other_variants": sum(v["type"] != "substitution" for v in var),
                "identity": identity, "variants": ",".join(variant_label(v) for v in var),
                "protein": cds["protein"]})
            variants += [{"strain": strain, "gene": gene} | v | {"label": variant_label(v)} for v in var]
            faa.append(f">{strain}|{gene} {contig}:{start}-{end}({strand}) {cds['status']}\n{cds['protein']}\n")
            fna.append(f">{strain}|{gene} {contig}:{start}-{end}({strand})\n{cds['cds']}\n")

    genes = pd.DataFrame(rows)
    present = genes["status"] != "absent"
    genes.loc[present, "haplotype"] = assign_haplotypes(genes[present], ref)
    save_tsv(pd.concat(all_hsps), out / "kinetochore_hsps.tsv")
    save_tsv(genes.drop(columns="protein"), out / "kinetochore_genes.tsv")
    save_tsv(pd.DataFrame(variants), out / "kinetochore_variants.tsv")
    (out / "kinetochore_genes.fasta").write_text("".join(faa))
    (out / "kinetochore_cds.fasta").write_text("".join(fna))

    print(genes.pivot(index="strain", columns="gene", values="status").to_string())
    print()
    print(genes.pivot(index="strain", columns="gene", values="haplotype").to_string())


if __name__ == "__main__":
    main()
