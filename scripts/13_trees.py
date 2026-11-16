#!/usr/bin/env python3
"""Analysis H: centromere tree vs genome tree (10 ScRAP assemblies + S288C).

Spec docs/specs/2026-11-13-flanks-diversity-trees.md, step 2. Needs aligned.tsv and the PAFs
from scripts/run_minimap2.sh. Outputs in `<output.dir>/trees/`:
  cen_concat.fasta          16 per-CEN MSAs concatenated (S288C = ref_aligned_seq)
  cen_blocks.tsv            column range of each CEN in the concatenation
  genome_snps.tsv.gz        variable single-copy core sites (S288C coordinates, 1-based)
  genome_core_stats.tsv     core length, single-copy coverage and SNPs vs S288C per strain
  dist_cen.tsv, dist_genome.tsv   K80 distance matrices
  strain_discordance.tsv    1 - Spearman rho of each strain's rows in the two matrices
  cen_contribution.tsv      per strain x CEN: discordance change when that CEN is dropped
  swap_check.tsv            spec check 8 on the real concatenation (configured pair)
  swap_scan.tsv             the same for every ordered strain/donor pair (power of the read-out)
  cen_nj.nwk, genome_nj.nwk, tree_compare.tsv   from scripts/trees_nj.R (NJ, RF, Mantel)
"""
import argparse
import gzip
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from centromere.io import load_config, save_tsv
from centromere.locate import CODE_BASE, assemblies, parse_paf, read_fasta, ref_alleles, s288c_cens
from centromere.trees import (cen_contribution, cen_matrix, concat_cen, discordance, dist_frame, k80,
                              default_swap_pair, swap_check, swap_scan)


def genome_matrix(asm, ref, cens, min_mapq, exclusion):
    """Allele codes of every assembly at single-copy core sites (all assemblies covered by
    exactly one primary block, ACGT, S288C ref ACGT), outside CEN +- exclusion."""
    blocks = {a: parse_paf(p, min_mapq) for a, p in zip(asm["accession"], asm["paf"])}
    snp_rows, stats = [], []
    cover = {s: 0 for s in asm["strain"]}
    n_core = 0
    for chrom, seq in ref.items():
        codes = np.stack([ref_alleles(blocks[a], chrom, seq) for a in asm["accession"]])
        for s, row in zip(asm["strain"], codes):
            cover[s] += int((row > 0).sum())
        r = np.zeros(len(seq), dtype=np.int8)
        for i, b in enumerate("ACGT", 1):
            r[np.frombuffer(seq.encode(), dtype=np.uint8) == ord(b)] = i
        core = (r > 0) & ((codes >= 1) & (codes <= 4)).all(axis=0)
        for _, c in cens[cens["chrom"] == chrom].iterrows():
            core[max(0, c["start"] - exclusion):c["end"] + exclusion] = False
        var = core & (codes != r).any(axis=0)
        n_core += int(core.sum())
        stats.append({"chrom": chrom, "length": len(seq), "core": int(core.sum()), "variable": int(var.sum())})
        idx = np.flatnonzero(var)
        snp_rows.append(pd.DataFrame({"chrom": chrom, "pos": idx + 1, "S288C": CODE_BASE[r[idx]],
                                      **{s: CODE_BASE[codes[k, idx]] for k, s in enumerate(asm["strain"])}}))
    snps = pd.concat(snp_rows, ignore_index=True)
    st = pd.DataFrame(stats)
    total = st["length"].sum()
    per = pd.DataFrame({"strain": list(cover), "single_copy_frac": [cover[s] / total for s in cover],
                        "snps_vs_s288c": [int((snps[s] != snps["S288C"]).sum()) for s in cover]})
    return snps, n_core, st, per


def write_dist(d: pd.DataFrame, path: Path):
    d.rename_axis("taxon").reset_index().to_csv(path, sep="\t", index=False, float_format="%.6g")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    cfg = load_config(ap.parse_args().config)
    out = Path(cfg["output"]["dir"])
    tdir = out / "trees"
    tdir.mkdir(parents=True, exist_ok=True)
    tc = cfg["trees"]
    order = cfg["visualize"]["centromere_order"]

    aligned = pd.read_csv(out / "aligned.tsv", sep="\t")
    seqs, blocks = concat_cen(aligned, order)
    (tdir / "cen_concat.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    save_tsv(blocks, tdir / "cen_blocks.tsv")
    taxa = list(seqs)
    d_cen = cen_matrix(seqs, blocks)
    print(f"CEN concatenation: {len(taxa)} taxa x {blocks['end'].max()} columns")

    asm = assemblies(cfg)
    ref = read_fasta(Path(cfg["locate"]["paf_dir"]) / "sacCer3.nuclear.fa")
    cens = s288c_cens(cfg)
    snps, n_core, chrom_stats, strain_stats = genome_matrix(asm, ref, cens, cfg["locate"]["min_mapq"],
                                                            tc["cen_exclusion"])
    with gzip.open(tdir / "genome_snps.tsv.gz", "wt") as fh:
        snps.to_csv(fh, sep="\t", index=False)
    core_stats = pd.concat([chrom_stats.assign(strain=""), strain_stats.assign(chrom="")], ignore_index=True)
    save_tsv(core_stats, tdir / "genome_core_stats.tsv")
    print(f"genome core: {n_core:,} single-copy sites in all {len(asm)} assemblies "
          f"({n_core / chrom_stats['length'].sum():.1%} of R64 nuclear); {len(snps):,} variable")
    print(strain_stats.to_string(index=False))

    codes = np.stack([np.searchsorted(list("ACGT"), snps[t].to_numpy()) + 1 for t in taxa])
    d_gen = dist_frame(k80(codes.astype(np.int8), extra_sites=n_core - len(snps))[0], taxa)
    write_dist(d_cen, tdir / "dist_cen.tsv")
    write_dist(d_gen, tdir / "dist_genome.tsv")

    disc = discordance(d_cen, d_gen).sort_values(ascending=False)
    save_tsv(disc.rename_axis("strain").reset_index(), tdir / "strain_discordance.tsv")
    contrib = cen_contribution(seqs, blocks, d_gen)
    save_tsv(contrib, tdir / "cen_contribution.tsv")
    print("strain discordance (1 - Spearman rho):")
    print(disc.round(3).to_string())
    top = contrib.sort_values("contribution", ascending=False).groupby("strain").head(2)
    print("top-2 contributing CENs for the 3 most discordant strains:")
    print(top[top["strain"].isin(disc.index[:3])].round(3).to_string(index=False))

    strain, donor = tc.get("swap_strain"), tc.get("swap_donor")
    if strain not in seqs or donor not in seqs:
        strain, donor = default_swap_pair(d_cen, d_gen)
        print(f"swap pair {tc.get('swap_strain')}<-{tc.get('swap_donor')} not in the data; using "
              f"{strain}<-{donor} (most concordant strain, most genome-distant donor)")
    sw = swap_check(seqs, blocks, d_gen, strain, donor, tc["swap_n"])
    save_tsv(pd.DataFrame([{k: (",".join(v) if isinstance(v, list) else v) for k, v in sw.items()
                            if k != "discordance"}]), tdir / "swap_check.tsv")
    print(f"swap check: {sw['strain']} <- {sw['donor']} at {sw['swapped']}: most discordant = "
          f"{sw['most_discordant']}, top contributors = {sw['top_contributors']}, detected = {sw['detected']}")

    scan = swap_scan(seqs, blocks, d_gen, tc["swap_n"], tc.get("swap_scan_max_pairs"), tc.get("seed", 1))
    save_tsv(scan, tdir / "swap_scan.tsv")
    print(f"swap scan over {len(scan)} pairs: swapped strain most discordant {scan['most_discordant'].mean():.0%}, "
          f"all {tc['swap_n']} swapped CENs on top {(scan['n_swapped_in_top'] == tc['swap_n']).mean():.0%}, "
          f"both (detected) {scan['detected'].mean():.0%}; "
          f"median swapped CENs in top {tc['swap_n']}: {scan['n_swapped_in_top'].median():.0f}")

    subprocess.run([tc.get("rscript", "Rscript"), "scripts/trees_nj.R", str(tdir)], check=True)
    print(pd.read_csv(tdir / "tree_compare.tsv", sep="\t").to_string(index=False))


if __name__ == "__main__":
    main()
