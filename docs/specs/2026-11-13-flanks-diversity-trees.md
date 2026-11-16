# Spec — L3 liftover, analysis F (π CEN vs flanks), analysis H (CEN vs genome tree) (2026-11-13)

Implements locator spec Phase 2 (L3 only) as the base, then sections F1 and H of
`2026-10-20-centromere-microevolution.md`, on the 10 ScRAP assemblies + S288C (n = 11).
Goal: (F) is nucleotide diversity at CDEI/CDEII/CDEIII lower, equal or higher than in the
±2 kb around it; (H) do the 16 centromeres tell the same strain history as the genome.

Choices taken (2026-11-13): L3 only (no L2, L4, reconcile); minimap2 2.30 in conda env
`~/miniconda3/envs/minimap2`; trees by NJ in R/ape (no IQ-TREE). F2 (CDEII grammar /
DNA shape) not in this round.

## Inputs
- Assemblies `data/assemblies/*.fasta.gz`, map `data/scrap_assemblies.tsv`.
- S288C R64 `~/Desktop/pyselectal-analysis/data/sacCer3.fa` (chrM excluded from
  everything), `data/s288c_cen.bed` (16 CENs, SGD coordinates, strand).
- `output/scrap/aligned.tsv`: per-CEN MSA (`aligned_seq` + `ref_aligned_seq` = S288C),
  CDE sequences/boundaries per locus; `filtered.tsv` loci (contig, cen_start/end, strand).

## Step 0 — L3: minimap2 + own liftover (`src/centromere/locate.py`)
1. `scripts/run_minimap2.sh <config>`: per assembly
   `minimap2 -x asm5 -c --cs -t 4 sacCer3.fa <asm> > <out>/minimap2/<accession>.paf`;
   skip if the PAF exists and is non-empty (same resume behaviour as `blast.sh`).
   ~3 s per genome. Preset from config (`asm5`; `asm10`/`asm20` selectable).
2. `parse_paf` → blocks (tp:A:P primary only, MAPQ ≥ `min_mapq`=20), with the `cs` string.
3. `liftover_paf(blocks, chrom, pos)` → (contig, query pos, strand) or reason
   (`unaligned`, `in_insertion_gap`, `multiple`), by walking `cs` (`:n`, `*xy`, `+ins`,
   `-del`); reverse-strand blocks mapped from query end. *Rejected:* `paftools.js
   liftover` (needs k8).
4. `lift_interval(chrom, start, end)` → query interval if both ends lift to the same
   contig and strand, else the reason (`split_contigs`, `one_end`, `strand_flip`).
5. `l3_loci.tsv`: liftover of each S288C CEN → compared with the L1 locus (same contig,
   overlap, Δstart/Δend). Locator check 2 (L1 vs L3 part): expect 160/160 agree;
   disagreements listed.

## Step 1 — flank windows (F)
1. Window = S288C CEN ± `flank_bp` (2000) lifted per assembly (step 0.4); a window that
   does not lift intact is shortened to the largest lifted span containing the CEN, and
   the row records how much of each flank was kept. Extracted in S288C CEN orientation
   → `output/scrap/flanks/<CEN>_windows.fasta` (11 sequences incl. S288C).
2. MAFFT (`--auto`, as `03_align.py`) per CEN → `flanks/<CEN>_windows.aln.fasta`.
3. Coordinates: every column gets an S288C coordinate relative to the CEN (columns
   that are gaps in S288C are assigned to the preceding S288C position). Regions:
   `CDEI`, `CDEII`, `CDEIII` from the S288C CDE boundaries (`aligned.tsv` /
   `cde.py`), `left_k` / `right_k` 100-bp windows (k = 1…20, counted outward from the CEN).
4. **π per region** = mean over sequence pairs of (differences / comparable sites),
   pairwise deletion of gaps; also `π_sub` (substitutions only) and `indel_rate`
   (gap-open events per comparable site, so CDEII length changes are visible
   separately from substitutions). Regions with < 50 % of pairs comparable → NA.
   → `diversity_windows.tsv` (CEN, region, offset, n_sites, n_pairs, π, π_sub, indel_rate).
5. **Comparison** per CEN: π of each CDE vs the distribution of its 40 flank windows
   (percentile); pooled over CENs: median π per element vs flanks, Wilcoxon signed-rank
   across the 16 CENs (CDE π vs mean flank π of that CEN) — 16 paired values, the only
   test reported. → `diversity_summary.tsv`.
6. Cross-check: π of the CEN region from the window MSA vs from `aligned.tsv` MSA
   (same 11 sequences) — must agree within 10 % or the difference is explained.

## Step 2 — trees (H)
1. **CEN alignment**: concatenate the 16 per-CEN MSAs from `aligned.tsv`
   (S288C row = `ref_aligned_seq`) → `trees/cen_concat.fasta` (~1.9 kb, 11 taxa).
2. **Genome SNP matrix** from the PAFs: for every S288C position on chrI–chrXVI, the
   allele each assembly carries (from `cs`), only where the position is covered by exactly
   one primary block (MAPQ ≥ 20) in that assembly. Core = positions covered in all 10;
   SNP-only (indels not used). Also positions inside CEN ± 2 kb are excluded so the two
   trees share no sites. → `trees/genome_snps.tsv.gz` (variable core sites only) and
   `trees/genome_core_stats.tsv` (core length, SNPs per strain).
3. **Distances** (Python, same code for both): K80 from transitions/transversions per pair
   on comparable sites (gaps / non-ACGT pairwise-deleted). → `trees/dist_cen.tsv`,
   `trees/dist_genome.tsv`.
4. **R (report helpers, ape + vegan, both installed):** `ape::nj` on each matrix,
   midpoint-rooted; `ape::dist.topo` (RF, also normalised by 2(n−3)); Mantel
   (`vegan::mantel`, Spearman, 9999 permutations); tanglegram `ape::cophyloplot`.
   `13_trees.py` calls `Rscript scripts/trees_nj.R` once, which writes
   `trees/cen_nj.nwk`, `trees/genome_nj.nwk` and `trees/tree_compare.tsv` (RF, nRF,
   Mantel r, p), so the results exist without rendering the report.
5. **Strain discordance** = 1 − Spearman ρ between a strain's row in the CEN and genome
   matrices. **Per-CEN contribution** = change in that strain's discordance when CEN k is
   dropped from the concatenation (positive = CEN k drives the discordance).
   → `trees/strain_discordance.tsv`, `trees/cen_contribution.tsv`.
6. **Swap check** (spec check 8): copy the CEN concatenation, give strain X strain Y's
   block for 4 CENs, recompute → X must become the most discordant strain and those 4
   CENs must top its contribution list. Run in the script (X, Y fixed in config) and in a test.

## Files
| File | Change |
|---|---|
| `scripts/run_minimap2.sh` | new: step 0.1 |
| `src/centromere/locate.py` | new: `parse_cs`, `parse_paf`, `liftover_paf`, `lift_interval`, `ref_alleles` (per-position alleles from cs) |
| `src/centromere/microevo.py` | add: `region_map`, `pairwise_pi`, `diversity_windows` |
| `src/centromere/trees.py` | new: concatenation, K80 distances, discordance, per-CEN contribution, swap simulation |
| `scripts/12_flanks_diversity.py` | new: steps 0.2–0.5 and 1 |
| `scripts/13_trees.py` | new: step 2 (calls `Rscript scripts/trees_nj.R`) |
| `scripts/trees_nj.R` | new: NJ, RF, Mantel, Newick |
| `config.scrap.yaml` | new section `locate:` (`minimap2`, `preset: asm5`, `min_mapq: 20`, `paf_dir`), `diversity:` (`flank_bp: 2000`, `window: 100`), `trees:` (`swap_strain`, `swap_donor`, `swap_n`) |
| `tests/test_locate.py`, `tests/test_trees.py`, `tests/test_microevo.py` (add) | below |
| `results/report.qmd`, `results/report_helpers.R` | sections `#sec-diversity` (F) and `#sec-trees` (H): π profile figure (offset vs π, CDEs marked, 16 CENs + median), element table; tanglegram, RF/Mantel table, discordance table |

## Decisions (rejected alternative)
- π from MAFFT of lifted windows, not directly from the minimap2 projection
  (rejected: projection drops insertions relative to S288C and asm5 may not align
  through AT-rich CDEII; the MSA keeps CDEII indels, which matter for F).
- Genome tree from SNPs on single-copy core sites (rejected: including multi-copy /
  partially covered positions — paralogous Ty and subtelomeric sites would add noise).
- CEN ± 2 kb removed from the genome matrix (rejected: keeping them — the trees would
  share sites and Mantel would be inflated).
- One test for F (Wilcoxon over 16 CENs), H descriptive + Mantel (rejected: per-CEN
  p-values at n = 11).
- K80 NJ in ape, no new packages (rejected: phangorn/phytools — not installed; ape
  `dist.topo` and `cophyloplot` cover RF and the tanglegram).
- Assemblies taken as one haplotype each, as in all earlier analyses (collapsed ones
  included); stated in the report.

## Checks
1. `tests/test_locate.py`: `cs` walk on hand-made blocks (match, SNP, insertion,
   deletion, reverse strand) returns exact coordinates; position in an insertion gap /
   unaligned gives the reason; real-data test: liftover of S288C CEN on sacCer3 aligned to
   itself returns identity.
2. Step 0.5: L1 vs L3 loci 160/160 (printed); disagreements listed.
3. `tests/test_microevo.py` (spec check 6): π on a synthetic alignment with known
   differences equals the exact value; gaps handled by pairwise deletion as specified.
4. `tests/test_trees.py`: K80 on a hand-computed pair; swap simulation detects the
   swapped strain and its CENs (spec check 8); also printed by `13_trees.py` on real data.
5. Step 1.6 cross-check printed.
6. `.venv/bin/pytest -q tests` (2 kaleido tests may fail as before); report renders;
   both new sections checked by screenshot.

No git writes; you stage and commit.
