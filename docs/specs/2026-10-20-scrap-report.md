# Spec — ScRAP centromere report rebuilt from real data (2026-10-20)

## Goal
Rewrite `results/report.qmd` so every figure is computed at render time from the
ScRAP pipeline outputs (now 10 T2T genomes, later all 142), drawn in ggplot2 /
plotly. Keep the report's explanatory voice; drop every figure from the old run.

## Data flow
```
data/assemblies/*.fasta.gz ─ blast.sh ─► output/scrap/centromeres.tsv
  ─ steps 01–03 (config.scrap.yaml) ─► output/scrap/{preprocessed,filtered,aligned}.tsv
  ─ quarto render results/report.qmd (R engine) ─► results/report.html
                                                  + results/figures/*.png|pdf
```
The report reads only `output/scrap/*.tsv`, `data/scrap_assemblies.tsv`,
`data/41586_2018_30_MOESM3_ESM.xls`. Data dir is a Quarto
param (`data_dir`, default `../output/scrap`) so the 142-genome cluster output
renders with `-P data_dir:...`.

## Changes

| File | Change |
|---|---|
| `config.scrap.yaml` (new) | input `output/scrap/centromeres.tsv`, output `output/scrap`, `sample_map`, `reference: data/sacer_centromere.fasta`, `cde.genome_dir: data/assemblies` |
| `src/centromere/cde.py` (new) | **per-strain CDE annotation** (method below); called from `01_preprocess.py` when `cde.genome_dir` is set |
| `src/centromere/preprocess.py` | fix: a split piece > `merge_max_distance` away no longer overwrites the merged row's end coordinate |
| `src/centromere/align.py`, `scripts/03_align.py` | optional `reference`: the S288C CEN is aligned together with each group; its aligned row is stored as a column `ref_aligned_seq` on every row (not as an extra sample, so steps 04/05 are unaffected) |
| `results/report.qmd` | rewritten, R/knitr engine; all numbers inline from data |
| `results/figures/` | 8 old PNGs deleted; new figures written by the render via `ggsave` (PNG 300 dpi + PDF) |
| `tests/test_align.py` | test for `ref_aligned_seq` |

## Figures

| # | Figure | Form | Interactive |
|---|---|---|---|
| 1 | Strains: assembly, strain, clade, origin, ploidy (join on `strain`) | table | – |
| 2 | Genome → table funnel: raw hits, cross-centromere hits dropped, split hits merged, final rows, per strain | table + compact bar | – |
| 3 | Length vs S288C: strain × CEN tile, fill = Δ bp from S288C (diverging, 0 = white), label = length | heatmap | plotly hover |
| 4 | GC per CEN: one point per assembly coloured by clade, S288C as a reference tick; switches automatically to a ggdist raincloud when n ≥ 20 per CEN | strip / raincloud | plotly |
| 5 | **16-row logo panel**, CDEI and CDEIII boxed, CDEII labelled (the coursework's Fig 13, rebuilt with ggseqlogo + facets) | logo | static |
| 6 | **Variant map**: per CEN, rows = assemblies, x = alignment column, only differences from S288C drawn (SNV coloured by alt base, insertion, deletion), CDE bands on top | tile, 16 facets | static + plotly for one CEN selector |
| 7 | Variant density by element: SNVs and indel events per bp in CDEI / CDEII / CDEIII, per CEN | dot plot | plotly |
| 8 | Substitution spectrum vs S288C: 12 classes stacked by CEN, plus a panel normalised by S288C base opportunity (per-site rate) | bar ×2 | plotly |

Split/duplicated hits: a table (sample, CEN, pieces, gap bp) instead of the old
repeat histogram — with 5 genomes there is one split hit.

## Decisions (and the rejected alternative)
- **Compare to S288C, not to the population consensus.** At n = 5 a majority
  consensus is unstable and ties are common; S288C is fixed and is the query
  BLAST used. Still not ancestral — the report keeps the "not polarised" caveat.
  *Rejected:* consensus (current step 04).
- **CDE boundaries called in every strain, not inherited from BLAST.** BLAST only
  locates the locus. The window hit ± 150 bp is cut from the assembly and oriented
  like the S288C CEN. PWMs are built from the 16 S288C elements (SGD convention:
  CDEI = first 8 bp, CDEIII = last 26 bp; `CCGAA` sits 13 bp from the end). CDEI
  is the best PWM hit within ± 20 bp of where the BLAST hit projects the CEN
  start, CDEIII likewise around the projected end; CDEII = everything between,
  so insertions are kept. `sequence` becomes CDEI start → CDEIII end
  (`blast_sequence` keeps the old one); new columns `strand`, `contig_name`,
  `cen_start/cen_end`, `cdei/cdeii/cdeiii_seq`, `cdeii_length`, PWM scores and a
  `cde_ok` flag (both scores ≥ the weakest S288C element − 2 bits).
  *Rejected:* strict regex (misses S288C's own CEN2 and CEN9); CDE coordinates
  only on S288C (gives no per-strain boundaries).
- **Only the variant map and logo carry per-position detail.** No entropy or
  "variable positions" bar: at n = 5 it saturates trivially, which the current
  report already criticised. *Rejected:* keeping Fig 2.
- **Colour.** Clade colours for strains (≤ 26 clades; at n = 5 only 5 used);
  centromeres are never colour-coded as 16 hues in one legend — stacked bars use
  a sequential ordering with direct labels or facets instead. Built following the
  `dataviz` skill.
- **Dropped figures:** global all-CEN logo (MAFFT cannot co-align the 16 CENs —
  already stated in the old report); BLAST/MUMmer/HMMER comparison (reference
  genome, not ScRAP).
- **Dependencies:** install CRAN `plotly`, `ggdist`, `readxl` (if missing).
  BLAST locally from the NCBI binary in the session scratchpad; on the cluster
  from conda as before.

## Check that proves it works
1. `pytest -q tests` — all pass, including the new align test.
2. Steps 01–03 on `config.scrap.yaml` → `aligned.tsv` has 160 rows and a
   `ref_aligned_seq` column with equal length to `aligned_seq` per CEN.
3. `quarto render results/report.qmd` exits 0; `results/figures/` contains only
   new files; no reference to the old figure names remains in the qmd.
4. Spot-checks printed by a hidden chunk: rows = 160, assemblies = 10, all `cde_ok`, CEN16 split
   in GCA_949125725.1 listed in the split table, CDE boxes land on `TCACRTG` /
   `CCGAA` columns for all 16 CENs.
5. I look at every rendered figure image before reporting.

No git writes; you stage and commit.
