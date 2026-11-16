# Spec — F2: CDEII grammar and DNA shape (2026-11-15)

## Goal
Test whether CDEII has an A-tract "grammar" (tract lengths, ~10-bp phasing) beyond what
its AT richness gives, and describe predicted DNA shape along it. Then ask whether natural
CDEII variants preserve shape more than random mutations of the same kind would.

Exploratory item F2 of `docs/specs/2026-10-20-centromere-microevolution.md`; follows F.

## Inputs
- ScRAP CDEII per locus: `output/scrap/preprocessed.tsv` (`cdeii_seq`), 176 loci.
- Helsen SD4 distinct CDEII alleles: `helsen_cdeii()` on `data/external/helsen2026/`
  (3048 alleles); divergent alleles (identity < 0.85 in `helsen_cdeii_steps.tsv`) dropped.
- Population modal allele per CEN: `output/scrap/cdeii_reference_alleles.tsv` (Helsen-based
  baseline, as in B; S288C is not the baseline).
- S288C genome (`sacCer3.fa`, path from config) for AT-matched null windows.

## Definitions
- **A-tract**: maximal run of A/T, length ≥ 4, containing no TpA step (AnTm form;
  `AAAATTTT` = one tract of 8, `TTTTAAAA` = two tracts of 4). Min length in config.
- **Phasing score**: fraction of tract-centre pairs whose distance is within ±1 bp of a
  multiple of 10.5 bp (10, 21, 31 …), up to 42 bp. Second measure: power at period
  10–11 bp in the AA/TT dinucleotide autocorrelation along the sequence.
- **Shape**: DNAshapeR pentamer model — MGW and ProT per bp, Roll and HelT per step.
- **Shape distance** between two equal-length alleles: mean |Δfeature| over positions,
  per feature.

## Analyses
### F2.1 Grammar vs composition
Per CDEII (S288C 16; ScRAP 176; Helsen modal allele per CEN): number of tracts, tract-length
distribution, phasing score, 10–11 bp power.
Nulls, per sequence:
- (a) 200 dinucleotide shuffles (`ltr_cde.dinuc_shuffle`);
- (b) 200 S288C genomic windows of the same length with AT within ±2 % of the CDEII,
  excluding ±2 kb around CENs (sampled with `ltr_cde.sample_windows`, then AT-filtered).
Per CEN: empirical p for each statistic. Across the 16 CENs (modal alleles): Wilcoxon
signed-rank of observed vs null median. This also tests the "CDEIIs resemble each other more
than random AT-rich DNA" item, using the same statistics instead of pairwise similarity.

### F2.2 Shape along CDEII
Profiles of the four features for the 16 modal alleles and S288C; summary per CEN (mean MGW,
fraction of positions with MGW < 4.0 Å, mean ProT) against the same two nulls.

### F2.3 Do variants preserve shape more than sequence?
- Alleles: Helsen alleles that differ from the modal allele of their CEN by substitutions
  only (same length, Hamming 1–10). Indel alleles are out (no position-wise comparison).
- Observed: shape distance allele vs modal, per feature.
- Null: 200 mutants per allele. Each observed substitution b1→b2 is replaced by b1→b2 at a
  random position of the modal allele that carries b1 (no position reused). The mutation
  spectrum is kept exactly; only positions move. So the test asks whether natural
  substitutions sit at shape-tolerant positions.
- Statistic: percentile of the observed distance in its null. Across alleles: one-sided
  Wilcoxon of percentile − 0.5, overall and per CEN. Alleles are treated as independent;
  shared ancestry is a stated caveat (same as B).
- ScRAP substitution-only loci are shown in the same figure as a consistency check.

### F2.3b Unique substitution events (added during implementation)
The allele-level test counts shared substitutions many times (pseudo-replication; the 309
single-substitution alleles alone give medians 0.44–0.48). Each unique event (CEN, pos,
b1>b2) among substitution-only Helsen alleles is therefore counted once. Its shape change
is ranked exhaustively among all single mutants b1>b2 at the other same-base positions,
and among the same-trinucleotide positions (mutational context), requiring >= 5 null
positions. Outputs `shape_events.tsv`, `shape_events_summary.tsv`, plus
`shape_variants_summary.tsv` for the allele test. This is the test the report interprets.

## Decisions (and rejected alternatives)
- **DNAshapeR (Bioconductor)**, called once per run from R on one FASTA with every sequence
  (real + nulls, about 1e5 × ~85 bp), like `trees_nj.R`. Install
  `BiocManager::install("DNAshapeR")` in the user R library.
  Rejected: Deep DNAshape (TensorFlow dependency, heavy); re-implementing the pentamer
  table in Python (step-feature averaging is easy to get subtly wrong).
- **Helsen alleles are the main data for F2.3**: ScRAP has too few substitution alleles.
- **Position-shuffled, spectrum-matched null** instead of uniformly random substitutions:
  A→T and A→G change shape very differently, so a uniform null would mostly measure
  the mutation spectrum.
- **AT-matched genomic windows + dinuc shuffles** as two nulls: shuffles keep composition
  but break real-DNA structure; genomic windows keep that structure but are rarer at
  CDEII's AT content (~90 %). If < 50 windows match for a CEN, widen to ±4 % and report it.
- Deferred: whether 10-bp vs 5-bp indel alleles (Helsen) keep phasing. This is a natural
  follow-up but needs alignment-aware shape comparison.

## Outputs
- `output/scrap/cdeii_grammar/atracts.tsv` — per sequence: set, CEN, allele, tracts, phasing, power.
- `output/scrap/cdeii_grammar/grammar_null.tsv` — per CEN × statistic × null: obs, null median, p.
- `output/scrap/cdeii_grammar/shape_profiles.tsv` — long: CEN, allele, pos, feature, value.
- `output/scrap/cdeii_grammar/shape_null.tsv` — per CEN shape summaries vs nulls.
- `output/scrap/cdeii_grammar/shape_variants.tsv` — per allele: n_subs, distance per feature, null percentile.
- Code: `src/centromere/grammar.py`, `scripts/14_cdeii_grammar.py`, `scripts/dnashape.R`,
  `tests/test_grammar.py`; config section `grammar:` in `config.scrap.yaml`.
- Report: section `#sec-cdeii-grammar` after `#sec-diversity`; helpers in
  `results/report_helpers.R`; figures `fig_f2_atract_phasing`, `fig_f2_shape_profiles`,
  `fig_f2_shape_variants`.

## Checks
1. Unit tests: tract caller on hand cases (above, plus `AAATAAA`, GC breaks); phasing score
   higher for a synthetic sequence with tracts every 10 bp than for its shuffles.
2. Shape sanity: MGW at the centre of `AAAAAA` is lower than in `GCGCGC` (A-tracts narrow
   the minor groove).
3. Null calibration: alleles simulated from the null itself give uniform percentiles
   (KS p > 0.05).
4. `.venv/bin/pytest` passes, apart from the two known kaleido tests.
5. Report renders.
