# Spec — analysis E: CDE-like motifs in Ty LTRs (2026-11-13)

Implements section E of `2026-10-20-centromere-microevolution.md` on the 10 ScRAP
assemblies + S288C. Goal: do Ty LTRs carry CDEI / CDEIII-like motifs, or full
CDEI–CDEII–CDEIII units, more often than composition-matched or genomic background
predicts? Any LTR with a full unit is listed as a "dormant proto-centromere" candidate
for the locator spec's CEN-like scan.

## Inputs
- `output/scrap/te_elements.tsv` (analysis D): 2,120 elements, 1,824 solo LTRs
  (median 339 bp), 193 full elements.
- `data/sacer_centromere.fasta` → PWMs via `cde.CDEModel` (CDEI 8 bp, CDEIII 26 bp).
- Assemblies (same paths `10_transposons.py` uses) + sacCer3.
- SGD GFF `data/external/sgd/saccharomyces_cerevisiae.gff.gz` (genes, for S288C intergenic).

## Definitions (from the 16 S288C CENs)
- **CDEI hit**: PWM score ≥ `CDEModel.cdei_min` (weakest S288C − 2 = 8.09; max 13.46).
- **CDEIII hit**: score ≥ `cdeiii_min` (24.17; max 38.22).
- **CDEII**: length 77–86 bp in S288C, AT 0.86–0.98 → criterion length 70–95 bp and
  AT ≥ 0.80 (`ltr_cde:` config; loose by design, so candidates are not missed).
- **Unit**: CDEI hit, then CDEII passing the criterion, then CDEIII hit, all on one
  strand. Both strands scanned. Unit score = CDEI + CDEIII score.
- Also kept per sequence: best CDEI score, best CDEIII score, number of hits, and
  E-box count (CACGTG / CACATG) so CDEI hits can be read as "mostly E-boxes" or not.

## Steps and outputs (`output/scrap/` unless noted)
1. **LTR set.** Solo LTRs = element interval. Full elements → two LTRs = terminal
   segments of the family's median S288C LTR length (Ty1/2 delta, Ty3 sigma, Ty4 tau,
   Ty5), refined to the extent of the best LTR HSP at each end if it exists.
   Sequences extracted in element orientation. → `ltr_seqs.fasta`.
   Identical sequences across strains collapsed to a `seq_id` (strain list kept),
   so pooled counts are not driven by the same LTR present in 11 genomes.
2. **Scan** every LTR (both strands) → `ltr_cde_scores.tsv` (one row per LTR: genome,
   element_id, family, class, seq_id, length, AT, best CDEI/CDEIII, n hits, n units,
   best unit coords and sequences).
3. **Controls**, each run on unique LTR sequences (`seq_id`):
   a. *Shuffled*: 100 dinucleotide-preserving shuffles per LTR (Altschul–Erickson);
      keeps AT-tract structure better than a mononucleotide shuffle.
   b. *Genomic*: per genome, windows with the LTR length distribution drawn at random
      outside TEs, tRNAs and CENs ±1 kb (1,000 resamples of the full LTR set).
   c. *S288C intergenic*: same as (b) but restricted to SGD intergenic regions
      (gene-free), as the spec asks; only possible for S288C.
   Statistic per family and pooled: hits per kb for CDEI, CDEIII, units; observed vs
   null mean, empirical two-sided p (BH across tests). → `ltr_cde_control.tsv`.
4. **Positive control**: the scanner finds a unit in each of the 16 S288C CEN
   sequences and in 160/160 ScRAP CEN loci from `filtered.tsv` (within ±2 bp of the
   `cde.py` calls). If it fails, the definitions above are wrong and E stops.
5. **Candidates**: every LTR with ≥ 1 unit → `ltr_cde_candidates.tsv` (seq_id, strains,
   family, unit position in LTR, the three element sequences, scores, nearest CEN
   distance from `te_near_cen.tsv`). This is the label file for the locator's CEN-like scan.

## Files
| File | Change |
|---|---|
| `src/centromere/ltr_cde.py` | new: unit scan (reuses `cde.build_pwm`/`score_positions`), dinucleotide shuffle, LTR extraction, random-window sampler, rate comparison |
| `scripts/11_ltr_cde.py` | new: steps 1–5; `PYTHONPATH=src .venv/bin/python scripts/11_ltr_cde.py --config config.scrap.yaml` |
| `config.scrap.yaml` | new section `ltr_cde:` (CDEII bounds, n_shuffles, n_resamples, seed) |
| `tests/test_ltr_cde.py` | new |
| `results/report.qmd`, `results/report_helpers.R` | section `#sec-ltr-cde` after D: control table, figure (best CDEI / CDEIII score distributions: LTR vs shuffled vs genomic, S288C CEN scores as ticks), candidate table |

## Decisions (rejected alternative)
- Thresholds = weakest S288C element − 2, as in `cde.py` (rejected: p-value cut-offs
  from a background model — would need a separate calibration and diverge from the
  annotation everyone else uses).
- LTRs of full elements included (rejected: solo LTRs only — throws away ~380 LTRs).
- Dedup by identical sequence (rejected: reporting per-genome only — the pooled
  question needs unique sequences; per-genome numbers are still in the scores table).
- Dinucleotide shuffle (rejected: mononucleotide — breaks A/T runs and makes CDEII
  look artificially rare in the null).
- Genomic control on all genomes, intergenic only for S288C (rejected: predicting genes
  in 10 assemblies — out of scope; a gene-free control is what the spec asked for, so
  S288C carries it).
- No HMM / DNA-shape scoring of CDEII (belongs to F2 / locator).

## Checks
- `.venv/bin/pytest tests/test_ltr_cde.py`: unit found in a synthetic CDEI+CDEII+CDEIII
  insert; not found when CDEII is 60 bp or AT 0.6; reverse-strand unit found; shuffle
  preserves dinucleotide counts exactly; LTR extraction coordinates on a toy element.
- Step 4 positive control: 16/16 and 160/160 printed by the script.
- Whole suite `.venv/bin/pytest` (2 kaleido tests may fail as before).
- Report renders.
