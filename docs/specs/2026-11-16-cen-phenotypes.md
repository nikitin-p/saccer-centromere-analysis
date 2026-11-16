# Spec — centromere alleles vs phenotypes (Loegler 2025 traits) (2026-11-16)

## Goal
Test whether natural centromere variation (rare/non-modal alleles, CDEII length) is
associated with fitness under segregation-relevant stress or with expression of nearby and
kinetochore genes, in ~850 strains, with a mixed model that corrects for relatedness.

## Data

| Input | Source | Strains |
|---|---|---|
| CEN alleles | Helsen 2026 SD4, Peter-genome rows (3-letter codes), variants A (+ B if heterozygous) | 1075 |
| Genotypes for kinship | Peter 2018 `1011GWASMatrix` (plink; 83,794 markers, MAF > 5 %; chr 17 = CNVs) in `data/external/peter2018/` (gitignored; md5 checked) | 1011 |
| Phenotypes | `Phenotypes_8391Traits/` (Growth 241, Expression 6003, Protein 2147; `ID_STRAIN` + value) | ~900 |
| Covariates | Peter 2018 Table S1: ploidy, aneuploidies (parsed by `aneuploidy.parse_aneuploidies`) | 1011 |

Intersection with alleles + kinship: 846 (growth), 814 (expression). **As built:** only strains
with all 16 CENs in SD4 are used (801 of 917), so one sample set serves every CEN of a trait
(one eigendecomposition per kinship); in practice n = 750 (growth), 726 (expression).

## Predictors (per strain)

| Name | Definition | Per |
|---|---|---|
| `non_modal` | share of the strain's variants (A, B) that are not the per-CEN population major allele: 0, 0.5, 1 | CEN |
| `cdeii_delta` | mean over variants of CDEII length − per-CEN median (CDEII = seq[8 : −27], as `microevo.helsen_cdeii`) | CEN |
| `off_step` | any variant with \|Δ\| ≥ 3 and not within ±2 of 10 or 20 bp (not a whole helical turn) | CEN |
| `rare_burden` | number of CENs carrying an allele with population frequency < 5 % | strain |

A per-CEN predictor is tested only if ≥ 20 strains are in the minority class (binary) or
≥ 20 strains have Δ ≠ 0 (`cdeii_delta`) in that trait's sample.

## Traits

**Primary (hypothesis-driven), reported with BH q within the primary family:**
1. Growth under segregation / genotoxic stress: all benomyl, 4-NQO, hydroxyurea,
   camptothecin traits (SC* and YPD* sets, all time points) — 32 traits × 16 CEN ×
   3 predictors + burden.
2. *Cis* expression and protein: genes whose ORF lies within ±10 kb of the CEN (SGD GFF,
   `gene` features with systematic IDs) vs that CEN's predictors only.
3. Kinetochore genes (the 10 in `data/kinetochore_s288c.faa` + NDC80, CSE4 already included):
   expression/protein vs `rare_burden` and vs each CEN's `non_modal`.

**Secondary:** all 241 growth traits × all predictors, BH within the secondary family;
reported as an exploratory table and a QQ plot, never as the headline.

Every trait is rank-inverse-normal transformed in its own sample; traits with n < 200 dropped.

## Model

EMMAX-style LMM per trait: `y = Xβ + g + e`, `g ~ N(0, σg² K)`.
- **K**: standardised GRM (VanRaden) from the plink SNPs, CNV markers (chr 17) excluded,
  missing genotypes mean-imputed. **LOCO**: for a CEN on chromosome c, K is built without
  chromosome c, so markers linked to the CEN do not absorb its effect. `rare_burden` uses the
  full K.
- Variance components fitted once per trait × K under the null (REML, Brent on
  δ = σe²/σg² over the eigen-decomposition), then each predictor tested by GLS Wald test
  (rotated OLS) — the standard EMMAX approximation.
- Covariates: intercept, ploidy, copy-number change of the CEN's chromosome (gain + / loss −,
  from Peter 2018); for `rare_burden`, total number of aneuploid chromosomes.
- Report β, SE, p, n, h² of the null model, and the genomic inflation λ per family.

Own implementation (numpy/scipy, ~150 lines) in `src/centromere/pheno.py`, own plink `.bed`
reader (2-bit SNP-major). *Rejected:* GEMMA / limix / fastlmm (new heavy dependencies, no
conda pkg cache for them offline); plain OLS with clade covariates (user chose LMM).

## Added during implementation
- **CNV markers** (chr 17): the `.bim` allele labels are the copy states (1, 2) and state 0 is
  plink's missing allele → `cnv_states` decodes to 0/1/2 (dosage of A1 gave the wrong sign).
- **Local-haplotype conditioning:** the 30 smallest-p tests of each family set are refitted with
  5 PCs of the SNPs within ±50 kb of the CEN as extra covariates (`p_local`).
- NDC10 is `CBF2` in the SGD GFF `gene=` field.

## Checks that prove it works
1. **Plink reader:** genotype counts of 5 markers equal those from a direct byte decode
   in a test fixture written by our own writer; fam/bim order preserved.
2. **Null calibration (simulation, test):** 600 simulated strains in 6 clades, trait =
   kinship-structured noise + clade shift, predictor correlated with clade but with no effect,
   200 replicates → LMM false-positive rate at 0.05 within 0.03–0.07; OLS clearly inflated.
3. **Power (simulation, test):** same setup + true effect → LMM recovers β within 2 SE.
4. **Real positive control:** CUP1 copy-number marker (chr 17 CNV covering YHR053C/YHR055C)
   vs copper growth (`YPDCUSO410MM_*`, `SCCuSO4*`) → p < 1e-10 with the same LMM.
   *Cis* control: an ORF's expression vs its own CNV marker, for ≥ 1 CNV gene, p < 1e-10.
5. λ of the primary family between 0.9 and 1.1 (reported either way).
6. `pytest -q tests` passes; report renders; new section checked by screenshot.

## Outputs (`output/peter748/phenotypes/`)

| File | Contents |
|---|---|
| `cen_predictors.tsv` | strain × CEN: modal/non-modal, Δ CDEII, off-step, het flag; strain burden |
| `trait_samples.tsv` | trait, class (growth/expression/protein), family, n, null h² |
| `assoc_primary.tsv` | trait, CEN, predictor, n, n_minor, β, SE, p, q, λ-family |
| `assoc_secondary.tsv` | same for the growth scan |
| `positive_controls.tsv` | CUP1 / cis-CNV results |
| `kinship_loco/` | cached eigendecompositions (npz), one per chromosome + full |

## Files

| File | Change |
|---|---|
| `src/centromere/pheno.py` | new: `read_plink`, `grm`, `fit_null`, `lmm_test`, `inverse_normal`, `cen_predictors`, `cis_genes` |
| `scripts/15_phenotypes.py` | new; `PYTHONPATH=src .venv/bin/python scripts/15_phenotypes.py --config config.peter748.yaml` |
| `config.peter748.yaml` | new section `phenotypes:` (paths, window 10 kb, min_minor 20, min_n 200, trait regexes) |
| `tests/test_pheno.py` | checks 1–3 + predictor definitions on toy alleles |
| `results/report.qmd`, `report_helpers.R` | section `#sec-phenotypes` after G: QQ plot (primary/secondary), forest plot of top hits, positive-control row |

## Decisions (and the rejected alternative)
- **Alleles from Helsen SD4, not our 748 BLAST set.** 846 vs 633 strains with phenotypes,
  full CEN sequence with CDE boundaries. *Rejected:* 748 set (no CDE, fewer strains). Helsen
  = SGD + 1 bp at 3′ was validated in analysis A.
- **Population major allele as reference, not S288C** (lesson from B and G).
- **LOCO kinship.** *Rejected:* one full K (would partly absorb a CEN effect via linked SNPs).
- **Two families, primary reported first.** *Rejected:* one pooled BH over all 8391 traits
  (the hypothesis tests would drown).
- **Aneuploidy and ploidy as covariates**, since aneuploidy changes growth and expression and
  could itself be caused by CEN alleles (G found no association, so adjustment is safe).
- **Expression/protein limited to cis + kinetochore genes.** *Rejected:* full transcriptome scan.

## Not in scope
Fine-mapping within the CEN (single-site tests); heterozygote vs homozygote effects beyond
the 0/0.5/1 dosage; the 142 ScRAP genomes.

No git writes; you stage and commit.
