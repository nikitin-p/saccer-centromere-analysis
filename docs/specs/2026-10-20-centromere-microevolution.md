# Spec — centromere microevolution in *S. cerevisiae* (2026-10-20)

## Goal
Describe how point centromeres change inside one species — CDEII length steps, the
mutations that make them, the transposons around them, and whether the kinetochore
proteins that read them vary with them — and place that against the literature on
centromere evolution between species (Helsen 2026, Haase 2026, Hession 2025,
Włodzimierz 2026; links in `README.md`).

## Framing (what is already known, what is ours)
- **Known:** CDEII-length variants in *S. cerevisiae* (~9.5 % of 1,493 strains, 13
  variants, mostly short tandem duplications; ~10-bp steps tolerated, 5-bp not; Cbf1
  coevolves with transitions between species) — Helsen. Point centromeres derive from
  Ty5 LTRs — Haase. Centromere *location* conserved ~200 Myr — Hession. The CDE is one
  centromeric particle: one DNA wrap on CDEII (~80 bp) plus Cbf1 and CBF3 — Krassovsky 2012.
- **Ours to add:** T2T, haplotype-aware data; substitution-level constraint in CDEI/CDEIII
  (Helsen studied CDEII length); kinetochore variation *within* the species; transposons
  at point centromeres across strains; diversity at the centromere vs its flanks;
  centromere vs genome trees; aneuploidy as a phenotype.
- **Tested and dropped:** CEN14–CEN15 CDEII exchange. Not significant under a
  run-preserving null (z = 2.5, p = 0.016, Bonferroni n.s.), rank 5–67/120 across strains,
  CEN14–CEN4 equally similar. Recorded in the report as a negative result; the general
  "CDEIIs resemble each other more than random AT-rich DNA" signal is kept as an
  exploratory item (F2).

## Data and dependencies

| Item | Status | Needed for |
|---|---|---|
| 10 ScRAP assemblies, pipeline outputs (`output/scrap/`) | have | all |
| S288C genome (R64, outside the repo) | have | B, D, E |
| Peter 2018 metadata incl. `Aneuploidies` | have (`data/`) | G |
| Loegler 2025 tables S3 (Ty SVs), S5 (SV hotspots) | have (`data/`) | D |
| 748-strain BLAST table (`~/Desktop/yeast_centromeres/sacer_figures/1000_centromeres.tsv`) | have; re-run through fixed pipeline | G |
| **Helsen Supplementary Data 4** (per-strain CEN annotations) | download, small table | A |
| **PCAn v1.0** (Zenodo 10.5281/zenodo.17293587) | download, small | A |
| **Kinetochore protein sequences** (S288C, SGD/UniProt; ~10 genes) | download, a few KB | C |
| **Ty1–Ty5 reference sequences and LTRs** (SGD) | download, a few KB | D |
| minimap2, HMMER | Homebrew (locator spec) | E (flank windows, alignments) |

Downloads are made only after your approval; all are small. Nothing large is fetched.

## Analyses

### A. Validation against the published atlas
1. Match our 10 assemblies to Helsen's Supplementary Data 4 (by strain code); compare
   CDEII length and flagged variants per centromere. Expected: CFS CEN16 (105 bp) flagged;
   ASN CEN3/CEN9 (78 bp) and S288C CEN4 (77 bp) checked against their calls.
2. Run PCAn on the 10 assemblies; compare CDEI/CDEIII coordinates with ours (`cde.py`).
   PCAn becomes a fifth locator in the locator spec's reconciliation if it adds anything.
- Output: `output/scrap/validation_helsen.tsv`, `validation_pcan.tsv`; report table.

### B. CDEII length in helical steps
1. Per locus: ΔCDEII vs S288C and vs the strain-population median for that centromere.
2. Classify: within ±2 bp; ≈ ±10 bp (one turn, 8–12); ≈ ±20 bp; off-step (3–7, 13–17).
   Off-step variants are the interesting ones under Helsen's model.
3. For every variant: is it a tandem duplication/deletion? Report unit length, copy
   orientation, flanking microhomology (exact repeat at the junction), position in CDEII.
   Include Helsen's 13 published variants (from A) to test whether units cluster near
   10 bp multiples.
4. Haplotypes: for phased assemblies, and for collapsed ones by mapping ScRAP reads
   (only if reads are already local — no download), check whether off-step variants
   (ASN CEN3/CEN9) are heterozygous.
- Figure: ΔCDEII per locus on a 10.2-bp grid, coloured by class; duplication diagrams.

### C. Kinetochore genes and coevolution
1. Genes: CBF1 (CDEI); CEP3, NDC10 (CBF2), CTF13, SKP1 (CBF3, CDEIII); CSE4, SCM3;
   MIF2; OKP1, AME1 (COMA).
2. Extraction: `tblastn` of S288C proteins against each assembly's existing BLAST DB;
   take the best hit spanning the full protein; extend to start/stop; translate; flag
   frameshifts/stops. `output/scrap/kinetochore_genes.fasta` + `kinetochore_variants.tsv`.
3. Per gene: protein alignment (MAFFT), haplotypes, variable positions mapped to domains —
   Cbf1 N-terminal disordered tail and bHLH-zip (Helsen); Cse4 N-terminal tail/END domain
   and loop 1 (Włodzimierz's high-entropy regions); Cep3 Zn2Cys6 domain (contacts CDEIII CCG).
4. Association: Cbf1 haplotype vs the strain's CDEI motif variants (e.g. ATCACGTG vs
   ATCACATG across its 16 CENs) and CDEII length variants; Cep3/Ndc10 haplotypes vs CDEIII
   variants. Descriptive at n = 10 (tables, haplotype × centromere-variant matrix);
   permutation tests only once ≥ 50 genomes are available.
- Figure: haplotype matrix (strains × kinetochore genes) beside the centromere variant map.

### D. Transposons around centromeres
1. Annotate full Ty1–Ty5 and solo LTRs (delta, sigma, tau, Ty5 LTR) with BLAST against
   each assembly; keep hits ≥ 80 % identity / ≥ 70 % length.
2. Per CEN: elements within ±1, ±5, ±20 kb; distance of nearest element; orientation.
   Null: same counts for 1,000 random positions per genome matched for distance to the
   nearest tRNA gene (Ty1/Ty3 target Pol III genes, so tRNA proximity is the confounder).
3. Polymorphic insertions: an element near a CEN in some strains but not others
   (presence decided at orthologous positions via the flank windows).
4. Cross-check against Loegler Table S3/S5 positions.
- Output: `te_near_cen.tsv`; figure: TE positions along ±20 kb of each CEN, per strain.

### E. CDE-like motifs in Ty LTRs (Haase, within *S. cerevisiae*)
1. Score every annotated LTR with our CDEI and CDEIII PWMs and the CDEII AT/length
   criteria; compare with composition-matched shuffled LTRs and with random intergenic
   windows. CDEI's CACGTG core is a common E-box, so the control is essential.
2. Any LTR scoring as a full CDEI–CDEII–CDEIII unit is listed as a "dormant
   proto-centromere" candidate; the locator spec's CEN-like scan uses this LTR annotation
   to label its hits instead of reporting them as neocentromere candidates.

### F. Diversity at the centromere vs its flanks
1. π (pairwise differences per site) in CDEI, CDEII, CDEIII and in 100-bp windows across
   ±2 kb flanks, from per-CEN alignments of the flank windows (locator spec L2/L3).
2. Question (Bensasson 2008, Haase 2026): does diversity peak at the centromere, dip, or
   follow the flanks? Per-element and per-centromere.
3. *Exploratory (F2):* CDEII "grammar" — A/T-tract lengths and their phasing (~10 bp),
   and predicted DNA shape (minor-groove width, roll) along CDEII; do variants preserve
   shape more than sequence? Only if F1 is done.

### G. Aneuploidy as a centromere phenotype (748-strain set)
1. Re-run the 748-strain BLAST table through the fixed pipeline (assign_loci, merge fix;
   CDE annotation only if those assemblies are local — otherwise BLAST loci only).
2. Join Peter 2018 `Aneuploidies` (e.g. `aneu;+1*9;`) per strain and chromosome.
3. Test: are strains aneuploid for chromosome *k* enriched for non-S288C alleles or
   length variants at CEN*k*? Logistic model per chromosome with clade as covariate;
   chromosomes with ≥ 10 aneuploid strains only (1, 3, 8, 9, 11, 14).
- Separate output dir `output/peter748/`; this does not change the ScRAP report.

### H. Centromere tree vs genome tree
1. Centromere tree: concatenated 16 CDE sequences per assembly (~1.9 kb), ML tree
   (IQ-TREE if installed, otherwise neighbour-joining in R/ape on K80 distances).
2. Genome tree: SNP distances from minimap2 whole-genome alignments to S288C (locator
   spec phase 2) for the same assemblies.
3. Compare: tanglegram (`ape::cophyloplot` / `phytools::cophylo`), Robinson–Foulds,
   Mantel test of the two distance matrices. Strains placed differently are listed with
   the centromeres that drive the difference (per-CEN distance contribution).
4. Per-CEN trees are not built (too few sites); across-species CEN trees wait for other
   *Saccharomyces* genomes.

## Files

| File | Change |
|---|---|
| `src/centromere/microevo.py` (new) | B (steps, duplication detection), F (π), E (LTR scoring) |
| `src/centromere/genes.py` (new) | C: tblastn extraction, translation, domain mapping |
| `scripts/06_microevolution.py` (new) | runs A–F on `config.scrap.yaml`; writes TSVs |
| `scripts/peter748.py` (new) | G, own config `config.peter748.yaml` |
| `data/kinetochore_s288c.faa`, `data/ty_reference.fasta` | downloaded (small), committed |
| `results/report.qmd`, `report_helpers.R` | sections: Validation, CDEII steps, Kinetochore, Transposons, Diversity, Trees, CEN14–CEN15 negative result |
| `tests/test_microevo.py`, `tests/test_genes.py` | below |

## Decisions (and the rejected alternative)
- **Within-species focus, literature for the between-species picture.** *Rejected:*
  re-deriving Helsen's atlas or Haase's origin analysis ourselves.
- **Validate against Helsen before claiming any CDEII result.** Their 1,493 strains are the
  reference; our contribution is T2T resolution and the analyses they did not do.
- **Kinetochore association is descriptive until ≥ 50 genomes.** *Rejected:* p-values at
  n = 10.
- **tRNA-matched null for transposons.** *Rejected:* uniform random positions (would
  confuse Ty1/Ty3 tRNA targeting with centromere attraction).
- **Aneuploidy on the 748-strain set, not ScRAP.** Only it has enough aneuploid strains.
- **CEN14–CEN15 dropped** after the run-preserving null; kept as a documented negative.

## Checks that prove it works
1. A: ≥ 90 % of our CDEII lengths equal Helsen's for the same strain/centromere (±1 bp);
   PCAn CDEI/CDEIII starts within 1 bp of ours for ≥ 155/160 loci; every disagreement listed.
2. B: CFS CEN16 detected as a 22-bp tandem duplication with junction microhomology;
   synthetic test — a 10-bp tandem duplication inserted into a CDEII is classified "one turn".
3. C: CBF1, CEP3, CSE4 extracted full-length in 10/10 assemblies; S288C proteins recovered
   100 % identical when the pipeline is run on the S288C genome itself.
4. D: on S288C, every SGD-annotated Ty element within ±20 kb of a CEN is recovered.
5. E: shuffled-LTR control gives ≤ 5 % "CDE units"; real-LTR rate reported against it.
6. F: π computed on a synthetic alignment with known differences returns the exact value.
7. G: aneuploidy parser reproduces the counts in the Peter 2018 table (793 euploid; 31 `+1*1`).
8. H: tanglegram and RF run on the 10 assemblies; a simulated swap of one strain's CEN
   block is detected as that strain's discordance.
9. `pytest -q tests` passes; report renders; new sections checked by screenshot.

## Order
A → B → C (no dependency on the locator spec) · D, E → F, H (need flank windows /
minimap2 from the locator spec) · G (independent; can run any time).

No git writes; you stage and commit.
