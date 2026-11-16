# Spec — independent centromere locators and an expanding HMM (2026-10-20)

## Goal
Find every centromere by several methods that fail in different ways — by its own
sequence, by its flanking DNA, by whole-genome synteny, and by a population profile
that grows as strains are added — and reconcile them, so a centromere too diverged
for BLAST, moved by a rearrangement, or carrying an insertion is caught and reported
instead of silently missing.

## Where we are
- **L1 — BLAST query (exists):** the 16 S288C CEN elements BLASTed against each
  assembly; cross-centromere hits removed; split hits merged; CDEI/CDEIII placed with
  S288C PWMs (`src/centromere/cde.py`). 160/160 in the 10 ScRAP genomes.
- **Blind spot:** everything depends on similarity to S288C's centromere sequence.
  A centromere BLAST cannot see is simply absent from the table.

## The four locators

| | Locator | Depends on | Fails when | Dependency |
|---|---|---|---|---|
| L1 | BLAST of S288C CENs | CEN sequence ≈ S288C | CEN diverged, large insertion | BLAST (have) |
| L2 | BLAST of 2 kb S288C flanks | flanks conserved and adjacent | SV breakpoint near the CEN | BLAST (have) |
| L3 | minimap2 whole-genome synteny | genome collinear with S288C around the CEN | flanks deleted/duplicated; distant species | minimap2 |
| L4 | Iterative per-CEN HMM | the population profile | the profile drifts (guarded, below) | HMMER |

L2 and L3 are sequence-agnostic about the centromere itself; L4 removes the S288C
bias of L1 by learning from the strains. Phases are independent: each lands with its
own tests and report rows.

## Data flow
```
S288C genome (R64, NC_001133–NC_001148, outside the repo)
  ├ make_cen_flanks.py ─► data/s288c_cen_flanks.fasta, data/s288c_cen.bed
assemblies
  ├ L1 blast.sh                  ─► centromeres.tsv          (unchanged)
  ├ L2 blast.sh MODE=flanks      ─► cen_flanks.tsv
  ├ L3 minimap2 -x asm5 -c --cs  ─► <sample>.paf ─► synteny_loci.tsv
  └ L4 hmm_iterate.py            ─► hmm/round_<k>/CEN*.hmm, hmm_loci.tsv, hmm_rounds.tsv
01_preprocess: reconcile(L1..L4) → main table (+ `locator`, `support`) → annotate_cde
               cen_like scan ─► cen_like_elements.tsv
report.qmd: section "Four locators", HMM growth figure, disagreement table
```

## Phase 1 — L2 flanks (no new dependency)
- `scripts/make_cen_flanks.py`: reads the S288C genome (path argument) and the SGD
  coordinates in the `data/sacer_centromere.fasta` headers; **asserts** the 16
  extracted CENs equal the FASTA (guards against a reference-version mismatch);
  writes the 64 kb flank FASTA and a CEN BED (both committed).
- `blast.sh MODE=flanks`: `blastn -task megablast` against the existing DBs.
- Locus: best left + best right flank (pident ≥ 90 %, ≥ 50 % of flank aligned) on the
  same contig in consistent orientation; locus = between their inner ends.

## Phase 2 — L3 minimap2 synteny
- `minimap2 -x asm5 -c --cs` each assembly to the S288C genome (`asm10`/`asm20`
  selectable for divergent genomes). PAF kept per sample in the run's output dir.
- **Own liftover from the `cs` tag** (`locate.liftover_paf`): for each CEN BED
  interval, find the alignment blocks covering its start and end and walk the `cs`
  string to the exact query coordinate. *Rejected:* `paftools.js liftover` (needs the
  k8 JavaScript runtime, a second dependency).
- A CEN whose ends lift to different contigs or only one end lifts → reported with
  the reason (this is where rearrangements show up).

## Phase 3 — L4 iterative, expanding HMM
**Round 0.** One HMM per centromere, built (`hmmbuild --dna`) from that centromere's
alignment in `aligned.tsv` — the current strains plus S288C. One profile per CEN, not
one for all 16 (the coursework's single profile was weak in CDEII and clipped the ends).

**Round k → k+1.**
1. `nhmmer` all 16 profiles against every assembly (both strands, E ≤ 1e-10).
2. **Accept** a new sequence for CEN *i* only if all hold:
   - its best-scoring profile is CEN *i*, by ≥ 10 bits over the runner-up
     (reciprocal check — this is what stops CEN14/CEN15 cross-capture);
   - CDEI and CDEIII pass the PWM `cde_ok` threshold inside the hit;
   - length within S288C ± 40 bp, or flagged as `insertion` if L2/L3 agree on the locus;
   - it is not already in the profile (identical sequences collapsed).
3. Realign the accepted sequences per CEN (MAFFT), rebuild the HMM. HMMER's default
   sequence weighting keeps a clade with many near-identical strains from dominating.
4. **Stop** when a round accepts nothing new, or after `max_rounds` (default 5).

Every round writes `hmm_rounds.tsv` (round, CEN, sequences in profile, accepted,
rejected by reason, mean bits of the S288C sequence against its own profile). A falling
S288C score with stable acceptance is the expected sign of the profile widening; a
rising rejection count for "runner-up within 10 bits" is the sign of drift, and stops
that CEN's iteration.

**Designed to grow with the data:** the same loop seeds from the 10 ScRAP genomes now,
absorbs the 142 when they arrive, and later other *Saccharomyces* species — each
batch widens the profiles that search the next.

## Phase 4 — reconcile and scan
**Per assembly × CEN**, loci from L1–L4 are merged when they overlap on one contig;
`support` = set of locators that found it. Categories:

| Category | Rule |
|---|---|
| `agree` | all locators that ran found the same locus |
| `rescued` | L1 missed, ≥ 1 other found it → enters the main table, CDEs placed in that window |
| `sequence_only` | L1/L4 found it, L2/L3 did not anchor (rearrangement near the CEN) |
| `conflict` | locators give different contigs or non-overlapping loci → listed, not resolved |
| `insertion` | L2/L3 locus longer than S288C by > 30 bp (e.g. a Ty element) |

**CEN-like scan:** the final HMMs plus CDEI/CDEIII PWMs over every contig; hits passing
the reciprocal and CDE checks that do not overlap a known locus → `cen_like_elements.tsv`
(candidates for neocentromeres, dicentric arms, CEN fragments), not calls.

## Files

| File | Change |
|---|---|
| `scripts/make_cen_flanks.py`, `data/s288c_cen_flanks.fasta`, `data/s288c_cen.bed` | new (Phase 1) |
| `scripts/blast.sh` | `MODE=flanks` |
| `scripts/run_minimap2.sh` | new (Phase 2), same resume/skip behaviour as `blast.sh` |
| `scripts/hmm_iterate.py` | new (Phase 3) |
| `src/centromere/locate.py` | new: `locate_by_flanks`, `liftover_paf`, `hmm_loci`, `reconcile`, `scan_cen_like` |
| `src/centromere/cde.py` | genome-wide PWM scoring helper |
| `scripts/01_preprocess.py` | runs reconcile when `locate:` is configured; adds `locator`, `support` |
| `config.scrap.yaml` | `locate: {flanks, paf_dir, hmm_dir, flank_bp: 2000, max_rounds: 5, rescue: true}` |
| `results/report.qmd`, `report_helpers.R` | "Four locators" section; HMM growth figure; workflow figure updated |
| `tests/test_locate.py` | below |

## Decisions (and the rejected alternative)
- **All four locators, not one replacement.** They fail in different ways; agreement
  is the evidence. *Rejected:* replacing BLAST with the HMM alone.
- **Reciprocal best-profile + CDE gate for HMM growth.** *Rejected:* E-value alone
  (would re-create cross-centromere capture as the profiles widen).
- **S288C genome stays outside the repo** (12 MB); only derived small files committed.
- **Rescue on by default for ScRAP**, every row tagged with `locator` so a rescued
  locus is always traceable. *Rejected:* report-only (keeps the blind spot in every figure).
- **Dependencies:** minimap2 and HMMER via Homebrew (a few MB each), installed only
  after your OK; Phase 1 needs neither.

## Checks that prove it works
1. Flank script assertion: 16/16 extracted CENs equal the FASTA.
2. Real data: 160/160 `agree` across all four locators; CFS CEN16 = 139 bp from L2, L3, L4.
3. **Rescue:** copy of one assembly with CEN5's CDEII replaced by random 90 %-AT
   sequence → L1 misses, L2/L3/L4 find it → `rescued`, CDEs correct, row in the table.
4. **Harder rescue:** CDEI/CDEIII also mutated (2–3 bp each) → shows which locators
   still find it; result tabulated, not required to pass.
5. **Insertion:** 6 kb inserted in a CDEII → `insertion`, gap ≈ 6 kb.
6. **Rearrangements:** CENs on `chrXV_chrIX` (CNT), `chrVI_chrXIII`, `chrXIII_chrIX`
   (BCE_3a) → `agree`.
7. **HMM generalisation:** leave-one-clade-out — build round 0 without one clade,
   iterate on the rest, then search the held-out genomes: recall 16/16 per genome.
8. **HMM drift guard:** on the 10 genomes no accepted sequence is a cross-centromere
   hit (every accepted locus agrees with L2/L3).
9. CEN-like scan on S288C: no hit outside the 16 loci, or each one listed and explained.
10. `pytest -q tests` passes; report renders; new section checked by screenshot.

## Out of scope
Polarising variants with *S. paradoxus*; phenotype association; cluster runs (same
code path, run when the 142 genomes are available).

No git writes; you stage and commit.
