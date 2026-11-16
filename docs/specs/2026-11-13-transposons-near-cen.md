# Spec — analysis D: transposons around centromeres (2026-11-13)

Implements section D of `2026-10-20-centromere-microevolution.md` on the 10 ScRAP
assemblies + S288C. Goal: annotate Ty1–Ty5 and solo LTRs genome-wide, describe what
sits within ±20 kb of each CEN, test enrichment/depletion against a tRNA-matched null,
and find insertions that are polymorphic between strains.

## Inputs
- `data/external/sgd/saccharomyces_cerevisiae.gff.gz` (SGD, downloaded 2026-11-13,
  gitignored): `LTR_retrotransposon` (50), `long_terminal_repeat` (384), `tRNA_gene`
  (299), `centromere`, and the R64 genome in its `##FASTA` block.
- Existing BLAST DBs `output/scrap/kinetochore_blastdb/` (10 assemblies + sacCer3).
- CEN loci: `output/scrap/filtered.tsv` (`sseqid`, `cen_start`, `cen_end`, `strand`).
- Loegler 2025 Table S5 (55 SV hotspots, S288C coordinates). S3 is only counts → not used.

## Steps and outputs (all under `output/scrap/` unless committed in `data/`)
1. `scripts/make_ty_reference.py` → `data/ty_reference.fasta` (committed, ~150 KB):
   one intact representative per family Ty1–Ty5 (longest SGD element) + **all** 384
   S288C LTRs (headers `LTR|<family>|<SGD name>`; family from the SGD Note:
   delta=Ty1/Ty2, sigma=Ty3, tau=Ty4, Ty5 LTR). Also `data/s288c_trna.fasta`
   (275 nuclear tRNAs) and `data/s288c_te_sgd.tsv` (SGD coordinates, for check 1).
   **Asserts** the 16 SGD CEN intervals cut from the GFF genome equal
   `data/sacer_centromere.fasta` and that the GFF genome equals sacCer3 on chrI–XVI.
2. `blastn -task blastn` of the library vs each genome. Hit kept if pident ≥ 80 % and
   aligned ≥ 70 % of the query (spec rule). Overlapping hits merged into elements.
   Element class: `full` if ≥ 70 % of a full-element *internal* region is covered,
   else `solo_LTR`; family from the best-scoring query. → `te_elements.tsv`
   (genome, contig, start, end, strand, class, family).
3. tRNAs per genome: blastn of `s288c_trna.fasta`, ≥ 90 % id, ≥ 90 % length → `trna_loci.tsv`.
4. Per CEN: elements within ±1, ±5, ±20 kb of the CEN edges, nearest element distance,
   orientation relative to the CEN strand → `te_near_cen.tsv` (one row per element in
   the window) + `te_cen_summary.tsv` (genome × CEN counts, nearest distance).
5. **Null:** per genome 1,000 random positions, ~62 matched to each CEN: drawn among
   positions on the same assembly whose nearest-tRNA distance is within ±10 % of that
   CEN's, ≥ 20 kb from contig ends and from every CEN. Count elements in the same
   windows; per window size: observed (16 CENs) vs null mean, empirical two-sided p,
   per genome and pooled. A uniform-position null is reported alongside, labelled
   secondary, to show how much tRNA proximity matters → `te_null.tsv`.
6. **Polymorphic insertions (orthology by junction flanks):** for each element in a
   ±20 kb CEN window, the 300 bp left and right of it are BLASTed (megablast, ≥ 90 % id,
   ≥ 80 % length, best hit unique: runner-up bitscore < 0.9 × best) against every
   other genome. Both flanks on one contig, consistent orientation: inner gap ≤ 50 bp
   → `absent`; gap overlaps an annotated element → `present`; otherwise
   `rearranged`; a flank not found uniquely → `unresolved`. Elements in different
   genomes that are each other's `present` hit are one site (union-find)
   → `te_sites.tsv` (site × strain matrix, with class/family per strain).
7. Loegler S5 cross-check: hotspots overlapping ±20 kb of a CEN in S288C coordinates,
   with their TE/LTR/tRNA flags, next to our S288C calls → `te_loegler_s5.tsv`.

## Files
| File | Change |
|---|---|
| `scripts/make_ty_reference.py` | new (step 1) |
| `src/centromere/transposons.py` | new: hit filtering, merging, classification, window counts, matched null, flank-orthology calls, site clustering |
| `scripts/10_transposons.py` | new: steps 2–7; `PYTHONPATH=src .venv/bin/python scripts/10_transposons.py --config config.scrap.yaml` |
| `config.scrap.yaml` | new section `transposons:` |
| `tests/test_transposons.py` | new |
| `results/report.qmd`, `results/report_helpers.R` | section `#sec-transposons` (after G, before length): summary table, null table, figure of TE positions along ±20 kb per CEN × strain (colour = family, shape = full/solo), polymorphic-site table |

## Decisions (rejected alternative)
- **All 384 S288C LTRs as queries**, not one consensus per family: solo deltas are
  diverged (~70–90 %), and the best of many copies passes 80 % where a consensus would not.
- **Element = merged hits, class decided by internal coverage.** *Rejected:* counting
  BLAST hits (a full Ty would count as 3: two LTRs + body).
- **Junction-flank BLAST for orthology** (user choice 2026-11-13). *Rejected for now:*
  minimap2 synteny (locator phase 2, new dependency).
- **tRNA-matched null as primary, uniform as secondary.** Both in one table.
- **Truncated fragments (< 70 % of any query) are dropped**, as the spec rule says;
  their count is logged so the loss is visible.

## Checks
1. S288C: every SGD Ty and LTR within ±20 kb of a CEN recovered (overlap ≥ 50 %);
   genome-wide recall of SGD LTRs and Ty elements printed.
2. Synthetic: a Ty1 copy inserted into a random sequence is called one `full` Ty1; a
   solo delta → one `solo_LTR`; flank calls `present`/`absent` on a genome with and
   without the insertion.
3. Null sampler: drawn positions' nearest-tRNA distances within ±10 % of their CEN's.
4. `.venv/bin/pytest -q tests` passes; report renders; new section checked by screenshot.

No git writes; you stage and commit.
