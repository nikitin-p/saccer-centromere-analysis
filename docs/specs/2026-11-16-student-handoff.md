# Student handoff: genome download + one-command run (2026-11-16)

## Goal
A student clones the repo on a Linux cluster, creates one conda env, downloads ScRAP 142
or Loegler 1086 assemblies with one script, and runs every assembly-based analysis plus
the Quarto report with one command. Nothing points at this Mac.

## Outputs
| Path | What |
|---|---|
| `environment.yml` | conda-forge + bioconda env `centromere` (python, pandas..., blast, mafft, minimap2, hmmer, r-base + tidyverse/ape/vegan/ggseqlogo/ggdist/patchwork/plotly/readxl, bioconductor-dnashaper, quarto) |
| `scripts/download_data.py` | `--set scrap142 \| loegler1086 \| reference` (see below) |
| `config.scrap142.yaml`, `config.loegler1086.yaml` | same sections as `config.scrap.yaml`; paths per set; tools from PATH |
| `scripts/run_all.sh CONFIG [--from STEP] [--only STEP]` | every step in order, resumable (`<output>/.done/<step>`), log per step |
| `scripts/slurm_run_all.sh` | example `sbatch` wrapper (cpus/mem/time as variables) |
| `README.md` section "Quick start (cluster)" | env -> download -> run -> report, expected sizes/times |
| `.gitignore` | un-ignore what a clone needs (below) |

### download_data.py
- `reference`: S288C R64 `sacCer3.fa` (UCSC goldenPath, chr names = chrI..chrM, as the code
  expects) -> `data/external/s288c/sacCer3.fa`; SGD GFF -> `data/external/sgd/`; Peter 2018
  and Loegler 2025 supplementary tables (springer static-content URLs) -> `data/`.
- `scrap142`: Zenodo 8210254 `nuclear.tar.gz` (1.7 GB). One nuclear assembly per strain:
  `collapsed` if present, else `HP0`; phased HP1/HP2 skipped unless `--haplotypes all`.
- `loegler1086`: Zenodo 15698884 `GenomeAssemblies.tar.gz` (16.8 GB). Only the final
  assembly per isolate (`*.Final.fasta`); SecondBest/annotations skipped. Phased isolates:
  HP1 only by default, `--haplotypes all` keeps HP2 as a separate sample.
- Archive is streamed (`tarfile` mode `r|gz`), never stored; `--archive PATH` uses a local
  copy instead. Each kept member is gzipped to `data/assemblies_<set>/<sample>.fasta.gz`
  (`.tmp` + rename, existing files skipped). `--strains A,B` / `--limit N` for test subsets.
- Writes `data/samples_<set>.tsv` (`accession assembly strain`, accession = file stem),
  the format `input.sample_map` already reads. `strain` = Standardized name (Peter 2018 key).
- md5 of the tarball is checked against Zenodo's when the archive is stored locally.

### run_all.sh step order
blast_cen -> pipeline 01-05 -> blast_flanks -> l2_flanks -> minimap2 -> 12_flanks_diversity (L3, F)
-> hmm_iterate (L4) -> pipeline 01-05 again (reconcile) -> cen_like_scan -> locator_synthetic
-> 06 (B) -> 07, 08, 09 (C) -> 10 (D) -> 11 (E) -> 13 (H) -> 14 (F2) -> report.
Off by default (`--with`): validate_external (A, needs Helsen/PCAn), peter748 (G), 15 (phenotypes).

## Code changes this needs
- Tool paths: `kinetochore.blast_bin`, `locate.minimap2`, `locate.hmmer_bin`, the hard-coded
  `BLAST_BIN` in `locator_synthetic.py`: `null` or absent -> found on PATH. One helper in
  `centromere/config.py`; existing configs keep working.
- `run_minimap2.sh` and the scripts that read `cde.genome_dir` expect `<accession>.fasta.gz`;
  download_data.py writes exactly that, so no change there.
- Strain-specific defaults: H swap pair (`ADS<-BAH`) and synthetic planting strain (`ADS`) fall
  back to a data-chosen strain with a warning when absent; nothing else names strains.
- Report: `quarto render results/report.qmd -P data_dir:../output/<set>`; narrative text is
  about the 10 ScRAP genomes (stated in the README; numbers are computed, prose is not).
- Fix whatever breaks or scales badly at 142 (see check). Expected hot spots: O(n^2) pairwise
  steps (F pi, H trees, D polymorphic-site BLAST all-vs-all), E shuffles, plotly/kaleido PNG
  export without Chrome on a cluster.

### Done during implementation (found on the 12-genome test)
- ScRAP ALH_1b: a 51-bp AT-rich partial CEN6 hit on chrXII (e 7e-6) was extended by the CDE
  annotation to 112 bp and passed the length filter -> two CEN6 loci, `compare_loci` failed.
  `preprocess.primary_loci`: one locus per sample x CEN (best e-value, then longest); others go to
  `secondary_loci.tsv` and the `secondary_loci` QC column.
- Scaling to 1086 (O(n^2) code paths): `trees.k80` as indicator-matrix products (same values,
  1086 x 50k sites in 10 s); H power scan capped by `trees.swap_scan_max_pairs` (200);
  `microevo.pairwise_pi` vectorised over j (same values; ~2.5 min per 1086-sequence window MSA);
  D junction flanks: one BLAST per target genome with all flanks (n runs instead of n^2), majority
  calls by groupby; `transposons.unique_best` vectorised.
- Plotly export without Chrome writes HTML with a warning instead of failing.
- Swap pair / planting strain fall back to data-chosen strains when ADS/BAH are absent.

## .gitignore
Un-ignore: `docs/specs/`, `results/report.qmd`, `results/report_helpers.R`, `results/report.css`,
`data/scrap_assemblies.tsv`, `data/samples_*.tsv`. Keep ignored: assemblies, `data/external/`,
supplementary xlsx (downloaded), `output/`, rendered html, `__pycache__/` (the tracked
`src/centromere/__pycache__/*.pyc` must be `git rm --cached` by you).

## Decisions (rejected alternative)
- One env file, conda (rejected: pip + system tools; the cluster may lack mafft/R).
- Stream archives (rejected: store 16.8 GB first; cluster quota / my 11 GB free disk).
- Collapsed/HP0 + HP1 default (rejected: all haplotypes, which double-counts diploids in
  population statistics).
- Python downloader (rejected: bash + curl | tar; selection rules and sample map are easier
  and testable in Python).
- Report prose stays as is (rejected: rewriting it for an unknown set).

## Check
1. `pytest` passes (new tests: haplotype selection + sample-map writing on a synthetic tar).
2. Linux solve: `CONDA_SUBDIR=linux-64 conda env create --dry-run -f environment.yml`.
Only test-sized downloads on this Mac (user, 2026-11-16). The ScRAP archive already fetched to
the session scratchpad is reused via `--archive` and deleted afterwards; Loegler is never
fetched in full here.
3. Clean checkout: copy `git ls-files` + the new files into a temp dir (no .venv/output),
   `download_data.py --set reference` and `--set scrap142 --archive <scratch> --limit 12`
   (the 10 known strains + 2), then `run_all.sh` to the rendered report. Paste per-step times;
   the full 142/1086 runs and their timing are the student's first cluster job.
4. Loegler: unit test on a synthetic tar with the real member names (`Assemblies_*_fasta/`,
   `X.Final.fasta`, `X.HP1...`); a live `--limit 2` stream stops once 2 members are written
   (only if Final members come early in the archive; otherwise skipped and said so).
5. Already confirmed: ScRAP collapsed/HP0 rule gives 142 strains; ASN collapsed = GCA_949124585.1
   16/16 chromosomes identical (only headers differ: `chrI` vs `ENA|...`).
6. Optional, if Docker Desktop is started: steps 2-3 in a linux/amd64 container.
