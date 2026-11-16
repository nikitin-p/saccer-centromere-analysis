# saccer-centromere-analysis

Analysis of *Saccharomyces cerevisiae* centromere variability and evolution across ~1000 yeast isolates.

## Background

Centromeres are chromosomal regions that recruit kinetochore proteins for chromosome segregation during cell division. *S. cerevisiae* has point centromeres (~125 bp) with three conserved elements:
- **CDEI** (~8 bp) — binds Cbf1
- **CDEII** (~80-90 bp) — AT-rich, wraps around Cse4 nucleosome  
- **CDEIII** (~25 bp) — binds CBF3 complex

This project examines centromere structure and variability in 1011 yeast genomes ([Peter et al. 2018](https://www.nature.com/articles/s41586-018-0030-5)) to understand the "centromeric paradox" — how highly variable sequences perform a conserved function.

## Quick start (cluster)

Everything below runs from the repository root on Linux. Two genome sets are supported:
**ScRAP 142** (O'Donnell et al. 2023, Zenodo [8210254](https://zenodo.org/records/8210254))
and **Loegler 1086** (Loegler et al. 2025, Zenodo [15698884](https://zenodo.org/records/15698884)).
Start with ScRAP: it is small and checks the whole chain.

```bash
# 1. Environment: Python, BLAST+, MAFFT, minimap2, HMMER, R (+ DNAshapeR), Quarto
conda env create -f environment.yml          # mamba is much faster if available
conda activate centromere
pip install -e . --no-deps
pytest -q                                    # unit tests (a few seconds)

# 2. Data
python scripts/download_data.py --set reference      # S288C genome, SGD GFF, Peter/Loegler/Helsen tables (~25 MB)
python scripts/download_data.py --set scrap142       # streams 1.7 GB, keeps 142 assemblies (~0.5 GB)
python scripts/download_data.py --set loegler1086    # streams 16.8 GB, keeps 1086 assemblies (~4 GB)

# 3. Run (or submit: sbatch scripts/slurm_run_all.sh config.scrap142.yaml)
bash scripts/run_all.sh config.scrap142.yaml
# -> output/scrap142/report/report.html
```

**Downloads.** Archives are streamed from Zenodo and never stored; only the selected assemblies are
written, gzipped, to `data/assemblies_<set>/<sample>.fasta.gz`, plus the sample table
`data/samples_<set>.tsv` (`accession assembly strain`). One assembly per strain: ScRAP `collapsed`
where it exists, else `HP0`; Loegler the final assembly from Table S2, `HP1` for phased isolates.
`--haplotypes all` adds the other haplotypes as extra samples (they double-count diploids in
population statistics). An interrupted download can be restarted; finished files are skipped. When
the whole archive was read, its md5 is checked against Zenodo. For a quick test:
`--limit 10` or `--strains ADS,CFS,...`; with a downloaded tarball: `--archive nuclear.tar.gz`.

**Run.** `run_all.sh` runs every step in order and logs each to `output/<set>/logs/<step>.log`;
`bash scripts/run_all.sh CONFIG --list` shows them. A finished step is skipped on rerun;
`--from STEP` redoes STEP and everything after it, `--only STEP` just that step. Optional steps that
need extra external data are off: `--with validate` (analysis A: Helsen 2026 SD4 + PCAn),
`--with peter748,phenotypes` (G and phenotypes: Peter 2018 set, `config.peter748.yaml`).
Programs come from PATH (config `blast_bin`, `minimap2`, `hmmer_bin: null`); `THREADS` sets BLAST
threads, `align.workers` and `locate.threads` in the config the rest.

**Report.** Rendered into `output/<set>/report/` from `results/report.qmd`. All numbers, tables
and figures are computed from that set, but the prose was written for the first 10 ScRAP genomes
(strain names such as ADS, CFA, BCE in the text) and has to be re-read against the new results.

**Known limits.**
- Plotly PNG/PDF export (step `pipeline`, `05_visualize`) needs Chrome for kaleido; without it
  the figures are written as HTML and the step continues (`plotly_get_chrome` installs one).
- Disk: BLAST DBs and uncompressed copies for analyses C/D/E take ~25 MB per genome
  (~30 GB for Loegler 1086).
- The H swap check uses ADS<-BAH; if a set lacks them, the most concordant strain and its most
  genome-distant donor are used. The H power scan samples at most `trees.swap_scan_max_pairs`
  ordered pairs (200).
- A second locus of a CEN in one genome (usually a short AT-rich partial BLAST hit) is set aside
  in `output/<set>/secondary_loci.tsv`; only the best locus per sample and CEN is analysed.

## Installation (local, without conda)

```bash
python -m venv .venv && .venv/bin/pip install -e ".[dev]" openpyxl
# plus BLAST+, MAFFT, minimap2, HMMER, R with the packages in environment.yml, Quarto
```

## Usage

### Full Pipeline

```bash
# Edit config.yaml to set your input data path
python scripts/run_pipeline.py --config config.yaml
```

### Individual Steps

```bash
python scripts/01_preprocess.py --config config.yaml  # Load, remove gaps, merge split hits
python scripts/02_filter.py --config config.yaml      # Filter by length and GC content
python scripts/03_align.py --config config.yaml       # MAFFT alignment per centromere
python scripts/04_analyze.py --config config.yaml     # Mutation patterns, consensus
python scripts/05_visualize.py --config config.yaml   # Generate figures
```

### As a Library

```python
from centromere.io import load_config, load_data
from centromere.preprocess import remove_gaps
from centromere.filter import by_length, by_gc

config = load_config('config.yaml')
df = load_data(config['input']['path'])
df = remove_gaps(df)
df = by_length(df, min_length=80)
df = by_gc(df, min_gc=0.05)
```

## Configuration

Edit `config.yaml`:

```yaml
input:
  path: "data/1000_centromeres.tsv"
  format: auto  # auto | tsv | fasta | bed

filter:
  min_length: 80
  min_gc: 0.05

align:
  mafft_args: "--auto"
  workers: 4
```

## Input Data Format

TSV with columns: `sample`, `centromere`, `qstart`, `qend`, `sseqid`, `sstart`, `send`, `sequence_gapped`, `evalue`

Generated by BLASTN alignment of centromere sequences to genome assemblies.

## Output

```
output/
├── preprocessed.tsv      # After gap removal and split-hit merging
├── filtered.tsv          # After length/GC filtering
├── aligned.tsv           # With aligned_seq column
├── mutation_stats.tsv    # Mutation type counts per centromere
├── consensus.tsv         # Consensus sequences and variable site counts
└── figures/
    ├── length_violin.pdf
    ├── gc_violin.pdf
    ├── logos/            # Sequence logos per centromere
    ├── gap_histograms/   # Gap distribution plots
    └── mutations/        # Mutation pattern barplots
```

## Testing

```bash
pytest tests/ -v
```

## Project Structure

```
├── src/centromere/       # Core library
│   ├── io.py             # Config and data I/O
│   ├── preprocess.py     # Gap removal, split-hit merging
│   ├── filter.py         # Length and GC filtering
│   ├── align.py          # MAFFT wrapper with multiprocessing
│   ├── mutations.py      # Consensus and mutation analysis
│   └── viz.py            # Visualization functions
├── scripts/              # Pipeline scripts
├── tests/                # Unit tests
├── notebooks/            # Interactive exploration
└── config.yaml           # Pipeline configuration
```

## Data Sources

### Current Data

The primary dataset (`1000_centromeres.tsv`) is derived from the [1011 yeast genomes](https://www.nature.com/articles/s41586-018-0030-5) (Peter et al. 2018).

### Sample Metadata

Two supplementary files in `data/` provide strain metadata:

| File | Source | Content |
|------|--------|---------|
| `41586_2018_30_MOESM3_ESM.xls` | Peter et al. 2018 | Ecological/geographical origins, phylogenetic clades (26 clades including Wine/European, Sake, Brazilian bioethanol) |
| `41586_2025_9637_MOESM3_ESM.xlsx` | Loegler et al. 2025 | Sequencing stats, assembly quality (MerquryQV, BUSCO scores), filenames |

**949 samples (87%)** can be linked between datasets via `Standardized name`, enabling analyses like:
- Comparing centromere variation across ecological niches (wine vs. clinical vs. wild)
- Geographic structure of centromere variants
- Clade-specific centromere signatures

### Potential Additional Sources

**T2T Assemblies for S. cerevisiae:**
- **1086 strains** — [Loegler et al. 2025](https://www.nature.com/articles/s41586-025-09637-0): Near-T2T assemblies with pangenome and structural variant atlas
- **142 strains (ScRAP)** — [O'Donnell et al. 2023](https://www.nature.com/articles/s41588-023-01459-y): Reference-quality T2T assemblies

**Other Saccharomyces Species:**
- [Scannell et al. 2011](https://academic.oup.com/g3journal/article/1/1/11/5986445): *S. sensu stricto* (5 species) — SRP006340
- [Kellis et al. 2003](https://www.nature.com/articles/nature01644): Comparative genomics of 4 species

**Centromere evolution — key literature:**
- [Helsen et al. 2026, *Nature*](https://doi.org/10.1038/s41586-025-09779-1): Progressive coevolution of the yeast centromere and kinetochore. Point-centromere atlas for 138 Saccharomycetaceae species and 1,493 *S. cerevisiae* strains; ~9.5% of strains carry CDEII-length variants (13 independent variants, spread through sex, mostly microhomology-mediated duplications); ~10-bp CDEII jumps are tolerated, 5-bp ones are not; Cbf1 coevolves with transitions. Annotation tool **PCAn**: [Zenodo](https://doi.org/10.5281/zenodo.17293587).
- [Haase et al. 2026, *Nature*](https://doi.org/10.1038/s41586-025-10092-0): Ancient co-option of LTR retrotransposons as yeast centromeres. "Proto-point" centromeres in Saccharomycodales; point centromeres (CDEI, CDEII, CDEIII) derived from Ty5 LTR sequences; CBF3 complex predates the tripartite CDE structure (110–140 Myr); rejects the 2µ-plasmid origin.
- [Włodzimierz et al. 2026, *bioRxiv* (preprint)](https://doi.org/10.64898/2026.09.17.752431): An atlas of eukaryotic centromere architecture. 325 Darwin Tree of Life genomes; centromeres cycle between satellite (α) and transposon-cluster (Ω) states; centrophilic transposons invade satellite arrays and are recombined back into tandem satellite-like repeats; satellite homogenization is mostly within-chromosome.
- [Hession et al. 2025, *PLoS Genet*](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011814): Centromeres in budding yeasts are conserved in chromosomal location but not in structure. Hi-C + Cse4 ChIP-seq in 9 Saccharomycodales/Phaffomycetales species; structures turn over in situ while synteny around centromeres is kept for ~200 Myr.

**ChIP-seq / Chromatin:**
- [Henikoff 2012](https://www.pnas.org/doi/10.1073/pnas.1118898109): Tripartite centromeric chromatin organization
- [Zhang et al. 2024](https://academic.oup.com/genetics/article/228/1/iyae108/7710407): Cse4/CENP-A localization

**Transposon Analysis:**
- [Yue et al. 2017](https://www.nature.com/articles/ng.3847): Domesticated vs wild yeast genomes, Ty elements
- [Dunn et al. 2012](https://genome.cshlp.org/content/22/5/908): Pan-genome copy number variants

**Phenotype Data:**
- `Phenotypes_8391Traits/` from [Loegler et al. 2025](https://www.nature.com/articles/s41586-025-09637-0): Growth traits (~241 conditions including DNA damage), RNA-seq (~6000 genes), protein abundance (~2148 proteins). Strain IDs match FASTA headers. Potential to correlate centromere variants with fitness under genotoxic stress or kinetochore gene expression.

## License

MIT
