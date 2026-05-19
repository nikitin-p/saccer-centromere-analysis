# Centromere Analysis Refactor Design

Refactor the monolithic `sacer_analysis.ipynb` notebook into a modular Python package with separate scripts and CLI.

## Project Structure

```
saccer-centromere-analysis/
├── src/centromere/
│   ├── __init__.py
│   ├── io.py              # load TSV/FASTA/BED, write outputs, load config
│   ├── preprocess.py      # gap removal, split-hit merging
│   ├── filter.py          # length/GC thresholds
│   ├── align.py           # MAFFT wrapper with multiprocessing
│   ├── mutations.py       # consensus, mutation patterns, base counts
│   └── viz.py             # violin plots, logos, barplots
├── scripts/
│   ├── 01_preprocess.py   # load raw data → merged sequences
│   ├── 02_filter.py       # apply length/GC filters
│   ├── 03_align.py        # run MAFFT per centromere (parallel)
│   ├── 04_analyze.py      # mutation patterns, statistics
│   ├── 05_visualize.py    # generate all figures
│   └── run_pipeline.py    # orchestrator: runs 01-05 in sequence
├── config.yaml            # all parameters
├── notebooks/
│   └── exploration.ipynb  # thin notebook importing from src/
├── requirements.txt
└── pyproject.toml
```

## Data Flow

```
input.tsv (or FASTA/BED)
    │
    ▼
01_preprocess.py
  • io.load_data() detects format, returns DataFrame
  • preprocess.remove_gaps() adds sequence column
  • preprocess.merge_split_hits() stitches overlapping hits
  → output/preprocessed.tsv
    │
    ▼
02_filter.py
  • filter.by_length(min_length)
  • filter.by_gc(min_gc)
  → output/filtered.tsv
    │
    ▼
03_align.py
  • align.run_mafft() per centromere group
  • multiprocessing.Pool for parallel execution
  → output/aligned.tsv
    │
    ▼
04_analyze.py
  • mutations.compute_consensus()
  • mutations.count_bases_per_position()
  • mutations.identify_variable_sites()
  → output/mutation_stats.tsv
  → output/consensus.tsv
    │
    ▼
05_visualize.py
  • viz.plot_length_violin()
  • viz.plot_gc_violin()
  • viz.plot_sequence_logos()
  • viz.plot_gap_histograms()
  • viz.plot_mutation_barplots()
  → output/figures/*.pdf, *.png
```

## Config File

```yaml
input:
  path: "data/1000_centromeres.tsv"
  format: auto  # auto | tsv | fasta | bed

output:
  dir: "output"

preprocess:
  merge_max_distance: 5000

filter:
  min_length: 80
  min_gc: 0.05

align:
  mafft_args: "--auto"
  workers: 4

visualize:
  formats: ["pdf", "png"]
  dpi: 300
  centromere_order: [CEN1, CEN2, CEN3, CEN4, CEN5, CEN6, CEN7, CEN8, CEN9, CEN10, CEN11, CEN12, CEN13, CEN14, CEN15, CEN16]
```

## Module Responsibilities

### io.py
- `load_config(path) -> dict`: Load YAML config
- `load_data(path, format) -> pd.DataFrame`: Detect format (TSV/FASTA/BED), return standardized DataFrame with columns: sample, centromere, qstart, qend, sseqid, sstart, send, sequence_gapped, evalue
- `save_tsv(df, path)`: Write DataFrame to TSV

### preprocess.py
- `remove_gaps(df) -> df`: Add `sequence` column (gaps removed) and `sequence_length`
- `find_split_hits(df) -> df`: Identify rows where same sample/centromere/contig appears multiple times
- `merge_split_hits(df, max_distance) -> df`: Use difflib.SequenceMatcher to stitch overlapping sequences; track non-overlapping hits within max_distance as potential repeats

### filter.py
- `by_length(df, min_length) -> df`: Keep sequences >= min_length
- `by_gc(df, min_gc) -> df`: Compute GC content, keep sequences > min_gc

### align.py
- `run_mafft(sequences, args) -> aligned_sequences`: Run MAFFT on a list of sequences
- `align_centromere_group(df, centromere, args) -> df`: Align all sequences for one centromere
- `align_all_parallel(df, args, workers) -> df`: Parallelize alignment across centromeres using multiprocessing.Pool

### mutations.py
- `compute_consensus(aligned_sequences) -> str`: Per-position majority base
- `count_bases_per_position(aligned_sequences) -> list[dict]`: Base counts at each position
- `identify_variable_sites(base_counts) -> list[int]`: Positions with >1 base present
- `compute_mutation_patterns(df) -> df`: Count substitution types (A>T, A>C, etc.)

### viz.py
- `plot_length_violin(df, output_path)`: Violin plot of sequence lengths per centromere
- `plot_gc_violin(df, output_path)`: Violin plot of GC content per centromere
- `plot_sequence_logos(df, output_dir)`: Per-centromere and global sequence logos using logomaker
- `plot_gap_histograms(df, output_dir)`: Gap frequency per position
- `plot_mutation_barplots(mutation_df, output_dir)`: Stacked barplots of mutation types

## Script Interface

All scripts use argparse with `--config` argument:

```bash
python scripts/01_preprocess.py --config config.yaml
python scripts/02_filter.py --config config.yaml
python scripts/03_align.py --config config.yaml
python scripts/04_analyze.py --config config.yaml
python scripts/05_visualize.py --config config.yaml

# Or run all steps:
python scripts/run_pipeline.py --config config.yaml
```

## Orchestrator

`run_pipeline.py` runs scripts sequentially via subprocess. If any step fails, the pipeline stops. User can fix and re-run from that step manually or re-run the full pipeline.

## Dependencies

```
pandas
numpy
pyyaml
biopython
seaborn
matplotlib
plotly
kaleido
logomaker
```

External: `mafft` must be installed and in PATH.

## Future Extensibility

Structure supports future additions:
- T2T assembly input: extend io.py format detection
- CDE motif analysis: add mutations.py functions or new module
- Flanking region analysis: new script 06_flanking.py
- Transposon detection: new module or integration with external tools
