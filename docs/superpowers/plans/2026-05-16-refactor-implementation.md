# Centromere Analysis Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the monolithic `sacer_analysis.ipynb` into a modular Python package with separate pipeline scripts.

**Architecture:** Six Python modules in `src/centromere/` provide reusable functions. Six scripts in `scripts/` form the pipeline stages, each reading previous output and writing to `output/`. Config via YAML.

**Tech Stack:** Python 3.10+, pandas, numpy, biopython, pyyaml, seaborn, matplotlib, plotly, logomaker, mafft (external)

**Note:** User handles git commits — skip commit steps.

---

## File Structure

```
src/centromere/
├── __init__.py          # Package init, version
├── io.py                # Config loading, data I/O (TSV/FASTA/BED)
├── preprocess.py        # Gap removal, split-hit merging
├── filter.py            # Length and GC filtering
├── align.py             # MAFFT wrapper with multiprocessing
├── mutations.py         # Consensus, base counts, mutation patterns
└── viz.py               # All visualization functions

scripts/
├── 01_preprocess.py     # Load → merge → save preprocessed.tsv
├── 02_filter.py         # Filter → save filtered.tsv
├── 03_align.py          # Align → save aligned.tsv
├── 04_analyze.py        # Analyze → save mutation_stats.tsv, consensus.tsv
├── 05_visualize.py      # Generate all figures
└── run_pipeline.py      # Orchestrator

tests/
├── __init__.py
├── test_io.py
├── test_preprocess.py
├── test_filter.py
├── test_align.py
├── test_mutations.py
└── test_viz.py
```

---

### Task 1: Project Setup

**Files:**
- Create: `pyproject.toml`
- Create: `requirements.txt`
- Create: `config.yaml`
- Create: `src/centromere/__init__.py`
- Create: `tests/__init__.py`

- [ ] **Step 1: Create pyproject.toml**

```toml
[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"

[project]
name = "centromere"
version = "0.1.0"
description = "S. cerevisiae centromere variability analysis"
requires-python = ">=3.10"
dependencies = [
    "pandas>=2.0",
    "numpy>=1.24",
    "pyyaml>=6.0",
    "biopython>=1.81",
    "seaborn>=0.12",
    "matplotlib>=3.7",
    "plotly>=5.15",
    "kaleido>=0.2",
    "logomaker>=0.8",
]

[project.optional-dependencies]
dev = ["pytest>=7.0"]

[tool.setuptools.packages.find]
where = ["src"]
```

- [ ] **Step 2: Create requirements.txt**

```
pandas>=2.0
numpy>=1.24
pyyaml>=6.0
biopython>=1.81
seaborn>=0.12
matplotlib>=3.7
plotly>=5.15
kaleido>=0.2
logomaker>=0.8
pytest>=7.0
```

- [ ] **Step 3: Create config.yaml**

```yaml
input:
  path: "data/1000_centromeres.tsv"
  format: auto

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
  centromere_order:
    - CEN1
    - CEN2
    - CEN3
    - CEN4
    - CEN5
    - CEN6
    - CEN7
    - CEN8
    - CEN9
    - CEN10
    - CEN11
    - CEN12
    - CEN13
    - CEN14
    - CEN15
    - CEN16
```

- [ ] **Step 4: Create directory structure and __init__.py files**

```bash
mkdir -p src/centromere tests scripts notebooks data output
```

`src/centromere/__init__.py`:
```python
__version__ = "0.1.0"
```

`tests/__init__.py`:
```python
# Test package
```

- [ ] **Step 5: Install package in editable mode**

```bash
pip install -e ".[dev]"
```

---

### Task 2: io.py Module

**Files:**
- Create: `src/centromere/io.py`
- Create: `tests/test_io.py`

- [ ] **Step 1: Write failing test for load_config**

`tests/test_io.py`:
```python
import tempfile
from pathlib import Path

import pytest
import yaml


def test_load_config_returns_dict():
    from centromere.io import load_config

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump({"input": {"path": "test.tsv"}}, f)
        f.flush()
        config = load_config(f.name)

    assert isinstance(config, dict)
    assert config["input"]["path"] == "test.tsv"


def test_load_config_missing_file_raises():
    from centromere.io import load_config

    with pytest.raises(FileNotFoundError):
        load_config("/nonexistent/config.yaml")
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/test_io.py -v
```
Expected: FAIL with "ModuleNotFoundError" or "ImportError"

- [ ] **Step 3: Implement load_config**

`src/centromere/io.py`:
```python
from pathlib import Path

import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with open(path) as f:
        return yaml.safe_load(f)
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tests/test_io.py::test_load_config_returns_dict -v
pytest tests/test_io.py::test_load_config_missing_file_raises -v
```

- [ ] **Step 5: Write failing test for load_data (TSV)**

Add to `tests/test_io.py`:
```python
def test_load_data_tsv():
    from centromere.io import load_data

    tsv_content = "sample1\tCEN1\t1\t100\tcontig1\t500\t600\tATCG--AT\t1e-10\n"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(tsv_content)
        f.flush()
        df = load_data(f.name, format="tsv")

    assert len(df) == 1
    assert list(df.columns) == [
        "sample", "centromere", "qstart", "qend",
        "sseqid", "sstart", "send", "sequence_gapped", "evalue"
    ]
    assert df.iloc[0]["sample"] == "sample1"
    assert df.iloc[0]["sequence_gapped"] == "ATCG--AT"
```

- [ ] **Step 6: Implement load_data**

Add to `src/centromere/io.py`:
```python
COLUMNS = [
    "sample", "centromere", "qstart", "qend",
    "sseqid", "sstart", "send", "sequence_gapped", "evalue"
]


def detect_format(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in (".tsv", ".txt"):
        return "tsv"
    elif suffix in (".fasta", ".fa", ".fna"):
        return "fasta"
    elif suffix == ".bed":
        return "bed"
    else:
        with open(path) as f:
            first_line = f.readline()
        if first_line.startswith(">"):
            return "fasta"
        elif "\t" in first_line:
            return "tsv"
        return "tsv"


def load_data(path: str | Path, format: str = "auto") -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")

    if format == "auto":
        format = detect_format(path)

    if format == "tsv":
        df = pd.read_csv(path, sep="\t", header=None, names=COLUMNS)
    elif format == "fasta":
        df = _load_fasta(path)
    elif format == "bed":
        df = _load_bed(path)
    else:
        raise ValueError(f"Unknown format: {format}")

    return df


def _load_fasta(path: Path) -> pd.DataFrame:
    from Bio import SeqIO

    records = []
    for record in SeqIO.parse(path, "fasta"):
        parts = record.id.split("_")
        sample = parts[0] if len(parts) > 0 else "unknown"
        centromere = parts[1] if len(parts) > 1 else "unknown"
        records.append({
            "sample": sample,
            "centromere": centromere,
            "qstart": 0,
            "qend": len(record.seq),
            "sseqid": "unknown",
            "sstart": 0,
            "send": len(record.seq),
            "sequence_gapped": str(record.seq),
            "evalue": 0.0,
        })
    return pd.DataFrame(records)


def _load_bed(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None)
    result = pd.DataFrame({
        "sample": df.iloc[:, 3] if df.shape[1] > 3 else "unknown",
        "centromere": df.iloc[:, 0],
        "qstart": df.iloc[:, 1],
        "qend": df.iloc[:, 2],
        "sseqid": df.iloc[:, 0],
        "sstart": df.iloc[:, 1],
        "send": df.iloc[:, 2],
        "sequence_gapped": df.iloc[:, 4] if df.shape[1] > 4 else "",
        "evalue": 0.0,
    })
    return result
```

- [ ] **Step 7: Run all io tests**

```bash
pytest tests/test_io.py -v
```

- [ ] **Step 8: Write test for save_tsv and implement**

Add to `tests/test_io.py`:
```python
def test_save_tsv():
    from centromere.io import save_tsv, load_data
    import pandas as pd

    df = pd.DataFrame({
        "sample": ["s1"],
        "centromere": ["CEN1"],
        "value": [42],
    })

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        save_tsv(df, f.name)

    df_read = pd.read_csv(f.name, sep="\t")
    assert len(df_read) == 1
    assert df_read.iloc[0]["sample"] == "s1"
```

Add to `src/centromere/io.py`:
```python
def save_tsv(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
```

- [ ] **Step 9: Run all io tests**

```bash
pytest tests/test_io.py -v
```

---

### Task 3: preprocess.py Module

**Files:**
- Create: `src/centromere/preprocess.py`
- Create: `tests/test_preprocess.py`

- [ ] **Step 1: Write failing test for remove_gaps**

`tests/test_preprocess.py`:
```python
import pandas as pd


def test_remove_gaps():
    from centromere.preprocess import remove_gaps

    df = pd.DataFrame({
        "sequence_gapped": ["ATCG--AT", "GC-TTA--A"],
    })

    result = remove_gaps(df)

    assert "sequence" in result.columns
    assert "sequence_length" in result.columns
    assert result.iloc[0]["sequence"] == "ATCGAT"
    assert result.iloc[0]["sequence_length"] == 6
    assert result.iloc[1]["sequence"] == "GCTTAA"
    assert result.iloc[1]["sequence_length"] == 6
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/test_preprocess.py::test_remove_gaps -v
```

- [ ] **Step 3: Implement remove_gaps**

`src/centromere/preprocess.py`:
```python
import difflib

import pandas as pd


def remove_gaps(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["sequence"] = df["sequence_gapped"].str.replace("-", "", regex=False)
    df["sequence_length"] = df["sequence"].str.len()
    return df
```

- [ ] **Step 4: Run test to verify it passes**

```bash
pytest tests/test_preprocess.py::test_remove_gaps -v
```

- [ ] **Step 5: Write failing test for find_split_hits**

Add to `tests/test_preprocess.py`:
```python
def test_find_split_hits():
    from centromere.preprocess import find_split_hits

    df = pd.DataFrame({
        "sample": ["s1", "s1", "s1", "s2"],
        "centromere": ["CEN1", "CEN1", "CEN2", "CEN1"],
        "sseqid": ["c1", "c1", "c1", "c1"],
        "sstart": [100, 200, 300, 100],
    })

    split_df, ok_df = find_split_hits(df)

    assert len(split_df) == 2  # s1/CEN1/c1 appears twice
    assert len(ok_df) == 2     # s1/CEN2/c1 and s2/CEN1/c1
```

- [ ] **Step 6: Implement find_split_hits**

Add to `src/centromere/preprocess.py`:
```python
def find_split_hits(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    mask = (
        df.groupby(["sample", "centromere", "sseqid"])["sseqid"]
        .transform("size") > 1
    )
    split_df = df[mask].copy()
    ok_df = df[~mask].copy()
    return split_df, ok_df
```

- [ ] **Step 7: Run test to verify it passes**

```bash
pytest tests/test_preprocess.py::test_find_split_hits -v
```

- [ ] **Step 8: Write failing test for merge_split_hits**

Add to `tests/test_preprocess.py`:
```python
def test_merge_split_hits_overlapping():
    from centromere.preprocess import merge_split_hits

    df = pd.DataFrame({
        "sample": ["s1", "s1"],
        "centromere": ["CEN1", "CEN1"],
        "sseqid": ["c1", "c1"],
        "qstart": [0, 20],
        "qend": [30, 50],
        "sstart": [100, 120],
        "send": [130, 150],
        "sequence": ["ATCACATGATATATTTTTTATTTTTAATTTT", "TTTAATTTTTTAATTATAAAAATAATTTTTTT"],
        "sequence_length": [31, 32],
        "sequence_gapped": ["ATCACATGATATATTTTTTATTTTTAATTTT", "TTTAATTTTTTAATTATAAAAATAATTTTTTT"],
        "evalue": [1e-10, 1e-10],
    })

    result = merge_split_hits(df, max_distance=5000)

    assert len(result) == 1
    assert result.iloc[0]["sample"] == "s1"
    assert len(result.iloc[0]["sequence"]) > 31  # merged sequence is longer
```

- [ ] **Step 9: Implement merge_split_hits**

Add to `src/centromere/preprocess.py`:
```python
def _overlap_sequences(s1: str, s2: str) -> str:
    matcher = difflib.SequenceMatcher(None, s1, s2)
    pos_a, pos_b, size = matcher.find_longest_match(0, len(s1), 0, len(s2))
    if size > 0:
        return s1[:pos_a] + s2[pos_b:]
    return s1 + s2


def merge_split_hits(df: pd.DataFrame, max_distance: int = 5000) -> pd.DataFrame:
    results = []

    for (sample, centromere, sseqid), group in df.groupby(
        ["sample", "centromere", "sseqid"]
    ):
        group = group.sort_values("qstart")

        if len(group) == 1:
            results.append(group.iloc[0].to_dict())
            continue

        merged_seq = group.iloc[0]["sequence"]
        for i in range(1, len(group)):
            row = group.iloc[i]
            prev_row = group.iloc[i - 1]
            distance = row["sstart"] - prev_row["send"]

            if distance <= max_distance:
                merged_seq = _overlap_sequences(merged_seq, row["sequence"])

        first_row = group.iloc[0].to_dict()
        first_row["sequence"] = merged_seq
        first_row["sequence_length"] = len(merged_seq)
        first_row["qend"] = group.iloc[-1]["qend"]
        first_row["send"] = group.iloc[-1]["send"]
        results.append(first_row)

    return pd.DataFrame(results)
```

- [ ] **Step 10: Run all preprocess tests**

```bash
pytest tests/test_preprocess.py -v
```

---

### Task 4: filter.py Module

**Files:**
- Create: `src/centromere/filter.py`
- Create: `tests/test_filter.py`

- [ ] **Step 1: Write failing test for by_length**

`tests/test_filter.py`:
```python
import pandas as pd


def test_by_length():
    from centromere.filter import by_length

    df = pd.DataFrame({
        "sequence_length": [50, 80, 100, 79],
    })

    result = by_length(df, min_length=80)

    assert len(result) == 2
    assert list(result["sequence_length"]) == [80, 100]
```

- [ ] **Step 2: Implement by_length**

`src/centromere/filter.py`:
```python
import pandas as pd


def by_length(df: pd.DataFrame, min_length: int) -> pd.DataFrame:
    return df[df["sequence_length"] >= min_length].copy()
```

- [ ] **Step 3: Run test to verify it passes**

```bash
pytest tests/test_filter.py::test_by_length -v
```

- [ ] **Step 4: Write failing test for by_gc**

Add to `tests/test_filter.py`:
```python
def test_by_gc():
    from centromere.filter import by_gc

    df = pd.DataFrame({
        "sequence": ["AAAA", "ATAT", "GCGC", "ATGC"],
    })

    result = by_gc(df, min_gc=0.05)

    assert len(result) == 3  # AAAA filtered out (GC=0)
    assert "GC" in result.columns


def test_gc_content_calculation():
    from centromere.filter import compute_gc

    assert compute_gc("AAAA") == 0.0
    assert compute_gc("GCGC") == 1.0
    assert compute_gc("ATGC") == 0.5
    assert compute_gc("atgc") == 0.5  # case insensitive
```

- [ ] **Step 5: Implement by_gc**

Add to `src/centromere/filter.py`:
```python
def compute_gc(seq: str) -> float:
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return gc_count / len(seq) if len(seq) > 0 else 0.0


def by_gc(df: pd.DataFrame, min_gc: float) -> pd.DataFrame:
    df = df.copy()
    df["GC"] = df["sequence"].apply(compute_gc)
    return df[df["GC"] > min_gc].copy()
```

- [ ] **Step 6: Run all filter tests**

```bash
pytest tests/test_filter.py -v
```

---

### Task 5: align.py Module

**Files:**
- Create: `src/centromere/align.py`
- Create: `tests/test_align.py`

- [ ] **Step 1: Write failing test for run_mafft**

`tests/test_align.py`:
```python
import shutil

import pytest


@pytest.mark.skipif(
    shutil.which("mafft") is None,
    reason="mafft not installed"
)
def test_run_mafft():
    from centromere.align import run_mafft

    sequences = ["ATCGATCG", "ATCGATCG", "ATCGTTCG"]

    aligned = run_mafft(sequences, args="--auto")

    assert len(aligned) == 3
    assert all(len(s) == len(aligned[0]) for s in aligned)
```

- [ ] **Step 2: Implement run_mafft**

`src/centromere/align.py`:
```python
import os
import subprocess
import tempfile
from multiprocessing import Pool

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_mafft(sequences: list[str], args: str = "--auto") -> list[str]:
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = os.path.join(tmpdir, "input.fasta")
        output_path = os.path.join(tmpdir, "output.fasta")

        records = [
            SeqRecord(Seq(seq), id=f"seq{i}", description="")
            for i, seq in enumerate(sequences)
        ]
        SeqIO.write(records, input_path, "fasta")

        cmd = f"mafft {args} {input_path} > {output_path}"
        subprocess.run(cmd, shell=True, check=True, capture_output=True)

        aligned = [str(rec.seq) for rec in SeqIO.parse(output_path, "fasta")]

    return aligned
```

- [ ] **Step 3: Run test (skip if mafft not installed)**

```bash
pytest tests/test_align.py::test_run_mafft -v
```

- [ ] **Step 4: Write test for align_centromere_group**

Add to `tests/test_align.py`:
```python
import pandas as pd


@pytest.mark.skipif(
    shutil.which("mafft") is None,
    reason="mafft not installed"
)
def test_align_centromere_group():
    from centromere.align import align_centromere_group

    df = pd.DataFrame({
        "sample": ["s1", "s2", "s3"],
        "centromere": ["CEN1", "CEN1", "CEN1"],
        "sequence": ["ATCGATCG", "ATCGATCG", "ATCGTTCG"],
    })

    result = align_centromere_group(df, "CEN1", "--auto")

    assert "aligned_seq" in result.columns
    assert len(result) == 3
```

- [ ] **Step 5: Implement align_centromere_group**

Add to `src/centromere/align.py`:
```python
def align_centromere_group(
    df: pd.DataFrame, centromere: str, args: str = "--auto"
) -> pd.DataFrame:
    group = df[df["centromere"] == centromere].copy()
    sequences = group["sequence"].tolist()

    if len(sequences) < 2:
        group["aligned_seq"] = sequences
        return group

    aligned = run_mafft(sequences, args)
    group["aligned_seq"] = aligned
    return group
```

- [ ] **Step 6: Write test for align_all_parallel**

Add to `tests/test_align.py`:
```python
@pytest.mark.skipif(
    shutil.which("mafft") is None,
    reason="mafft not installed"
)
def test_align_all_parallel():
    from centromere.align import align_all_parallel

    df = pd.DataFrame({
        "sample": ["s1", "s2", "s3", "s4"],
        "centromere": ["CEN1", "CEN1", "CEN2", "CEN2"],
        "sequence": ["ATCGATCG", "ATCGTTCG", "GCTAGCTA", "GCTTGCTA"],
    })

    result = align_all_parallel(df, args="--auto", workers=2)

    assert "aligned_seq" in result.columns
    assert len(result) == 4
```

- [ ] **Step 7: Implement align_all_parallel**

Add to `src/centromere/align.py`:
```python
def _align_group_wrapper(args):
    df, centromere, mafft_args = args
    return align_centromere_group(df, centromere, mafft_args)


def align_all_parallel(
    df: pd.DataFrame, args: str = "--auto", workers: int = 4
) -> pd.DataFrame:
    centromeres = df["centromere"].unique()

    tasks = [(df, cen, args) for cen in centromeres]

    with Pool(workers) as pool:
        results = pool.map(_align_group_wrapper, tasks)

    return pd.concat(results, ignore_index=True)
```

- [ ] **Step 8: Run all align tests**

```bash
pytest tests/test_align.py -v
```

---

### Task 6: mutations.py Module

**Files:**
- Create: `src/centromere/mutations.py`
- Create: `tests/test_mutations.py`

- [ ] **Step 1: Write failing test for count_bases_per_position**

`tests/test_mutations.py`:
```python
def test_count_bases_per_position():
    from centromere.mutations import count_bases_per_position

    sequences = ["ATCG", "ATCG", "ATGG"]

    counts = count_bases_per_position(sequences)

    assert len(counts) == 4
    assert counts[0] == {"A": 3}
    assert counts[1] == {"T": 3}
    assert counts[2] == {"C": 2, "G": 1}
    assert counts[3] == {"G": 3}
```

- [ ] **Step 2: Implement count_bases_per_position**

`src/centromere/mutations.py`:
```python
from collections import Counter

import pandas as pd


def count_bases_per_position(sequences: list[str]) -> list[dict]:
    if not sequences:
        return []

    length = len(sequences[0])
    counts = []

    for i in range(length):
        column = [seq[i].upper() for seq in sequences if i < len(seq)]
        counts.append(dict(Counter(column)))

    return counts
```

- [ ] **Step 3: Run test to verify it passes**

```bash
pytest tests/test_mutations.py::test_count_bases_per_position -v
```

- [ ] **Step 4: Write failing test for compute_consensus**

Add to `tests/test_mutations.py`:
```python
def test_compute_consensus():
    from centromere.mutations import compute_consensus

    sequences = ["ATCG", "ATCG", "ATGG"]

    consensus = compute_consensus(sequences)

    assert consensus == "ATCG"  # majority at each position
```

- [ ] **Step 5: Implement compute_consensus**

Add to `src/centromere/mutations.py`:
```python
def compute_consensus(sequences: list[str]) -> str:
    counts = count_bases_per_position(sequences)
    consensus = ""

    for position_counts in counts:
        if position_counts:
            most_common = max(position_counts, key=position_counts.get)
            consensus += most_common
        else:
            consensus += "N"

    return consensus
```

- [ ] **Step 6: Run test to verify it passes**

```bash
pytest tests/test_mutations.py::test_compute_consensus -v
```

- [ ] **Step 7: Write failing test for identify_variable_sites**

Add to `tests/test_mutations.py`:
```python
def test_identify_variable_sites():
    from centromere.mutations import identify_variable_sites

    counts = [
        {"A": 3},           # position 0: not variable
        {"T": 2, "C": 1},   # position 1: variable
        {"G": 3},           # position 2: not variable
        {"A": 1, "T": 1, "G": 1},  # position 3: variable
    ]

    variable = identify_variable_sites(counts)

    assert variable == [1, 3]
```

- [ ] **Step 8: Implement identify_variable_sites**

Add to `src/centromere/mutations.py`:
```python
def identify_variable_sites(base_counts: list[dict]) -> list[int]:
    variable = []
    for i, counts in enumerate(base_counts):
        non_gap_bases = {k: v for k, v in counts.items() if k != "-"}
        if len(non_gap_bases) > 1:
            variable.append(i)
    return variable
```

- [ ] **Step 9: Run test to verify it passes**

```bash
pytest tests/test_mutations.py::test_identify_variable_sites -v
```

- [ ] **Step 10: Write failing test for compute_mutation_patterns**

Add to `tests/test_mutations.py`:
```python
import pandas as pd


def test_compute_mutation_patterns():
    from centromere.mutations import compute_mutation_patterns

    df = pd.DataFrame({
        "centromere": ["CEN1", "CEN1"],
        "aligned_seq": ["ATCG", "ATGG"],
    })

    result = compute_mutation_patterns(df)

    assert "mutation" in result.columns
    assert "count" in result.columns
    assert len(result) > 0
```

- [ ] **Step 11: Implement compute_mutation_patterns**

Add to `src/centromere/mutations.py`:
```python
def compute_mutation_patterns(df: pd.DataFrame) -> pd.DataFrame:
    results = []

    for centromere, group in df.groupby("centromere"):
        sequences = group["aligned_seq"].tolist()
        consensus = compute_consensus(sequences)
        counts = count_bases_per_position(sequences)

        mutation_counts = Counter()

        for i, position_counts in enumerate(counts):
            ref = consensus[i].upper()
            for alt, count in position_counts.items():
                alt = alt.upper()
                if alt != ref and alt != "-" and ref != "-":
                    mutation = f"{ref}>{alt}"
                    mutation_counts[mutation] += count

        for mutation, count in mutation_counts.items():
            results.append({
                "centromere": centromere,
                "mutation": mutation,
                "count": count,
            })

    return pd.DataFrame(results)
```

- [ ] **Step 12: Run all mutations tests**

```bash
pytest tests/test_mutations.py -v
```

---

### Task 7: viz.py Module

**Files:**
- Create: `src/centromere/viz.py`
- Create: `tests/test_viz.py`

- [ ] **Step 1: Write test for plot_length_violin**

`tests/test_viz.py`:
```python
import os
import tempfile

import pandas as pd


def test_plot_length_violin():
    from centromere.viz import plot_length_violin

    df = pd.DataFrame({
        "centromere": ["CEN1"] * 10 + ["CEN2"] * 10,
        "sequence_length": list(range(80, 90)) + list(range(90, 100)),
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "length_violin.png")
        plot_length_violin(df, output_path)
        assert os.path.exists(output_path)
```

- [ ] **Step 2: Implement plot_length_violin**

`src/centromere/viz.py`:
```python
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import seaborn as sns


def plot_length_violin(
    df: pd.DataFrame,
    output_path: str | Path,
    centromere_order: list[str] | None = None,
) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    if centromere_order:
        df = df[df["centromere"].isin(centromere_order)]
        df["centromere"] = pd.Categorical(
            df["centromere"], categories=centromere_order, ordered=True
        )

    fig = px.violin(
        df,
        y="sequence_length",
        x="centromere",
        color="centromere",
        box=False,
        points="all",
    )
    fig.update_traces(orientation="v", side="positive", width=0.5)
    fig.write_image(str(output_path))
    plt.close("all")
```

- [ ] **Step 3: Run test to verify it passes**

```bash
pytest tests/test_viz.py::test_plot_length_violin -v
```

- [ ] **Step 4: Write test for plot_gc_violin and implement**

Add to `tests/test_viz.py`:
```python
def test_plot_gc_violin():
    from centromere.viz import plot_gc_violin

    df = pd.DataFrame({
        "centromere": ["CEN1"] * 10 + ["CEN2"] * 10,
        "GC": [0.3 + i * 0.01 for i in range(20)],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "gc_violin.png")
        plot_gc_violin(df, output_path)
        assert os.path.exists(output_path)
```

Add to `src/centromere/viz.py`:
```python
def plot_gc_violin(
    df: pd.DataFrame,
    output_path: str | Path,
    centromere_order: list[str] | None = None,
) -> None:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    if centromere_order:
        df = df[df["centromere"].isin(centromere_order)]
        df["centromere"] = pd.Categorical(
            df["centromere"], categories=centromere_order, ordered=True
        )

    fig = px.violin(
        df,
        y="GC",
        x="centromere",
        color="centromere",
        box=False,
        points="all",
    )
    fig.update_traces(orientation="v", side="positive", width=0.5)
    fig.write_image(str(output_path))
    plt.close("all")
```

- [ ] **Step 5: Write test for plot_sequence_logos and implement**

Add to `tests/test_viz.py`:
```python
def test_plot_sequence_logos():
    from centromere.viz import plot_sequence_logos

    df = pd.DataFrame({
        "centromere": ["CEN1"] * 5,
        "aligned_seq": ["ATCG", "ATCG", "ATGG", "ATCG", "ATCG"],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        plot_sequence_logos(df, tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "CEN1_logo.png"))
```

Add to `src/centromere/viz.py`:
```python
import logomaker


def plot_sequence_logos(
    df: pd.DataFrame,
    output_dir: str | Path,
    centromere_order: list[str] | None = None,
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    centromeres = centromere_order or df["centromere"].unique()

    for centromere in centromeres:
        group = df[df["centromere"] == centromere]
        if group.empty:
            continue

        sequences = group["aligned_seq"].tolist()
        if not sequences:
            continue

        counts_df = _sequences_to_counts_df(sequences)
        if counts_df.empty:
            continue

        fig, ax = plt.subplots(figsize=(max(10, len(counts_df) * 0.1), 3))
        logomaker.Logo(counts_df, ax=ax)
        ax.set_title(centromere)
        fig.savefig(output_dir / f"{centromere}_logo.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


def _sequences_to_counts_df(sequences: list[str]) -> pd.DataFrame:
    if not sequences:
        return pd.DataFrame()

    length = len(sequences[0])
    bases = ["A", "C", "G", "T", "-"]

    data = {base: [0] * length for base in bases}

    for seq in sequences:
        for i, char in enumerate(seq.upper()):
            if i < length and char in data:
                data[char][i] += 1

    df = pd.DataFrame(data)
    df = df.div(df.sum(axis=1), axis=0).fillna(0)
    return df
```

- [ ] **Step 6: Write test for plot_gap_histograms and implement**

Add to `tests/test_viz.py`:
```python
def test_plot_gap_histograms():
    from centromere.viz import plot_gap_histograms
    from centromere.mutations import count_bases_per_position

    df = pd.DataFrame({
        "centromere": ["CEN1"] * 3,
        "aligned_seq": ["AT-G", "ATCG", "AT-G"],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        plot_gap_histograms(df, tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "CEN1_gaps.png"))
```

Add to `src/centromere/viz.py`:
```python
def plot_gap_histograms(
    df: pd.DataFrame,
    output_dir: str | Path,
) -> None:
    from centromere.mutations import count_bases_per_position

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for centromere, group in df.groupby("centromere"):
        sequences = group["aligned_seq"].tolist()
        counts = count_bases_per_position(sequences)

        gaps = [c.get("-", 0) for c in counts]
        positions = list(range(len(gaps)))

        fig, ax = plt.subplots(figsize=(max(10, len(gaps) * 0.05), 4))
        ax.bar(positions, gaps, width=1.0)
        ax.set_xlabel("Position")
        ax.set_ylabel("Gap count")
        ax.set_title(f"{centromere} - Gap distribution")
        fig.savefig(output_dir / f"{centromere}_gaps.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
```

- [ ] **Step 7: Write test for plot_mutation_barplots and implement**

Add to `tests/test_viz.py`:
```python
def test_plot_mutation_barplots():
    from centromere.viz import plot_mutation_barplots

    mutation_df = pd.DataFrame({
        "centromere": ["CEN1", "CEN1", "CEN2"],
        "mutation": ["A>T", "A>G", "A>T"],
        "count": [10, 5, 8],
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        plot_mutation_barplots(mutation_df, tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "stacked_mutation_barplot.png"))
```

Add to `src/centromere/viz.py`:
```python
def plot_mutation_barplots(
    mutation_df: pd.DataFrame,
    output_dir: str | Path,
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pivot = mutation_df.pivot_table(
        index="mutation",
        columns="centromere",
        values="count",
        fill_value=0,
    )

    fig, ax = plt.subplots(figsize=(12, 6))
    pivot.plot(kind="bar", stacked=True, ax=ax)
    ax.set_xlabel("Mutation type")
    ax.set_ylabel("Count")
    ax.set_title("Mutation patterns by centromere")
    ax.legend(title="Centromere", bbox_to_anchor=(1.05, 1))
    fig.savefig(
        output_dir / "stacked_mutation_barplot.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig)
```

- [ ] **Step 8: Run all viz tests**

```bash
pytest tests/test_viz.py -v
```

---

### Task 8: Pipeline Scripts

**Files:**
- Create: `scripts/01_preprocess.py`
- Create: `scripts/02_filter.py`
- Create: `scripts/03_align.py`
- Create: `scripts/04_analyze.py`
- Create: `scripts/05_visualize.py`
- Create: `scripts/run_pipeline.py`

- [ ] **Step 1: Create 01_preprocess.py**

`scripts/01_preprocess.py`:
```python
#!/usr/bin/env python3
import argparse
from pathlib import Path

from centromere.io import load_config, load_data, save_tsv
from centromere.preprocess import find_split_hits, merge_split_hits, remove_gaps


def main():
    parser = argparse.ArgumentParser(description="Preprocess centromere data")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(config["output"]["dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data from {config['input']['path']}...")
    df = load_data(config["input"]["path"], config["input"]["format"])
    print(f"  Loaded {len(df)} rows")

    print("Removing gaps...")
    df = remove_gaps(df)

    print("Finding split hits...")
    split_df, ok_df = find_split_hits(df)
    print(f"  Found {len(split_df)} split hits, {len(ok_df)} ok hits")

    if len(split_df) > 0:
        print("Merging split hits...")
        max_dist = config["preprocess"]["merge_max_distance"]
        merged_df = merge_split_hits(split_df, max_distance=max_dist)
        print(f"  Merged into {len(merged_df)} rows")

        import pandas as pd
        df = pd.concat([ok_df, merged_df], ignore_index=True)
    else:
        df = ok_df

    output_path = output_dir / "preprocessed.tsv"
    save_tsv(df, output_path)
    print(f"Saved preprocessed data to {output_path} ({len(df)} rows)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Create 02_filter.py**

`scripts/02_filter.py`:
```python
#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.filter import by_gc, by_length
from centromere.io import load_config, save_tsv


def main():
    parser = argparse.ArgumentParser(description="Filter centromere data")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(config["output"]["dir"])

    input_path = output_dir / "preprocessed.tsv"
    print(f"Loading preprocessed data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")
    print(f"  Loaded {len(df)} rows")

    min_length = config["filter"]["min_length"]
    print(f"Filtering by length >= {min_length}...")
    df = by_length(df, min_length)
    print(f"  {len(df)} rows remaining")

    min_gc = config["filter"]["min_gc"]
    print(f"Filtering by GC > {min_gc}...")
    df = by_gc(df, min_gc)
    print(f"  {len(df)} rows remaining")

    output_path = output_dir / "filtered.tsv"
    save_tsv(df, output_path)
    print(f"Saved filtered data to {output_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Create 03_align.py**

`scripts/03_align.py`:
```python
#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.align import align_all_parallel
from centromere.io import load_config, save_tsv


def main():
    parser = argparse.ArgumentParser(description="Align centromere sequences")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(config["output"]["dir"])

    input_path = output_dir / "filtered.tsv"
    print(f"Loading filtered data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")
    print(f"  Loaded {len(df)} rows")

    mafft_args = config["align"]["mafft_args"]
    workers = config["align"]["workers"]
    print(f"Aligning sequences with MAFFT ({workers} workers)...")

    df = align_all_parallel(df, args=mafft_args, workers=workers)
    print(f"  Alignment complete")

    output_path = output_dir / "aligned.tsv"
    save_tsv(df, output_path)
    print(f"Saved aligned data to {output_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Create 04_analyze.py**

`scripts/04_analyze.py`:
```python
#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.io import load_config, save_tsv
from centromere.mutations import (
    compute_consensus,
    compute_mutation_patterns,
    count_bases_per_position,
    identify_variable_sites,
)


def main():
    parser = argparse.ArgumentParser(description="Analyze mutation patterns")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(config["output"]["dir"])

    input_path = output_dir / "aligned.tsv"
    print(f"Loading aligned data from {input_path}...")
    df = pd.read_csv(input_path, sep="\t")
    print(f"  Loaded {len(df)} rows")

    print("Computing mutation patterns...")
    mutation_df = compute_mutation_patterns(df)

    print("Computing consensus sequences...")
    consensus_data = []
    for centromere, group in df.groupby("centromere"):
        sequences = group["aligned_seq"].tolist()
        consensus = compute_consensus(sequences)
        counts = count_bases_per_position(sequences)
        variable = identify_variable_sites(counts)
        consensus_data.append({
            "centromere": centromere,
            "consensus": consensus,
            "length": len(consensus),
            "variable_sites": len(variable),
            "total_sequences": len(sequences),
        })
    consensus_df = pd.DataFrame(consensus_data)

    mutation_path = output_dir / "mutation_stats.tsv"
    save_tsv(mutation_df, mutation_path)
    print(f"Saved mutation stats to {mutation_path}")

    consensus_path = output_dir / "consensus.tsv"
    save_tsv(consensus_df, consensus_path)
    print(f"Saved consensus data to {consensus_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Create 05_visualize.py**

`scripts/05_visualize.py`:
```python
#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd

from centromere.io import load_config
from centromere.viz import (
    plot_gap_histograms,
    plot_gc_violin,
    plot_length_violin,
    plot_mutation_barplots,
    plot_sequence_logos,
)


def main():
    parser = argparse.ArgumentParser(description="Generate visualizations")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = Path(config["output"]["dir"])
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    centromere_order = config["visualize"].get("centromere_order")
    formats = config["visualize"].get("formats", ["png"])

    filtered_path = output_dir / "filtered.tsv"
    aligned_path = output_dir / "aligned.tsv"
    mutation_path = output_dir / "mutation_stats.tsv"

    print("Loading data...")
    filtered_df = pd.read_csv(filtered_path, sep="\t")
    aligned_df = pd.read_csv(aligned_path, sep="\t")
    mutation_df = pd.read_csv(mutation_path, sep="\t")

    for fmt in formats:
        print(f"Generating {fmt.upper()} figures...")

        print("  Length violin plot...")
        plot_length_violin(
            filtered_df,
            figures_dir / f"length_violin.{fmt}",
            centromere_order,
        )

        print("  GC violin plot...")
        plot_gc_violin(
            filtered_df,
            figures_dir / f"gc_violin.{fmt}",
            centromere_order,
        )

    print("Generating sequence logos...")
    plot_sequence_logos(aligned_df, figures_dir / "logos", centromere_order)

    print("Generating gap histograms...")
    plot_gap_histograms(aligned_df, figures_dir / "gap_histograms")

    print("Generating mutation barplots...")
    plot_mutation_barplots(mutation_df, figures_dir / "mutations")

    print(f"All figures saved to {figures_dir}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 6: Create run_pipeline.py**

`scripts/run_pipeline.py`:
```python
#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path


STEPS = [
    "01_preprocess.py",
    "02_filter.py",
    "03_align.py",
    "04_analyze.py",
    "05_visualize.py",
]


def main():
    parser = argparse.ArgumentParser(description="Run full analysis pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument(
        "--start-from",
        type=int,
        default=1,
        help="Start from step N (1-5)",
    )
    args = parser.parse_args()

    scripts_dir = Path(__file__).parent

    for i, script in enumerate(STEPS, start=1):
        if i < args.start_from:
            print(f"Skipping step {i}: {script}")
            continue

        print(f"\n{'='*60}")
        print(f"Step {i}: {script}")
        print(f"{'='*60}\n")

        result = subprocess.run(
            [sys.executable, str(scripts_dir / script), "--config", args.config],
            check=False,
        )

        if result.returncode != 0:
            print(f"\nError: Step {i} failed with exit code {result.returncode}")
            print(f"Fix the issue and re-run with --start-from {i}")
            sys.exit(result.returncode)

    print(f"\n{'='*60}")
    print("Pipeline complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 7: Make scripts executable**

```bash
chmod +x scripts/*.py
```

- [ ] **Step 8: Test pipeline with sample data**

Create a small test TSV and run the pipeline:

```bash
# Create test data
mkdir -p data
echo -e "sample1\tCEN1\t1\t100\tcontig1\t500\t600\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t1e-10" > data/test.tsv
echo -e "sample2\tCEN1\t1\t100\tcontig1\t500\t600\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATTGATCGATCGATCG\t1e-10" >> data/test.tsv

# Update config to use test data
# Then run:
python scripts/run_pipeline.py --config config.yaml
```

---

### Task 9: Exploration Notebook

**Files:**
- Create: `notebooks/exploration.ipynb`

- [ ] **Step 1: Create minimal exploration notebook**

`notebooks/exploration.ipynb`:
```json
{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Centromere Analysis Exploration\n",
    "\n",
    "Interactive exploration using the centromere package."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "from centromere.io import load_config, load_data\n",
    "from centromere.preprocess import remove_gaps, find_split_hits, merge_split_hits\n",
    "from centromere.filter import by_length, by_gc\n",
    "from centromere.align import align_all_parallel\n",
    "from centromere.mutations import compute_consensus, compute_mutation_patterns\n",
    "from centromere.viz import plot_length_violin, plot_gc_violin"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Load config\n",
    "config = load_config('../config.yaml')\n",
    "config"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Load and preprocess data\n",
    "df = load_data(config['input']['path'], config['input']['format'])\n",
    "df = remove_gaps(df)\n",
    "df.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Filter\n",
    "df = by_length(df, config['filter']['min_length'])\n",
    "df = by_gc(df, config['filter']['min_gc'])\n",
    "print(f'Filtered: {len(df)} rows')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Quick stats\n",
    "df.groupby('centromere')['sequence_length'].describe()"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "name": "python",
   "version": "3.10.0"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
```

---

### Task 10: Run All Tests

- [ ] **Step 1: Run full test suite**

```bash
pytest tests/ -v
```

- [ ] **Step 2: Fix any failing tests**

Review output and fix any issues.

- [ ] **Step 3: Verify pipeline runs end-to-end**

```bash
python scripts/run_pipeline.py --config config.yaml
ls -la output/
ls -la output/figures/
```
