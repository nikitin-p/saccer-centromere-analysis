from pathlib import Path

import pandas as pd
import yaml


def load_config(path: str | Path) -> dict:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with open(path) as f:
        return yaml.safe_load(f)


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


def save_tsv(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
