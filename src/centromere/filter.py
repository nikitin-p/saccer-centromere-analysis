import pandas as pd


def by_length(df: pd.DataFrame, min_length: int) -> pd.DataFrame:
    return df[df["sequence_length"] >= min_length].copy()


def compute_gc(seq: str) -> float:
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return gc_count / len(seq) if len(seq) > 0 else 0.0


def by_gc(df: pd.DataFrame, min_gc: float) -> pd.DataFrame:
    df = df.copy()
    df["GC"] = df["sequence"].apply(compute_gc)
    return df[df["GC"] > min_gc].copy()
