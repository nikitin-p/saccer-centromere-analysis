import difflib

import pandas as pd


def remove_gaps(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["sequence"] = df["sequence_gapped"].str.replace("-", "", regex=False)
    df["sequence_length"] = df["sequence"].str.len()
    return df


def find_split_hits(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    mask = (
        df.groupby(["sample", "centromere", "sseqid"])["sseqid"]
        .transform("size") > 1
    )
    split_df = df[mask].copy()
    ok_df = df[~mask].copy()
    return split_df, ok_df


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
