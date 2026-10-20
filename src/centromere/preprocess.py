import difflib

import pandas as pd


def remove_gaps(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["sequence"] = df["sequence_gapped"].str.replace("-", "", regex=False)
    df["sequence_length"] = df["sequence"].str.len()
    return df


def assign_loci(df: pd.DataFrame) -> pd.DataFrame:
    """Keep, at each genomic locus, only the hits of the best-matching centromere.

    CDEI/CDEIII are shared between centromeres, so a CEN15 query also hits the
    CEN14 locus weakly. Overlapping hits on the same sample/contig form one
    locus, and the locus goes to the centromere with the lowest e-value.
    """
    df = df.copy()
    df["_lo"] = df[["sstart", "send"]].min(axis=1)
    df["_hi"] = df[["sstart", "send"]].max(axis=1)
    df = df.sort_values(["sample", "sseqid", "_lo"]).reset_index(drop=True)

    new_locus = (
        (df["sample"] != df["sample"].shift())
        | (df["sseqid"] != df["sseqid"].shift())
        | (df["_lo"] > df.groupby(["sample", "sseqid"])["_hi"].cummax().shift())
    )
    df["_locus"] = new_locus.cumsum()

    best = df.loc[df.groupby("_locus")["evalue"].idxmin(), ["_locus", "centromere"]]
    df = df.merge(best, on=["_locus", "centromere"])
    return df.drop(columns=["_lo", "_hi", "_locus"]).reset_index(drop=True)


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
        last_merged = group.iloc[0]
        for i in range(1, len(group)):
            row = group.iloc[i]
            prev_row = group.iloc[i - 1]
            # Strand-independent gap between the two hits (negative = overlap)
            distance = (
                max(min(row["sstart"], row["send"]), min(prev_row["sstart"], prev_row["send"]))
                - min(max(row["sstart"], row["send"]), max(prev_row["sstart"], prev_row["send"]))
            )

            if distance <= max_distance:
                merged_seq = _overlap_sequences(merged_seq, row["sequence"])
                last_merged = row

        # Extend coordinates only to pieces that were actually merged; a distant
        # piece (a separate copy) must not stretch this hit across the contig
        first_row = group.iloc[0].to_dict()
        first_row["sequence"] = merged_seq
        first_row["sequence_length"] = len(merged_seq)
        first_row["qend"] = last_merged["qend"]
        first_row["send"] = last_merged["send"]
        results.append(first_row)

    return pd.DataFrame(results)
