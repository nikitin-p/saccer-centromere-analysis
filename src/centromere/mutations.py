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


def identify_variable_sites(base_counts: list[dict]) -> list[int]:
    variable = []
    for i, counts in enumerate(base_counts):
        non_gap_bases = {k: v for k, v in counts.items() if k != "-"}
        if len(non_gap_bases) > 1:
            variable.append(i)
    return variable


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
