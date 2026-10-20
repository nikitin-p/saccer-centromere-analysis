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


def align_centromere_group(
    df: pd.DataFrame, centromere: str, args: str = "--auto", reference: str | None = None
) -> pd.DataFrame:
    group = df[df["centromere"] == centromere].copy()
    sequences = group["sequence"].tolist()

    # The reference is aligned with the strains but stored as a column, not a
    # row, so per-strain statistics downstream never count it as a sample
    if reference is not None:
        sequences = [reference] + sequences

    if len(sequences) < 2:
        aligned = sequences
    else:
        aligned = run_mafft(sequences, args)

    if reference is not None:
        group["ref_aligned_seq"] = aligned[0]
        aligned = aligned[1:]
    group["aligned_seq"] = aligned
    return group


def _align_group_wrapper(args):
    df, centromere, mafft_args, reference = args
    return align_centromere_group(df, centromere, mafft_args, reference)


def align_all_parallel(
    df: pd.DataFrame, args: str = "--auto", workers: int = 4,
    references: dict[str, str] | None = None,
) -> pd.DataFrame:
    centromeres = df["centromere"].unique()
    references = references or {}

    tasks = [(df, cen, args, references.get(cen)) for cen in centromeres]

    with Pool(workers) as pool:
        results = pool.map(_align_group_wrapper, tasks)

    return pd.concat(results, ignore_index=True)
