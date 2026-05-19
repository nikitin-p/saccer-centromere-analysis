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
