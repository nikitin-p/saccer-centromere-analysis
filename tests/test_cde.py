import gzip

import pandas as pd

from centromere.cde import CDEI_LEN, CDEIII_LEN, CDEModel, annotate_cde, read_reference

REFERENCE = read_reference("data/sacer_centromere.fasta")


def _revcomp(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def test_model_recovers_every_s288c_element():
    model = CDEModel(REFERENCE)
    flank = "GCGCGCGCGCGCGCGCGCGC" * 2
    for cen, seq in REFERENCE.items():
        window = flank + seq + flank
        cde = model.annotate(window, len(flank), len(flank) + len(seq))
        assert cde["cdei_seq"] == seq[:CDEI_LEN], cen
        assert cde["cdeiii_seq"] == seq[-CDEIII_LEN:], cen
        assert cde["cde_ok"], cen


def _write_genome(tmp_path, contig_seq):
    path = tmp_path / "s1.fasta.gz"
    with gzip.open(path, "wt") as f:
        f.write(f">ctg1 test assembly, contig: chrXVI\n{contig_seq}\n")
    return tmp_path


def test_annotate_keeps_cdeii_insertion_on_minus_strand(tmp_path):
    cen16 = REFERENCE["CEN16"]
    dup = cen16[36:58]  # duplicate 22 bp inside CDEII, as in CFS CEN16
    strain_cen = cen16[:58] + dup + cen16[58:]
    left, right = "GC" * 300, "CG" * 300
    genome_dir = _write_genome(tmp_path, left + _revcomp(strain_cen) + right)

    # BLAST-like hit on the minus strand covering only query 1..117 minus the dup
    hit_lo = len(left) + 1
    hit_hi = len(left) + len(strain_cen)
    df = pd.DataFrame([{
        "sample": "s1", "centromere": "CEN16", "sseqid": "ctg1",
        "qstart": 1, "qend": len(cen16), "sstart": hit_hi, "send": hit_lo,
        "sequence": cen16, "sequence_length": len(cen16),
    }])

    result = annotate_cde(df, genome_dir, REFERENCE).iloc[0]

    assert result["sequence"] == strain_cen
    assert result["cdeii_length"] == len(cen16) - CDEI_LEN - CDEIII_LEN + len(dup)
    assert result["strand"] == "-"
    assert result["contig_name"] == "chrXVI"
    assert (result["cen_start"], result["cen_end"]) == (hit_lo, hit_hi)
