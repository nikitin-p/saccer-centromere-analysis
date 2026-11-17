"""scripts/download_data.py: assembly selection, strain names and tar extraction (no network)."""
import gzip
import importlib.util
import io
import tarfile
from pathlib import Path

import pandas as pd
import pytest

SPEC = importlib.util.spec_from_file_location(
    "download_data", Path(__file__).resolve().parent.parent / "scripts" / "download_data.py")
dd = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dd)

# Real member-name shapes (ScRAP nuclear.tar.gz, 2026-11-16 listing)
SCRAP = [
    "nuclear/ASN.asm01.collapsed.nuclear_genome.tidy.fa.gz",
    "nuclear/ASN.asm01.HP1.nuclear_genome.tidy.fa.gz",
    "nuclear/ASN.asm01.HP2.nuclear_genome.tidy.fa.gz",
    "nuclear/ANM.asm01.HP.nuclear_genome.tidy.fa.gz",
    "nuclear/ANM.asm01.collapsed.nuclear_genome.tidy.fa.gz",
    "nuclear/CCC_1a.asm01.HP0.nuclear_genome.tidy.fa.gz",
    "nuclear/CCC_1a.asm01.HP0.nuclear_genome.tidy.gff3.gz",
    "nuclear/CCC_1a.asm01.HP0.nuclear_genome.tidy.cds.fa.gz",
    "nuclear/CGH_1.asm01.HP0.nuclear_genome.tidy.fa.gz",
    "nuclear/CGH_3.asm01.HP0.nuclear_genome.tidy.fa.gz",
]


def test_scrap_one_per_strain():
    t = dd.scrap_targets(SCRAP)
    assert {v["sample"] for v in t.values()} == {"ASN_collapsed", "ANM_collapsed", "CCC_1a", "CGH_1", "CGH_3"}
    assert "nuclear/CCC_1a.asm01.HP0.nuclear_genome.tidy.cds.fa.gz" not in t   # CDS/pep/gff skipped


def test_scrap_all_haplotypes():
    t = dd.scrap_targets(SCRAP, "all")
    assert {"ASN_HP1", "ASN_HP2"} <= {v["sample"] for v in t.values()}
    assert "ANM_HP" not in {v["sample"] for v in t.values()}   # unnamed "HP" files never used


def test_strain_names():
    s = dd.strain_names(dd.scrap_targets(SCRAP, "all"), "scrap142")
    assert s["CCC_1a"] == "CCC"                       # Peter 2018 key
    assert s["ASN_collapsed"] == "ASN"
    assert s["CGH_1"] == "CGH_1" and s["CGH_3"] == "CGH_3"   # suffix kept: would collide
    assert s["ASN_HP2"] == "ASN_HP2"                  # extra haplotype = own strain key
    assert len(set(s.values())) == len(s)


def test_strain_names_avoid_reference_name():
    t = dd.scrap_targets(["x/S288C.asm01.HP0.nuclear_genome.tidy.fa.gz"])
    assert dd.strain_names(t, "scrap142") == {"S288C": "S288C_asm"}   # "S288C" = R64 control


def loegler_table():
    return pd.DataFrame({"Standardized_name": ["AAB", "AAL", "AAL", "AGA", "AGA_1a"],
                         "Haplotype": [None, "HP1", "HP2", None, None],
                         "Filename": ["AAB.Final.fasta", "AAL.HP1.Final.fasta", "AAL.HP2.Final.fasta",
                                      "AGA.Final.fasta", "AGA_1a.Final.fasta"]})


def test_loegler_targets_and_strains():
    t = dd.loegler_targets(loegler_table())
    assert set(t) == {"AAB.Final.fasta", "AAL.HP1.Final.fasta", "AGA.Final.fasta", "AGA_1a.Final.fasta"}
    s = dd.strain_names(t, "loegler1086")
    assert s == {"AAB": "AAB", "AAL_HP1": "AAL", "AGA": "AGA", "AGA_1a": "AGA_1a"}   # AGA != AGA_1a
    assert "AAL.HP2.Final.fasta" in dd.loegler_targets(loegler_table(), "all")


def make_tar(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


def test_extract_loegler_layout(tmp_path):
    fa = b">chrI\nACGT\n"
    tar = make_tar({"GenomeAssemblies/Assemblies_SecondBest_fasta/AAB.SecondBest.fasta": b">x\nN\n",
                    "GenomeAssemblies/Assemblies_Final_fasta/AAB.Final.fasta": fa,
                    "GenomeAssemblies/Assemblies_Final_fasta/AAL.HP1.Final.fasta": fa,
                    "GenomeAssemblies/Assemblies_Final_fasta/AAL.HP2.Final.fasta": fa})
    targets = dd.loegler_targets(loegler_table())
    key = lambda n: Path(n).name if "/Assemblies_SecondBest" not in n else None  # noqa: E731
    done, full = dd.extract(io.BytesIO(tar), targets, key, tmp_path, None)
    assert {t["sample"] for t in done} == {"AAB", "AAL_HP1"}
    assert full                                                   # AGA missing -> read to the end
    assert gzip.decompress((tmp_path / "AAB.fasta.gz").read_bytes()) == fa
    assert not list(tmp_path.glob("*.tmp"))
    # restart: nothing rewritten, both reported as present
    done2, _ = dd.extract(io.BytesIO(tar), targets, key, tmp_path, None)
    assert {t["sample"] for t in done2} == {"AAB", "AAL_HP1"}


def test_extract_limit_stops_early(tmp_path):
    fa = gzip.compress(b">chrI\nACGT\n")
    tar = make_tar({m: fa for m in SCRAP if m.endswith("tidy.fa.gz")})
    targets = dd.scrap_targets(SCRAP)
    done, full = dd.extract(io.BytesIO(tar), targets, lambda n: n, tmp_path, 2)
    assert len(done) == 2 and not full
    assert len(list(tmp_path.glob("*.fasta.gz"))) == 2
    assert gzip.decompress((tmp_path / f"{done[0]['sample']}.fasta.gz").read_bytes()) == b">chrI\nACGT\n"


def test_apply_filter():
    t = dd.scrap_targets(SCRAP)
    assert {v["sample"] for v in dd.apply_filter(t, ["ASN", "CCC_1a"]).values()} == {"ASN_collapsed", "CCC_1a"}
    with pytest.raises(SystemExit):
        dd.apply_filter(t, ["NOPE"])
