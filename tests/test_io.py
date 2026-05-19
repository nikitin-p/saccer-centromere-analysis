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


def test_save_tsv():
    from centromere.io import save_tsv
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
