"""Tests for config loading."""
import json
from pathlib import Path

import pytest

from oligoclaude.config import load_config


def _write_cfg(tmp_path: Path, **overrides) -> Path:
    base = {
        "gene_symbol": "SETD5",
        "assembly": "hg38",
        "gtf_url": "https://example.com/gtf.feather",
        "fasta_path": "data/genome.fa",
        "results_dir": "results",
        "data_dir": "data",
        "dna_api_key": "test_key",
        "track_filter": "polyA plus RNA-seq",
        "ontology_terms": ["CL:0000127"],
        "requested_outputs": ["RNA_SEQ"],
        "exon_intervals": [9429826, 9430051],
        "flank": [200, 200],
        "strand": "+",
        "ASO_length": 18,
    }
    base.update(overrides)
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(base))
    return p


def test_load_basic(tmp_path: Path):
    p = _write_cfg(tmp_path)
    cfg = load_config(p)
    assert cfg.gene_symbol == "SETD5"
    assert cfg.exon_intervals == (9429826, 9430051)
    assert cfg.flank == (200, 200)
    assert cfg.ASO_length == 18
    assert cfg.strand == "+"
    assert cfg.aso_step == 1
    assert cfg.target_mode == "exclude"
    assert cfg.config_name == "cfg"


def test_load_defaults(tmp_path: Path):
    p = _write_cfg(tmp_path)
    cfg = load_config(p)
    assert cfg.experimental_data is None
    assert cfg.spliceai_batch == 12
    assert cfg.spliceai_threads is None


def test_path_resolution_relative(tmp_path: Path):
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "genome.fa").touch()
    p = _write_cfg(tmp_path, fasta_path="data/genome.fa")
    cfg = load_config(p)
    assert cfg.fasta_path.is_absolute()
    assert cfg.fasta_path.name == "genome.fa"


def test_invalid_target_mode_raises(tmp_path: Path):
    p = _write_cfg(tmp_path, target_mode="bogus")
    with pytest.raises(ValueError):
        load_config(p)
