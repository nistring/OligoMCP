"""Tests for config loading."""
import json
from pathlib import Path

import pytest

from oligomcp.config import load_config


def _write_cfg(tmp_path: Path, **overrides) -> Path:
    base = {
        "gene_symbol": "SETD5",
        "assembly": "hg38",
        "gtf_url": "https://example.com/gtf.feather",
        "results_dir": "results",
        "data_dir": "data",
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
    # 0 is the "auto" sentinel — `score_asos_spliceai` expands this to
    # 256 on CUDA / 64 on CPU at runtime.
    assert cfg.spliceai_batch == 0
    assert cfg.spliceai_threads is None
    assert cfg.alphagenome_workers == 16
    assert cfg.fasta_path is None
    assert cfg.dna_api_key is None


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


def test_results_dir_resolves_relative_to_cwd(tmp_path: Path, monkeypatch):
    """Output dirs should be relative to CWD, not the config file's parent,
    so `config/foo.json` with `results_dir=results` writes to ./results/,
    not ./config/results/."""
    cwd = tmp_path / "workdir"
    cwd.mkdir()
    cfg_dir = tmp_path / "config"
    cfg_dir.mkdir()
    p = _write_cfg(cfg_dir, results_dir="results", data_dir="data")

    monkeypatch.chdir(cwd)
    cfg = load_config(p)

    assert cfg.results_dir == (cwd / "results").resolve()
    assert cfg.data_dir == (cwd / "data").resolve()


def test_results_dir_absolute_passthrough(tmp_path: Path):
    abs_results = tmp_path / "elsewhere" / "out"
    p = _write_cfg(tmp_path, results_dir=str(abs_results))
    cfg = load_config(p)
    assert cfg.results_dir == abs_results.resolve()


def test_missing_opinionated_fields_lists_absent_keys(tmp_path: Path):
    from oligomcp.config import missing_opinionated_fields

    raw = {"gene_symbol": "SETD5", "ASO_length": 20}  # missing aso_step, flank, target_mode
    out = missing_opinionated_fields(raw)
    names = {e["name"] for e in out}
    assert "ASO_length" not in names  # explicit, not missing
    assert {"aso_step", "flank", "target_mode"}.issubset(names)
    # Every entry carries a default and a human-readable description
    for e in out:
        assert "default" in e and "description" in e


def test_missing_opinionated_fields_empty_when_all_present(tmp_path: Path):
    from oligomcp.config import missing_opinionated_fields

    raw = {
        "gene_symbol": "SETD5", "ASO_length": 18, "aso_step": 1,
        "flank": [200, 200], "target_mode": "exclude",
    }
    assert missing_opinionated_fields(raw) == []


def test_variants_field_round_trips(tmp_path: Path):
    """`variants` is opt-in: when present it must round-trip verbatim
    through load_config (parsing is deferred to workflow where
    gene_symbol / assembly are in scope)."""
    p = _write_cfg(
        tmp_path,
        variants=[
            "chr5:70070740:G:A",
            {"chrom": "chr5", "position": 70070700, "ref": "AAG", "alt": "", "id": "demo_del"},
        ],
    )
    cfg = load_config(p)
    assert cfg.variants is not None
    assert len(cfg.variants) == 2
    assert cfg.variants[0] == "chr5:70070740:G:A"
    assert cfg.variants[1]["id"] == "demo_del"
    assert cfg.variants[1]["position"] == 70070700


def test_variants_default_none(tmp_path: Path):
    p = _write_cfg(tmp_path)
    cfg = load_config(p)
    assert cfg.variants is None


def test_variants_wrong_type_rejected(tmp_path: Path):
    p = _write_cfg(tmp_path, variants="chr5:70070740:G:A")   # must be a list
    with pytest.raises(ValueError):
        load_config(p)


def test_variants_items_must_be_str_or_dict(tmp_path: Path):
    p = _write_cfg(tmp_path, variants=[123])
    with pytest.raises(ValueError):
        load_config(p)


def test_variants_not_in_missing_opinionated_fields(tmp_path: Path):
    """variants is purely opt-in and must never appear in the
    needs-info gate (otherwise Claude would be prompted to set it for
    every ASO run)."""
    from oligomcp.config import missing_opinionated_fields

    raw = {"gene_symbol": "SETD5"}
    names = {e["name"] for e in missing_opinionated_fields(raw)}
    assert "variants" not in names


def test_input_paths_still_resolve_relative_to_config(tmp_path: Path, monkeypatch):
    """fasta_path / experimental_data must stay config-dir-relative so input
    bundles shipped next to the JSON keep working regardless of CWD."""
    cfg_dir = tmp_path / "config"
    cfg_dir.mkdir()
    (cfg_dir / "genome.fa").touch()
    (cfg_dir / "exp.csv").touch()
    p = _write_cfg(cfg_dir, fasta_path="genome.fa", experimental_data="exp.csv")

    monkeypatch.chdir(tmp_path)
    cfg = load_config(p)

    assert cfg.fasta_path == (cfg_dir / "genome.fa").resolve()
    assert cfg.experimental_data == (cfg_dir / "exp.csv").resolve()
