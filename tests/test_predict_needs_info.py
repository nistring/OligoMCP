"""Tests for the needs_info branch of predict_aso_efficacy.

We don't exercise the full workflow here (that requires network + heavy
deps); we only verify that the tool short-circuits with a structured
needs_info response when opinionated fields are missing from the config."""
from __future__ import annotations

import asyncio
import json
from pathlib import Path

import pytest

from oligomcp.mcp_server import mcp


def _predict_tool():
    return asyncio.run(mcp.get_tools())["predict_aso_efficacy"].fn


def _write_cfg(path: Path, **overrides) -> Path:
    base = {
        "gene_symbol": "SETD5",
        "assembly": "hg38",
        "exon_intervals": [9429826, 9430051],
    }
    base.update(overrides)
    path.write_text(json.dumps(base))
    return path


def test_needs_info_when_fields_omitted(tmp_path: Path):
    cfg = _write_cfg(tmp_path / "cfg.json")  # missing ASO_length, aso_step, flank, target_mode
    r = _predict_tool()(config_path=str(cfg))
    assert r["status"] == "needs_info"
    names = {e["name"] for e in r["missing_fields"]}
    assert {"ASO_length", "aso_step", "flank", "target_mode"}.issubset(names)
    assert r["config_path"].endswith("cfg.json")


def test_confirm_defaults_skips_the_check(tmp_path: Path, monkeypatch):
    """When confirm_defaults=True the tool must NOT short-circuit; it should
    proceed to the workflow instead (we intercept run_workflow to verify
    without actually running it)."""
    cfg = _write_cfg(tmp_path / "cfg.json")
    called = {"n": 0}

    def fake_run_workflow(*args, **kwargs):
        called["n"] += 1
        raise RuntimeError("short-circuit for test")

    monkeypatch.setattr("oligomcp.mcp_server.run_workflow", fake_run_workflow)
    with pytest.raises(RuntimeError, match="short-circuit"):
        _predict_tool()(config_path=str(cfg), confirm_defaults=True)
    assert called["n"] == 1


def test_no_needs_info_when_all_fields_present(tmp_path: Path, monkeypatch):
    cfg = _write_cfg(
        tmp_path / "cfg.json",
        ASO_length=18, aso_step=1, flank=[200, 200], target_mode="exclude",
    )

    def fake_run_workflow(*args, **kwargs):
        raise RuntimeError("short-circuit")

    monkeypatch.setattr("oligomcp.mcp_server.run_workflow", fake_run_workflow)
    with pytest.raises(RuntimeError, match="short-circuit"):
        _predict_tool()(config_path=str(cfg))
