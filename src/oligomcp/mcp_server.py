"""MCP server exposing OligoMCP as a Claude / MCP-compatible connector.

Runs locally as a stdio MCP server, registered via:
    pip install -e .
    claude mcp add oligomcp -- oligomcp-mcp

(or the equivalent in Claude Desktop, Cursor, Cline, Continue, Zed,
Goose, etc. — see README.md for per-client configuration.)

Tools exposed:

- `list_gene_exons(gene_symbol, assembly)` — discover exons of a gene's
  canonical transcript (mygene.info, no AlphaGenome key needed).

- `predict_aso_efficacy_inline(gene_symbol, exon_intervals, ...)` — run
  the full workflow from tool arguments (no config file, no disk).
  Returns CSV + BED content inline. Convenient whenever you'd rather
  describe the run in tool arguments than maintain a JSON config.

- `predict_aso_efficacy(config_path, ...)` — traditional file-based
  workflow; reproducible local runs driven from a JSON config.

Both scoring tools emit BED9 files suitable for manual upload to UCSC
Genome Browser's `hgCustom` page, IGV, JBrowse, or any other viewer.

AlphaGenome API key resolution: set $ALPHAGENOME_API_KEY in the server's
environment before launch. SpliceAI-only runs (skip_alphagenome=True)
need no key.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
import threading
from pathlib import Path
from typing import Optional

# Use the `fastmcp` package's FastMCP (v2) rather than `mcp.server.fastmcp`
# (v1): fastmcp's `run` CLI has a bug in its v1-wrapper path where it calls
# the v1 `.run()` (which does `asyncio.run()` internally) from inside an
# already-running async context, producing "Already running asyncio in
# this thread". v2 instances take the native code path and work everywhere.
from fastmcp import FastMCP

from .config import missing_opinionated_fields
from .resources import ENV_VAR, get_alphagenome_api_key
from .workflow import ExonIntervalsRequired, WorkflowResult, run_workflow

_DEFAULT_GTF_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/hg38/"
    "gencode.v46.annotation.gtf.gz.feather"
)
_MAX_CANDIDATES_INLINE = 300

mcp = FastMCP("oligomcp")


def _warm_spliceai_background() -> None:
    """Load SpliceAI models in a background thread at server startup.

    Model load takes a second or two and is the main contributor to
    first-request latency on the SpliceAI path. Running it in a daemon
    thread at import time means:

      - The MCP handshake returns immediately; the host doesn't wait.
      - When the first `predict_aso_efficacy_*` call arrives, the cache
        is usually warm, so scoring kicks off without delay.
      - If the load fails (e.g. bundled weights are missing and no
        network is available), the server still starts; only
        SpliceAI-dependent calls surface the error later. `list_gene_exons`
        and `skip_spliceai=True` runs are unaffected.

    Disable with `OLIGOMCP_PRELOAD_SPLICEAI=0` (useful for tests or
    AlphaGenome-only use where torch import is unwanted overhead).
    """
    def _load() -> None:
        try:
            from .predict import setup_spliceai
            setup_spliceai()
        except Exception as e:  # noqa: BLE001 — log and swallow, don't crash server
            sys.stderr.write(
                f"[oligomcp-mcp] SpliceAI preload failed: {e!r}. "
                "The server is up; SpliceAI-dependent calls will retry "
                "the load on demand.\n"
            )

    threading.Thread(target=_load, daemon=True, name="spliceai-preload").start()


if os.environ.get("OLIGOMCP_PRELOAD_SPLICEAI", "1") != "0":
    _warm_spliceai_background()


# ---------- helpers ----------


def _stats_to_json(stats: Optional[dict]) -> Optional[dict]:
    """Convert nested stats dict (tuples) to JSON-serializable nested dicts."""
    if not stats:
        return None
    out: dict[str, dict] = {}
    for exon, by_src in stats.items():
        out[exon] = {}
        for src, vals in by_src.items():
            r, p, rho, rp = vals
            out[exon][src] = {
                "pearson_r": r,
                "pearson_p": p,
                "spearman_rho": rho,
                "spearman_p": rp,
            }
    return out


def _estimate_candidate_count(
    exon_intervals: list[int], flank: list[int], aso_length: int, aso_step: int
) -> int:
    scan_width = (exon_intervals[1] - exon_intervals[0]) + flank[0] + flank[1]
    return max(0, (scan_width - aso_length) // max(1, aso_step) + 1)


def _scores_from_csv(csv_path: Optional[Path]) -> list[dict]:
    """Read the scores CSV into a list of plain dicts (JSON-serializable)."""
    if not csv_path or not csv_path.exists():
        return []
    import pandas as pd

    df = pd.read_csv(csv_path)
    df = df.where(df.notna(), None)
    return df.to_dict(orient="records")


# ---------- tools ----------


@mcp.tool()
def list_gene_exons(gene_symbol: str, assembly: str = "hg38") -> dict:
    """List all exons of a gene's canonical transcript via mygene.info.

    Use this when the user wants to target a specific gene but hasn't said
    which exon. Present the returned exons to the user (index, coordinates,
    UTR/CDS annotation) and ask which one to exclude/include with an ASO.

    Args:
        gene_symbol: HGNC-style symbol, e.g. "SETD5".
        assembly: Genome assembly ("hg38" default, or "hg19").

    Returns a dict with: gene_symbol, assembly, chrom, strand, transcript,
    cdsstart, cdsend, and `exons` — a list of
    `{index, start, end, length, annotation}` rows.
    """
    from .resources import canonical_transcript_exons, lookup_gene_info

    info = lookup_gene_info(gene_symbol, assembly)
    transcript_id, exons, cdsstart, cdsend = canonical_transcript_exons(info)
    annotated = []
    for i, (a, b) in enumerate(exons, start=1):
        if cdsstart is not None and cdsend is not None:
            if b < cdsstart or a > cdsend:
                anno = "UTR"
            elif a < cdsstart or b > cdsend:
                anno = "CDS-edge"
            else:
                anno = "CDS"
        else:
            anno = "unknown"
        if i == 1:
            anno = f"first/{anno}"
        elif i == len(exons):
            anno = f"last/{anno}"
        annotated.append({
            "index": i,
            "start": a,
            "end": b,
            "length": b - a,
            "annotation": anno,
        })
    return {
        "gene_symbol": gene_symbol,
        "assembly": assembly,
        "chrom": info["chrom"],
        "strand": info["strand"],
        "transcript": transcript_id,
        "cdsstart": cdsstart,
        "cdsend": cdsend,
        "exons": annotated,
    }


@mcp.tool()
def predict_aso_efficacy_inline(
    gene_symbol: str,
    exon_intervals: list[int],
    assembly: str = "hg38",
    strand: str = "+",
    aso_length: int = 18,
    aso_step: int = 5,
    flank: Optional[list[int]] = None,
    target_mode: str = "exclude",
    skip_alphagenome: bool = True,
    skip_spliceai: bool = False,
    samples_max: int = 20,
    requested_outputs: Optional[list[str]] = None,
    ontology_terms: Optional[list[str]] = None,
    track_filter: str = "polyA plus RNA-seq",
) -> dict:
    """Run the full ASO scoring pipeline from arguments, return results inline.

    No config file required: describe the run via tool arguments and the
    server writes CSV/BEDs to a temp directory, reads them back into
    strings, and returns everything in the response payload (including a
    pre-loaded UCSC Genome Browser URL). Use this when driving the tool
    from natural language, or any time maintaining a JSON config isn't
    worth it.

    Safety cap: at most 300 sliding-window ASO candidates per call. If
    you'd get more than that with the given `aso_step` and `flank`, the
    tool returns status="too_many_candidates" with a suggested
    `aso_step` — raise the step (or narrow the flank) and retry.

    Args:
        gene_symbol: HGNC-style symbol, e.g. "SETD5".
        exon_intervals: `[start, end]` of the target exon (genomic, 0-based).
        assembly: "hg38" or "hg19".
        strand: "+" or "-".
        aso_length: ASO length in bp (default 18).
        aso_step: Sliding-window step (default 5; use larger for big exons).
        flank: `[upstream, downstream]` flank around the exon (default 200/200).
        target_mode: "exclude" (skip the exon) or "include".
        skip_alphagenome: Default True — SpliceAI-only is self-contained.
            Set False to also score via AlphaGenome (hosted API, counts
            against your quota; requires $ALPHAGENOME_API_KEY in the
            server environment).
        skip_spliceai: Skip SpliceAI scoring (CPU-heavy ~0.25s/candidate).
        samples_max: Top/bottom ASOs to include in the compact BED (default 20).
        requested_outputs: AlphaGenome output types; defaults to
            ["RNA_SEQ", "SPLICE_SITE_USAGE"].
        ontology_terms: AlphaGenome ontology terms (cell types). Optional.
        track_filter: AlphaGenome track-name substring filter.

    Returns:
        status: "ok" | "too_many_candidates" | "exon_intervals_required"
        n_candidates: int
        scores: list of per-ASO dicts (aso_id, position, length, sequences, scores)
        top_candidates: dict {source: [{aso_id, position, score, sequence}, ...]}
            — samples_max top hits per source.
        bed_tracks: dict {filename: bed_text} (compact BED, for manual upload
            to UCSC / IGV / JBrowse).
    """
    if flank is None:
        flank = [200, 200]
    if requested_outputs is None:
        requested_outputs = ["RNA_SEQ", "SPLICE_SITE_USAGE"]
    if ontology_terms is None:
        ontology_terms = []
    if len(exon_intervals) != 2:
        return {
            "status": "error",
            "message": "exon_intervals must be a 2-element list [start, end].",
        }

    n_est = _estimate_candidate_count(exon_intervals, flank, aso_length, aso_step)
    if n_est > _MAX_CANDIDATES_INLINE:
        scan_width = (exon_intervals[1] - exon_intervals[0]) + flank[0] + flank[1]
        suggested_step = max(1, scan_width // _MAX_CANDIDATES_INLINE + 1)
        return {
            "status": "too_many_candidates",
            "n_estimated": n_est,
            "limit": _MAX_CANDIDATES_INLINE,
            "hint": (
                f"With aso_step={aso_step} you'd scan {n_est} candidates. "
                f"Inline runs are capped at {_MAX_CANDIDATES_INLINE} to keep "
                f"the request under a minute. Retry with aso_step>={suggested_step}."
            ),
        }

    with tempfile.TemporaryDirectory(prefix="oligomcp_") as tmp:
        tmpdir = Path(tmp)
        cfg_data = {
            "gene_symbol": gene_symbol,
            "exon_intervals": list(exon_intervals),
            "assembly": assembly,
            "strand": strand,
            "ASO_length": aso_length,
            "aso_step": aso_step,
            "flank": list(flank),
            "target_mode": target_mode,
            "results_dir": str(tmpdir / "results"),
            "data_dir": str(tmpdir),
            "gtf_url": _DEFAULT_GTF_URL,
            "ontology_terms": ontology_terms,
            "requested_outputs": requested_outputs,
            "track_filter": track_filter,
        }
        cfg_path = tmpdir / "inline_config.json"
        cfg_path.write_text(json.dumps(cfg_data))

        try:
            result: WorkflowResult = run_workflow(
                cfg_path,
                skip_alphagenome=skip_alphagenome,
                skip_spliceai=skip_spliceai,
                samples_max=samples_max,
            )
        except ExonIntervalsRequired as e:
            return {
                "status": "exon_intervals_required",
                "message": str(e),
            }

        scores = _scores_from_csv(result.scores_csv)
        compact_beds = {
            p.name: p.read_text()
            for p in result.bed_files
            if p.exists() and not p.name.endswith("_full.bed")
        }

        score_cols = [
            c for c in (scores[0].keys() if scores else [])
            if c not in {"aso_id", "ASO_sequence", "ASO_antisense", "position",
                         "length", "Measured (RT-PCR)", "Region (Exon)"}
        ]
        top: dict[str, list[dict]] = {}
        for col in score_cols:
            valid = [s for s in scores if s.get(col) is not None]
            ranked = sorted(valid, key=lambda s: -abs(float(s[col])))
            top[col] = [
                {
                    "aso_id": s["aso_id"],
                    "position": s["position"],
                    "score": float(s[col]),
                    "ASO_antisense": s.get("ASO_antisense"),
                }
                for s in ranked[:samples_max]
            ]

        return {
            "status": "ok",
            "n_candidates": result.n_candidates,
            "scores": scores,
            "top_candidates": top,
            "bed_tracks": compact_beds,
        }


@mcp.tool()
def predict_aso_efficacy(
    config_path: str,
    skip_alphagenome: bool = False,
    skip_spliceai: bool = False,
    samples_max: int = 20,
    confirm_defaults: bool = False,
) -> dict:
    """Predict ASO efficacy from a JSON config file.

    Reproducible local workflow — the config file must live on the same
    filesystem as the server (which is always the case for stdio-based
    MCP hosts like Claude Code / Claude Desktop / Cursor). If you'd
    rather pass arguments directly without maintaining a config, use
    `predict_aso_efficacy_inline` instead.

    Args:
        config_path: Absolute path to an OligoMCP JSON config.
        skip_alphagenome: Skip AlphaGenome scoring (no API key needed).
        skip_spliceai: Skip SpliceAI scoring.
        samples_max: Top/bottom ASOs per compact BED track.
        confirm_defaults: When False (default), the tool returns
            `status="needs_info"` if the JSON omits any design-critical
            field (ASO_length, aso_step, flank, target_mode) so Claude can
            ask the user before silently applying a default. Set True to
            skip the check and accept the baked-in defaults.

    Returns:
        status: "ok" | "needs_info" | "exon_intervals_required"
        On "needs_info": `missing_fields` lists `{name, default, description}`
            entries and `config_path` points at the JSON to update.
        On "ok": scores_csv, bed_files, correlation_plot (server-side
            paths), stats, n_candidates.
    """
    cfg_path = Path(config_path).expanduser().resolve()
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config file not found: {cfg_path}")

    if not confirm_defaults:
        raw = json.loads(cfg_path.read_text(encoding="utf-8"))
        missing = missing_opinionated_fields(raw)
        if missing:
            return {
                "status": "needs_info",
                "config_path": str(cfg_path),
                "missing_fields": missing,
                "hint": (
                    "These fields aren't set in the config. Ask the user to "
                    "confirm the default or provide a different value for each "
                    "one, then add them to the JSON and re-invoke. To accept "
                    "all defaults in one shot, pass confirm_defaults=true."
                ),
            }

    try:
        result: WorkflowResult = run_workflow(
            cfg_path,
            skip_alphagenome=skip_alphagenome,
            skip_spliceai=skip_spliceai,
            samples_max=samples_max,
        )
    except ExonIntervalsRequired as e:
        return {
            "status": "exon_intervals_required",
            "message": str(e),
            "hint": (
                "The config is missing `exon_intervals`. Ask the user which "
                "exon to target, then add `\"exon_intervals\": [start, end]` "
                "to the config JSON. You can also call `list_gene_exons(gene_symbol)` "
                "to get the list."
            ),
        }
    return {
        "status": "ok",
        "scores_csv": str(result.scores_csv) if result.scores_csv else None,
        "bed_files": [str(p) for p in result.bed_files],
        "correlation_plot": (
            str(result.correlation_plot) if result.correlation_plot else None
        ),
        "stats": _stats_to_json(result.stats),
        "n_candidates": result.n_candidates,
    }


def _check_startup_credentials() -> None:
    """Print a one-line stderr advisory if no AlphaGenome key is configured."""
    if get_alphagenome_api_key() is None:
        sys.stderr.write(
            "[oligomcp-mcp] WARNING: no AlphaGenome API key in env/file.\n"
            f"  Set ${ENV_VAR} or run `oligomcp set-api-key <KEY>`.\n"
            "  SpliceAI-only runs (skip_alphagenome=true) work without a key.\n"
        )


def main() -> None:
    """Console-script entrypoint: run as a local stdio MCP server."""
    _check_startup_credentials()
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
