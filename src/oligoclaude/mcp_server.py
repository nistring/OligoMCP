"""MCP server exposing OligoClaude as a Claude connector.

Two deployment modes, same tools:

1. Local stdio server (for Claude Desktop / Claude Code):
       pip install -e .
       claude mcp add oligoclaude -- oligoclaude-mcp

2. Remote HTTP server on Prefect Horizon (FastMCP Cloud):
       Sign in at https://horizon.prefect.io with GitHub → select this repo.
       Entrypoint: src/oligoclaude/mcp_server.py:mcp
       Auto-redeploys on push to main. URL: https://<name>.fastmcp.app/mcp

Tools exposed:

- `list_gene_exons(gene_symbol, assembly)` — discover exons of a gene's
  canonical transcript (mygene.info, no AlphaGenome key needed).

- `predict_aso_efficacy_inline(gene_symbol, exon_intervals, ...)` — run the
  full workflow from tool arguments (no config file, no disk). Accepts an
  optional `alphagenome_api_key` per call; returns CSV + BED content +
  UCSC URL inline. Intended for remote HTTP deployments.

- `predict_aso_efficacy(config_path, ...)` — the original file-based tool,
  for local stdio use where both caller and server share a filesystem.

Credential resolution:
    1. `alphagenome_api_key` tool argument (remote — per-request)
    2. $ALPHAGENOME_API_KEY environment variable
    3. ~/.oligoclaude/credentials.json (local only)
    4. Legacy `dna_api_key` field inside the config JSON
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Optional

from mcp.server.fastmcp import FastMCP

from .resources import ENV_VAR, get_alphagenome_api_key
from .workflow import ExonIntervalsRequired, WorkflowResult, run_workflow

_DEFAULT_GTF_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/hg38/"
    "gencode.v46.annotation.gtf.gz.feather"
)
_MAX_CANDIDATES_REMOTE = 300

mcp = FastMCP("oligoclaude")


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
    alphagenome_api_key: Optional[str] = None,
    samples_max: int = 20,
    requested_outputs: Optional[list[str]] = None,
    ontology_terms: Optional[list[str]] = None,
    track_filter: str = "polyA plus RNA-seq",
) -> dict:
    """Run the full ASO scoring pipeline from arguments, return results inline.

    Designed for remote HTTP deployments: no config file required, no
    filesystem access on the caller side. Takes each config field as a
    tool argument; writes CSV/BEDs to a server-side temp directory and
    returns their contents as strings plus a UCSC URL.

    Safety cap: at most 300 sliding-window ASO candidates. If you'd get
    more than that with the given `aso_step` and `flank`, the tool
    returns status="too_many_candidates" with a suggested `aso_step`.

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
            Set False and pass `alphagenome_api_key` to also score via
            AlphaGenome (hosted API, counts against your quota).
        skip_spliceai: Skip SpliceAI scoring (CPU-heavy ~0.25s/candidate).
        alphagenome_api_key: Your AlphaGenome key, used only for this request.
            Ignored when skip_alphagenome=True.
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
        bed_tracks: dict {filename: bed_text} for UCSC upload (compact only).
        ucsc_url: Full URL that preloads the compact BED tracks on UCSC.
        gene_info: {chrom, strand, transcript} derived from mygene.info.
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
    if n_est > _MAX_CANDIDATES_REMOTE:
        scan_width = (exon_intervals[1] - exon_intervals[0]) + flank[0] + flank[1]
        suggested_step = max(1, scan_width // _MAX_CANDIDATES_REMOTE + 1)
        return {
            "status": "too_many_candidates",
            "n_estimated": n_est,
            "limit": _MAX_CANDIDATES_REMOTE,
            "hint": (
                f"With aso_step={aso_step} you'd scan {n_est} candidates. "
                f"Remote runs are capped at {_MAX_CANDIDATES_REMOTE} to keep "
                f"the request under a minute. Retry with aso_step>={suggested_step}."
            ),
        }

    restore_key: Optional[str] = None
    had_key = ENV_VAR in os.environ
    if alphagenome_api_key and not skip_alphagenome:
        restore_key = os.environ.get(ENV_VAR)
        os.environ[ENV_VAR] = alphagenome_api_key

    try:
        with tempfile.TemporaryDirectory(prefix="oligoclaude_") as tmp:
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
                    open_browser=False,
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
                "ucsc_url": result.ucsc_url,
                "ucsc_instructions": result.ucsc_instructions,
            }
    finally:
        if alphagenome_api_key and not skip_alphagenome:
            if restore_key is None and not had_key:
                os.environ.pop(ENV_VAR, None)
            elif restore_key is not None:
                os.environ[ENV_VAR] = restore_key


@mcp.tool()
def predict_aso_efficacy(
    config_path: str,
    skip_alphagenome: bool = False,
    skip_spliceai: bool = False,
    samples_max: int = 20,
) -> dict:
    """Predict ASO efficacy from a JSON config file (local-stdio only).

    For remote HTTP use, call `predict_aso_efficacy_inline` instead — this
    tool requires the config file to exist on the server's filesystem.

    Args:
        config_path: Absolute path to an OligoClaude JSON config on the
            server's filesystem.
        skip_alphagenome: Skip AlphaGenome scoring (no API key needed).
        skip_spliceai: Skip SpliceAI scoring.
        samples_max: Top/bottom ASOs per compact BED track.

    Returns:
        status: "ok" | "exon_intervals_required"
        scores_csv, bed_files, correlation_plot: server-side paths.
        stats, ucsc_instructions, ucsc_url, n_candidates.
    """
    cfg_path = Path(config_path).expanduser().resolve()
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config file not found: {cfg_path}")

    try:
        result: WorkflowResult = run_workflow(
            cfg_path,
            skip_alphagenome=skip_alphagenome,
            skip_spliceai=skip_spliceai,
            samples_max=samples_max,
            open_browser=False,
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
        "ucsc_instructions": result.ucsc_instructions,
        "ucsc_url": result.ucsc_url,
        "n_candidates": result.n_candidates,
    }


def _check_startup_credentials() -> None:
    """Print a one-line stderr advisory if no AlphaGenome key is configured.

    Remote callers can still pass `alphagenome_api_key` per request; this
    check only catches the common local-misconfig case.
    """
    if get_alphagenome_api_key() is None:
        sys.stderr.write(
            "[oligoclaude-mcp] WARNING: no AlphaGenome API key in env/file.\n"
            f"  Local: set ${ENV_VAR} or run `oligoclaude set-api-key <KEY>`.\n"
            "  Remote: pass `alphagenome_api_key` per request to "
            "`predict_aso_efficacy_inline`.\n"
            "  SpliceAI-only runs (skip_alphagenome=true) work without a key.\n"
        )


def main() -> None:
    """Console-script entrypoint: run as a local stdio MCP server."""
    _check_startup_credentials()
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
