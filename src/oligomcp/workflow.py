"""End-to-end OligoMCP workflow — shared by CLI and MCP server."""
from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .config import OligoConfig, load_config
from .core import (
    AsoCandidate,
    _strip_match_suffix,
    aggregate_experimental_candidates,
    enumerate_from_experimental,
    enumerate_sliding,
    load_experimental,
    load_reference_sequence,
)
from .output import (
    correlation_plot,
    export_all,
    print_correlation_table,
    write_experimental_bed,
)
from .resources import resolve_fasta_path


@dataclass
class WorkflowResult:
    scores_csv: Optional[Path]
    bed_files: list[Path]
    correlation_plot: Optional[Path]
    stats: Optional[dict]
    n_candidates: int = 0


def run_workflow(
    config_path: Path,
    *,
    skip_alphagenome: bool = False,
    skip_spliceai: bool = False,
    verbose: bool = False,
) -> WorkflowResult:
    """Run the full OligoMCP pipeline and return paths + stats."""
    cfg = load_config(Path(config_path))
    # Nest each run under a timestamped subdirectory so reruns never
    # overwrite prior output. Format: YYYYMMDD_HHMMSS (local time).
    run_id = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    cfg.results_dir = Path(cfg.results_dir) / cfg.gene_symbol / run_id
    cfg.results_dir.mkdir(parents=True, exist_ok=True)

    cfg.fasta_path = resolve_fasta_path(cfg.fasta_path, verbose=verbose)
    _require_exon_intervals(cfg)

    if verbose:
        print(f"Config: {config_path}")
        print(
            f"Gene: {cfg.gene_symbol} | Exon: {cfg.exon_intervals} | "
            f"Strand: {cfg.strand} | ASO length: {cfg.ASO_length}"
        )
        if cfg.fasta_path:
            print(f"FASTA: {cfg.fasta_path}")
        else:
            print("FASTA: online (UCSC API — no local genome download needed)")

    all_scores: dict[str, np.ndarray] = {}
    ag_ctx = None
    chrom: Optional[str] = None
    ref_anchor_genomic: Optional[int] = None
    candidates: list[AsoCandidate] = []
    ref_seq: Optional[str] = None
    scan_start_rel: int = 0
    scan_end_rel: int = 0
    exp_df: Optional[pd.DataFrame] = None

    if not skip_alphagenome:
        from .predict import setup_alphagenome

        ag_ctx = setup_alphagenome(cfg)
        chrom = ag_ctx.interval.chromosome
        ref_seq = ag_ctx.ref_seq
        ref_anchor_genomic = ag_ctx.interval.start
        scan_start_rel = ag_ctx.start_rel
        scan_end_rel = ag_ctx.end_rel
    else:
        exon_start, exon_end = cfg.exon_intervals
        scan_genomic_start = exon_start - cfg.flank[0]
        scan_genomic_end = exon_end + cfg.flank[1]
        chrom = _infer_chrom(cfg)
        ref_seq = load_reference_sequence(
            cfg.fasta_path, chrom, scan_genomic_start, scan_genomic_end,
            assembly=cfg.assembly,
        )
        ref_anchor_genomic = scan_genomic_start
        scan_start_rel = 0
        scan_end_rel = len(ref_seq)

    if cfg.experimental_data:
        exp_df = load_experimental(cfg.experimental_data)
        n_exp_rows = len(exp_df)
        candidates = enumerate_from_experimental(
            exp_df,
            ref_seq=ref_seq,
            variant_interval_start_rel=0,
            variant_interval_end_rel=len(ref_seq),
        )
        matched_ids = {_strip_match_suffix(c.aso_id) for c in candidates}
        n_matched = len(matched_ids)
        if verbose:
            print(
                f"Enumerated {len(candidates)} candidate(s) from "
                f"{n_matched}/{n_exp_rows} experimental row(s) "
                f"(searched entire ref_seq, {len(ref_seq):,} bp)"
            )
        if n_matched < n_exp_rows:
            dropped = n_exp_rows - n_matched
            print(
                f"WARNING: {dropped} experimental ASO(s) did not match the "
                "loaded reference sequence. On `--skip-alphagenome` the scan "
                "region is just exon±flank; run with AlphaGenome to load the "
                "full gene interval."
            )
    else:
        candidates = enumerate_sliding(
            ref_seq=ref_seq,
            start_rel=scan_start_rel,
            end_rel=scan_end_rel,
            aso_length=cfg.ASO_length,
            step=cfg.aso_step,
        )
        if verbose:
            print(
                f"Enumerated {len(candidates)} candidates "
                f"(sliding window, step={cfg.aso_step})"
            )

    if not candidates:
        raise RuntimeError(
            "No ASO candidates produced. Check experimental data and scan region."
        )

    if not skip_alphagenome and ag_ctx is not None:
        from .predict import score_asos_alphagenome

        ag_results = score_asos_alphagenome(ag_ctx, cfg, candidates)
        for name, arr in ag_results.items():
            all_scores[f"AlphaGenome_{name}"] = arr

    if not skip_spliceai:
        from .predict import score_asos_spliceai, setup_spliceai

        sai_models, _ = setup_spliceai(cfg.spliceai_threads)
        sai_results = score_asos_spliceai(
            cfg,
            candidates,
            chrom=chrom,
            variant_interval_start_genomic=ref_anchor_genomic,
            models=sai_models,
        )
        all_scores.update(sai_results)

    scores_df = pd.DataFrame(
        {
            "aso_id": [c.aso_id for c in candidates],
            "ASO_sequence": [c.genomic_target_seq for c in candidates],
            "ASO_antisense": [c.aso_sequence_antisense for c in candidates],
            "position": [c.position for c in candidates],
            "length": [c.length for c in candidates],
        }
    )
    for source, arr in all_scores.items():
        scores_df[source] = arr
    if any(c.measured is not None for c in candidates):
        scores_df["Measured (RT-PCR)"] = [c.measured for c in candidates]
    if any(c.exon_label is not None for c in candidates):
        scores_df["Region (Exon)"] = [c.exon_label for c in candidates]

    scores_csv = cfg.results_dir / f"{cfg.config_name}_scores.csv"
    scores_df.to_csv(scores_csv, index=False)
    if verbose:
        print(f"Wrote {scores_csv}")

    bed_files = export_all(
        cfg=cfg,
        chrom=chrom,
        variant_interval_start=ref_anchor_genomic,
        candidates=candidates,
        all_scores=all_scores,
    )

    correlation_png: Optional[Path] = None
    stats: Optional[dict] = None
    if cfg.experimental_data and len(all_scores) > 0:
        matched = aggregate_experimental_candidates(
            candidates=candidates, scores=all_scores
        )
        if not matched.empty:
            matched_csv = (
                cfg.results_dir / f"{cfg.config_name}_experimental_matched.csv"
            )
            matched.to_csv(matched_csv, index=False)
            if verbose:
                print(
                    f"Wrote {matched_csv} ({len(matched)} experimental rows)"
                )

            # Plot the primary (fractional) columns only — the `_raw`
            # siblings rank identically, so adding them just crowds the
            # panel with collinear lines.
            score_cols = [k for k in all_scores if not k.endswith("_raw")]
            correlation_png = cfg.results_dir / f"{cfg.config_name}_correlation.png"
            stats = correlation_plot(
                matched_df=matched,
                score_columns=score_cols,
                measured_col="Measured (RT-PCR)",
                out_png=correlation_png,
            )
            print_correlation_table(stats)

    exp_bed = write_experimental_bed(
        results_dir=cfg.results_dir,
        config_name=cfg.config_name,
        chrom=chrom,
        strand=cfg.strand,
        candidates=candidates,
        variant_interval_start=ref_anchor_genomic,
    )
    if exp_bed:
        bed_files.append(exp_bed)

    return WorkflowResult(
        scores_csv=scores_csv,
        bed_files=bed_files,
        correlation_plot=correlation_png,
        stats=stats,
        n_candidates=len(candidates),
    )


class ExonIntervalsRequired(RuntimeError):
    """Raised when a run is started without exon_intervals and the user
    must pick one. The message embeds the list of candidate exons from
    mygene.info so the caller (CLI user, or MCP-calling Claude) can present
    choices and retry with an explicit exon."""


def _require_exon_intervals(cfg: OligoConfig) -> None:
    """Validate that the user has specified an exon, listing options if not.

    Does NOT auto-pick an exon — only the user can say which exon to target
    for skipping/inclusion. If exon_intervals is missing, this raises
    `ExonIntervalsRequired` with the list of exons from the canonical
    transcript so the caller can surface them to the user.
    """
    if cfg.exon_intervals is not None:
        return

    from .resources import canonical_transcript_exons, lookup_gene_info

    info = lookup_gene_info(cfg.gene_symbol, cfg.assembly)
    transcript_id, exons, cdsstart, cdsend = canonical_transcript_exons(info)
    lines = [
        f"`exon_intervals` is required but missing.",
        f"Gene:       {cfg.gene_symbol}",
        f"Assembly:   {cfg.assembly}",
        f"Chromosome: {info['chrom']}   Strand: {info['strand']}",
        f"Transcript: {transcript_id}   ({len(exons)} exons)",
    ]
    if cdsstart is not None and cdsend is not None:
        lines.append(f"CDS:        {cdsstart}-{cdsend}")
    lines.append("")
    lines.append("Candidate exons (1-based index, start-end):")
    for i, (a, b) in enumerate(exons, start=1):
        flags = []
        if cdsstart is not None and cdsend is not None:
            if b < cdsstart or a > cdsend:
                flags.append("UTR")
            elif a < cdsstart or b > cdsend:
                flags.append("CDS-edge")
        if i == 1:
            flags.append("first")
        if i == len(exons):
            flags.append("last")
        tag = f"  [{','.join(flags)}]" if flags else ""
        lines.append(f"  {i:>3}. {a}-{b}  (length {b - a} bp){tag}")
    lines.append("")
    lines.append(
        "Add one of these to your config as:  \"exon_intervals\": [start, end]"
    )
    raise ExonIntervalsRequired("\n".join(lines))


def _infer_chrom(cfg: OligoConfig) -> str:
    """Resolve the target chromosome when AlphaGenome is skipped.

    Tries mygene.info first (no heavy deps, handles bare-symbol NL input),
    then falls back to reading the GTF feather via AlphaGenome if installed.
    """
    from .resources import lookup_gene_chromosome

    try:
        return lookup_gene_chromosome(cfg.gene_symbol, cfg.assembly)
    except Exception as mygene_err:
        try:
            from alphagenome.data import gene_annotation

            gtf = pd.read_feather(cfg.gtf_url)
            return gene_annotation.get_gene_interval(
                gtf, gene_symbol=cfg.gene_symbol
            ).chromosome
        except Exception as gtf_err:
            raise RuntimeError(
                f"Could not resolve chromosome for {cfg.gene_symbol!r}. "
                f"mygene.info lookup failed: {mygene_err!r}. "
                f"GTF fallback failed: {gtf_err!r}."
            ) from gtf_err


_infer_chrom_from_gtf = _infer_chrom
