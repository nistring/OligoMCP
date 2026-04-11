"""End-to-end OligoClaude workflow — shared by CLI and MCP server."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .aso_enum import (
    AsoCandidate,
    enumerate_from_experimental,
    enumerate_sliding,
)
from .bed_export import export_all
from .config import OligoConfig, load_config
from .experimental import aggregate_experimental_candidates, load_experimental
from .plot import correlation_plot, print_correlation_table
from .sequence_utils import load_reference_sequence


@dataclass
class WorkflowResult:
    scores_csv: Optional[Path]
    bed_files: list[Path]
    correlation_plot: Optional[Path]
    stats: Optional[dict]
    ucsc_instructions: str
    n_candidates: int = 0


def _build_ucsc_instructions(assembly: str, bed_files: list[Path]) -> str:
    url = f"https://genome.ucsc.edu/cgi-bin/hgCustom?db={assembly}"
    lines = [
        "=== UCSC Custom Track Upload ===",
        f"1. Open: {url}",
        "2. Click 'Choose File' and select one of the BED files below.",
        "3. Click 'Submit' to load the track.",
        "",
        "Generated BED files:",
    ]
    for p in bed_files:
        lines.append(f"  - {p}")
    return "\n".join(lines)


def run_workflow(
    config_path: Path,
    *,
    skip_alphagenome: bool = False,
    skip_spliceai: bool = False,
    samples_max: int = 20,
    verbose: bool = False,
) -> WorkflowResult:
    """Run the full OligoClaude pipeline and return paths + stats."""
    cfg = load_config(Path(config_path))
    cfg.results_dir = Path(cfg.results_dir) / "ASO"
    cfg.results_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print(f"Config: {config_path}")
        print(
            f"Gene: {cfg.gene_symbol} | Exon: {cfg.exon_intervals} | "
            f"Strand: {cfg.strand} | ASO length: {cfg.ASO_length}"
        )

    all_scores: dict[str, np.ndarray] = {}
    ag_ctx = None
    chrom: Optional[str] = None
    variant_interval_start_genomic: Optional[int] = None
    candidates: list[AsoCandidate] = []
    ref_seq: Optional[str] = None

    if not skip_alphagenome:
        from .alphagenome_predict import score_asos_alphagenome, setup_alphagenome

        ag_ctx = setup_alphagenome(cfg)
        chrom = ag_ctx.interval.chromosome
        variant_interval_start_genomic = ag_ctx.variant_interval.start
        ref_seq = ag_ctx.ref_seq
        start_rel_in_ref = ag_ctx.start_rel
        end_rel_in_ref = ag_ctx.end_rel
    else:
        # Without AlphaGenome context we still need genomic coordinates + ref seq.
        # Load just the scan region from FASTA.
        exon_start, exon_end = cfg.exon_intervals
        variant_interval_start_genomic = exon_start - cfg.flank[0]
        variant_interval_end_genomic = exon_end + cfg.flank[1]
        chrom_guess = _infer_chrom_from_gtf(cfg)
        chrom = chrom_guess
        ref_seq = load_reference_sequence(
            cfg.fasta_path,
            chrom,
            variant_interval_start_genomic,
            variant_interval_end_genomic,
        )
        start_rel_in_ref = 0
        end_rel_in_ref = variant_interval_end_genomic - variant_interval_start_genomic

    if cfg.experimental_data:
        exp_df = load_experimental(cfg.experimental_data)
        candidates = enumerate_from_experimental(
            exp_df,
            ref_seq=ref_seq,
            variant_interval_start_rel=start_rel_in_ref,
            variant_interval_end_rel=end_rel_in_ref,
        )
        if verbose:
            print(
                f"Enumerated {len(candidates)} candidates from experimental data "
                f"({len(exp_df)} rows)"
            )
    else:
        exp_df = None
        candidates = enumerate_sliding(
            ref_seq=ref_seq,
            start_rel=start_rel_in_ref,
            end_rel=end_rel_in_ref,
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
        from .alphagenome_predict import score_asos_alphagenome

        ag_results = score_asos_alphagenome(ag_ctx, cfg, candidates)
        for name, arr in ag_results.items():
            all_scores[f"AlphaGenome_{name}"] = arr

    if not skip_spliceai:
        from .spliceai_predict import score_asos_spliceai, setup_spliceai

        sai_models, _ = setup_spliceai(cfg.spliceai_threads)
        sai_scores = score_asos_spliceai(
            cfg,
            candidates,
            chrom=chrom,
            variant_interval_start_genomic=variant_interval_start_genomic,
            models=sai_models,
        )
        all_scores["SpliceAI"] = sai_scores

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

    scores_csv = cfg.results_dir / f"{cfg.config_name}_ASO_scores.csv"
    scores_df.to_csv(scores_csv, index=False)
    if verbose:
        print(f"Wrote {scores_csv}")

    bed_files = export_all(
        cfg=cfg,
        chrom=chrom,
        variant_interval_start=variant_interval_start_genomic,
        candidates=candidates,
        all_scores=all_scores,
        samples_max=samples_max,
    )

    correlation_png: Optional[Path] = None
    stats: Optional[dict] = None
    if cfg.experimental_data and len(all_scores) > 0:
        # Experimental mode: aggregate candidate scores per experimental row.
        matched = aggregate_experimental_candidates(
            candidates=candidates, scores=all_scores
        )
        if not matched.empty:
            matched_csv = (
                cfg.results_dir / f"{cfg.config_name}_experimental_matched.csv"
            )
            matched.to_csv(matched_csv, index=False)
            if verbose:
                print(f"Wrote {matched_csv}")

            score_cols = list(all_scores.keys())
            correlation_png = cfg.results_dir / f"{cfg.config_name}_correlation.png"
            stats = correlation_plot(
                matched_df=matched,
                score_columns=score_cols,
                measured_col="Measured (RT-PCR)",
                out_png=correlation_png,
            )
            print_correlation_table(stats)

    ucsc = _build_ucsc_instructions(cfg.assembly, bed_files)

    return WorkflowResult(
        scores_csv=scores_csv,
        bed_files=bed_files,
        correlation_plot=correlation_png,
        stats=stats,
        ucsc_instructions=ucsc,
        n_candidates=len(candidates),
    )


def _infer_chrom_from_gtf(cfg: OligoConfig) -> str:
    """Fallback to load chromosome from GTF when AlphaGenome is skipped."""
    try:
        from alphagenome.data import gene_annotation

        gtf = pd.read_feather(cfg.gtf_url)
        gene_interval = gene_annotation.get_gene_interval(
            gtf, gene_symbol=cfg.gene_symbol
        )
        return gene_interval.chromosome
    except Exception as e:
        raise RuntimeError(
            "Could not infer chromosome without AlphaGenome. Add an explicit "
            "`chromosome` field to the config or do not use --skip-alphagenome."
        ) from e
