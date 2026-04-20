"""End-to-end OligoMCP workflow — shared by CLI and MCP server."""
from __future__ import annotations

import datetime as _dt
import json as _json
from concurrent.futures import ThreadPoolExecutor
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
    applied_variants: Optional[list[dict]] = None
    applied_variants_json: Optional[Path] = None


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
    # Populated below when cfg.variants is set; used downstream by
    # SpliceAI's window builder and by BED coord remapping.
    coord_map = None
    parsed_variants: list = []

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

    # --- patient-baseline swap -------------------------------------------
    # If the user has supplied variants, edit them into ref_seq BEFORE ASO
    # enumeration so the downstream pipeline designs against the patient
    # genome. AlphaGenome's `predict_sequence` requires len(sequence) ==
    # interval.width, so indels are re-padded/trimmed from downstream
    # reference; AG's WT baseline prediction is then re-fetched on the
    # patient sequence so `diff_mean_frac`'s ref-vs-alt deltas don't
    # absorb the variant's own splicing effect into every ASO score.
    if cfg.variants:
        from .variants import (
            ExonDeletedByVariant,
            apply_variants_to_ref,
            pad_or_trim_to_length,
            parse_variant,
        )

        reference_seq_length = len(ref_seq)
        parsed_variants = [
            parse_variant(v, gene_symbol=cfg.gene_symbol, assembly=cfg.assembly)
            for v in cfg.variants
        ]
        ref_seq, coord_map = apply_variants_to_ref(
            ref_seq, parsed_variants,
            anchor_genomic=ref_anchor_genomic, chrom=chrom,
        )
        if verbose:
            print(
                f"Applied {len(coord_map.applied)} variant(s) to ref_seq "
                f"(total Δ={coord_map.total_delta():+d} bp):"
            )
            for av in coord_map.applied:
                v = av.variant
                print(
                    f"  - {v.notation[:40]:40s}  "
                    f"{v.chrom}:{v.position} {v.ref or '-'}>{v.alt or '-'}  "
                    f"(Δ={av.delta:+d}, ref_offset={av.ref_offset})"
                )

        ex_start_g, ex_end_g = cfg.exon_intervals
        ex_start_patient = coord_map.ref_to_patient(ex_start_g)
        ex_end_patient = coord_map.ref_to_patient(ex_end_g)
        if ex_start_patient is None or ex_end_patient is None:
            raise ExonDeletedByVariant(
                f"Applied variants remove part of the target exon "
                f"[{ex_start_g}, {ex_end_g}). Nothing to design ASOs against."
            )

        if ag_ctx is not None:
            def _fetch(chrom_, a, b):
                return load_reference_sequence(cfg.fasta_path, chrom_, a, b, assembly=cfg.assembly)

            required = ag_ctx.interval.width
            ref_seq = pad_or_trim_to_length(
                ref_seq,
                target=required,
                fetcher=_fetch,
                chrom=chrom,
                anchor_genomic=ref_anchor_genomic,
                original_length=reference_seq_length,
            )
            ag_ctx.ref_seq = ref_seq
            # Re-fetch WT prediction on the patient baseline.
            if verbose:
                print("Re-fetching AlphaGenome baseline on patient sequence...")
            ag_ctx.ref_output = ag_ctx.model.predict_sequence(
                interval=ag_ctx.interval,
                sequence=ref_seq,
                requested_outputs=ag_ctx.requested_outputs,
                ontology_terms=cfg.ontology_terms,
            )
            # Remap coordinate offsets that were computed in reference
            # space by `setup_alphagenome`. gene_start_rel / gene_end_rel
            # are almost always not at the exon boundaries, so use the
            # coord map; scan_start/end are `variant_interval.start` ±
            # flank and track the exon.
            def _remap(g):
                pat = coord_map.ref_to_patient(g)
                if pat is None:
                    return None
                return pat  # already 0-based patient offset, anchored at ref_anchor_genomic

            ag_ctx.exon_start_rel = ex_start_patient
            ag_ctx.exon_end_rel = ex_end_patient
            gs = _remap(ref_anchor_genomic + ag_ctx.gene_start_rel)
            ge = _remap(ref_anchor_genomic + ag_ctx.gene_end_rel)
            if gs is not None:
                ag_ctx.gene_start_rel = gs
            if ge is not None:
                ag_ctx.gene_end_rel = ge
            # Update scan range (start_rel/end_rel) to match new patient
            # coordinates. variant_interval is flanked around the exon.
            ag_ctx.start_rel = max(0, ex_start_patient - cfg.flank[0])
            ag_ctx.end_rel = min(len(ref_seq), ex_end_patient + cfg.flank[1])
            scan_start_rel = ag_ctx.start_rel
            scan_end_rel = ag_ctx.end_rel
        else:
            # skip_alphagenome path — ref_seq is exon±flank and its
            # length equals scan_end_rel. Re-anchor to patient coords.
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

    # Run AlphaGenome (HTTP-bound) and SpliceAI (torch compute-bound) on
    # parallel threads — their workloads don't compete for the same
    # resource, so wall time drops to roughly max(t_ag, t_sai). Both code
    # paths release the GIL during their heavy work (urllib3 I/O and
    # torch C++ kernels respectively), so plain threading is enough.
    ag_results: dict[str, np.ndarray] = {}
    sai_scores = None

    sai_models = None
    sai_device = None
    if not skip_spliceai:
        from .predict import setup_spliceai

        # Load (or reuse cached) models before dispatching so cache init
        # happens once under the module-level lock rather than racing two
        # worker threads against the same double-checked setup.
        sai_models, sai_device = setup_spliceai(cfg.spliceai_threads)

    def _run_ag() -> dict[str, np.ndarray]:
        if skip_alphagenome or ag_ctx is None:
            return {}
        from .predict import score_asos_alphagenome

        return score_asos_alphagenome(ag_ctx, cfg, candidates)

    def _run_sai():
        if skip_spliceai:
            return None
        from .predict import score_asos_spliceai

        return score_asos_spliceai(
            cfg,
            candidates,
            chrom=chrom,
            variant_interval_start_genomic=ref_anchor_genomic,
            models=sai_models,
            device=sai_device,
            applied_variants=parsed_variants if cfg.variants else None,
            coord_map=coord_map,
        )

    with ThreadPoolExecutor(max_workers=2, thread_name_prefix="oligomcp-score") as ex:
        ag_future = ex.submit(_run_ag)
        sai_future = ex.submit(_run_sai)
        ag_results = ag_future.result()
        sai_scores = sai_future.result()

    for name, arr in ag_results.items():
        all_scores[f"AlphaGenome_{name}"] = arr
    if sai_scores is not None:
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
        coord_map=coord_map,
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

            score_cols = list(all_scores.keys())
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
        coord_map=coord_map,
    )
    if exp_bed:
        bed_files.append(exp_bed)

    applied_variant_records: Optional[list[dict]] = None
    applied_variants_json: Optional[Path] = None
    if coord_map is not None:
        from .variants import applied_variants_to_records

        applied_variant_records = applied_variants_to_records(coord_map)
        applied_variants_json = cfg.results_dir / f"{cfg.config_name}_applied_variants.json"
        applied_variants_json.write_text(
            _json.dumps(
                {
                    "run_id": run_id,
                    "baseline": "patient",
                    "gene_symbol": cfg.gene_symbol,
                    "assembly": cfg.assembly,
                    "applied_variants": applied_variant_records,
                    "patient_seq_length": len(ref_seq),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        if verbose:
            print(f"Wrote {applied_variants_json}")

    return WorkflowResult(
        scores_csv=scores_csv,
        bed_files=bed_files,
        correlation_plot=correlation_png,
        stats=stats,
        n_candidates=len(candidates),
        applied_variants=applied_variant_records,
        applied_variants_json=applied_variants_json,
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
