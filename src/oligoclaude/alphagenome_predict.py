"""AlphaGenome prediction — ports helpers from AlphaGenome_ASO/aso.ipynb.

Original helpers (parse_output_types, filter_td, diff_mean, etc.) are copied
from `/home/nistring/AlphaGenome_ASO/aso.ipynb` to keep OligoClaude self-contained.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from .aso_enum import AsoCandidate
from .config import OligoConfig
from .sequence_utils import load_reference_sequence


@dataclass
class AGContext:
    model: Any
    gene_interval: Any
    interval: Any
    ref_seq: str
    variant_interval: Any
    ref_output: Any
    requested_outputs: list[Any]
    exon_start_rel: int
    exon_end_rel: int
    gene_start_rel: int
    gene_end_rel: int
    start_rel: int
    end_rel: int


def parse_output_types(output_type_names: list[str]) -> list[Any]:
    """Parse output type names from config into AlphaGenome OutputType enums."""
    from alphagenome.models.dna_client import OutputType

    output_types: list[Any] = []
    for name in output_type_names:
        name_upper = name.upper()
        if name_upper in {"SPLICE_JUNCTIONS", "SPLICE_SITES"}:
            raise ValueError(
                f"{name_upper} is not supported. Use SPLICE_SITE_USAGE instead."
            )
        if name_upper in OutputType.__members__:
            output_types.append(OutputType[name_upper])
        else:
            print(f"Warning: Unknown output type '{name}', skipping...")
    return output_types


def name_filter(names: list[str], substring: str) -> list[bool]:
    return [substring in n for n in names]


def filter_by_strand(track_data, strand: str):
    """Apply strand selection to a TrackData object."""
    s = str(strand).lower()
    mapping = {
        "+": track_data.filter_to_positive_strand,
        "-": track_data.filter_to_negative_strand,
        "nonnegative": track_data.filter_to_nonnegative_strand,
        "nonpositive": track_data.filter_to_nonpositive_strand,
        "stranded": track_data.filter_to_stranded,
        "unstranded": track_data.filter_to_unstranded,
    }
    if s not in mapping:
        return track_data
    return mapping[s]()


def filter_td(td, cfg: OligoConfig):
    if td is None:
        return None
    if cfg.track_filter:
        td = td.filter_tracks(name_filter(td.names, cfg.track_filter))
    if cfg.strand:
        td = filter_by_strand(td, cfg.strand)
    return td


def diff_mean(
    ref: np.ndarray,
    alt: np.ndarray,
    start: int,
    end: int,
    gene_start: int,
    gene_end: int,
) -> np.ndarray:
    """Difference of means in [start:end], normalized by gene body mean.

    Port of `diff_mean` from AlphaGenome_ASO/aso.ipynb. Both `ref` and `alt`
    are arrays of shape (seq_length, num_variants) (or (seq_length, 1) for ref).
    """
    return (
        alt[start:end].mean(0) / ref[gene_start:gene_end].mean(0)
        * alt[gene_start:gene_end].mean(0)
        - ref[start:end].mean(0)
    )


def get_optimal_resize_width(
    interval_width: int, config_resize: Optional[int] = None
) -> int:
    """Select an optimal resize width from AlphaGenome supported lengths."""
    from alphagenome.models.dna_client import SUPPORTED_SEQUENCE_LENGTHS

    if config_resize is not None:
        return int(config_resize)

    supported = sorted(SUPPORTED_SEQUENCE_LENGTHS.values())
    for length in supported:
        if length >= interval_width:
            return length

    max_len = supported[-1]
    print(
        f"WARNING: Interval ({interval_width:,} bp) exceeds max supported "
        f"({max_len:,} bp). Using {max_len:,} bp; predictions will be truncated."
    )
    return max_len


def setup_alphagenome(cfg: OligoConfig) -> AGContext:
    """Create model, load GTF, resize interval, fetch reference sequence,
    and run the initial wildtype prediction.
    """
    from alphagenome.data import gene_annotation, genome
    from alphagenome.models import dna_client

    if not cfg.dna_api_key or cfg.dna_api_key == "REPLACE_WITH_YOUR_KEY":
        raise ValueError(
            "dna_api_key is missing or placeholder in config. Set a valid "
            "AlphaGenome API key or use --skip-alphagenome."
        )

    model = dna_client.create(cfg.dna_api_key)
    gtf = pd.read_feather(cfg.gtf_url)
    gene_interval = gene_annotation.get_gene_interval(gtf, gene_symbol=cfg.gene_symbol)
    optimal_resize = get_optimal_resize_width(gene_interval.width, cfg.resize_width)
    interval = gene_interval.resize(optimal_resize)

    print(f"Gene: {cfg.gene_symbol}")
    print(
        f"Interval: {interval.chromosome}:{interval.start}-{interval.end} "
        f"(width={interval.width:,} bp, gene width={gene_interval.width:,} bp)"
    )

    requested_outputs = parse_output_types(cfg.requested_outputs)

    ref_seq = load_reference_sequence(
        cfg.fasta_path, interval.chromosome, interval.start, interval.end
    )

    exon_start = int(cfg.exon_intervals[0])
    exon_end = int(cfg.exon_intervals[1])
    variant_interval = genome.Interval(
        chromosome=interval.chromosome,
        start=exon_start - cfg.flank[0],
        end=exon_end + cfg.flank[1],
    )

    start_rel = variant_interval.start - interval.start
    end_rel = variant_interval.end - interval.start
    exon_start_rel = exon_start - interval.start
    exon_end_rel = exon_end - interval.start
    gene_start_rel = gene_interval.start - interval.start
    gene_end_rel = gene_interval.end - interval.start

    print(
        f"Target exon: {exon_start}-{exon_end} | "
        f"Scan region: {variant_interval.start}-{variant_interval.end} "
        f"({variant_interval.width:,} bp)"
    )

    print("Fetching wildtype prediction from AlphaGenome...")
    ref_output = model.predict_interval(
        interval=interval,
        requested_outputs=requested_outputs,
        ontology_terms=cfg.ontology_terms,
    )

    return AGContext(
        model=model,
        gene_interval=gene_interval,
        interval=interval,
        ref_seq=ref_seq,
        variant_interval=variant_interval,
        ref_output=ref_output,
        requested_outputs=requested_outputs,
        exon_start_rel=exon_start_rel,
        exon_end_rel=exon_end_rel,
        gene_start_rel=gene_start_rel,
        gene_end_rel=gene_end_rel,
        start_rel=start_rel,
        end_rel=end_rel,
    )


def _build_variant_sequence(
    ref_seq: str, position_rel: int, length: int
) -> str:
    """Return ref_seq with `length` bases at position_rel replaced by 'N'."""
    n_block = "N" * length
    return ref_seq[:position_rel] + n_block + ref_seq[position_rel + length :]


def score_asos_alphagenome(
    ctx: AGContext, cfg: OligoConfig, candidates: list[AsoCandidate]
) -> dict[str, np.ndarray]:
    """Generate masked variants for each candidate, batch-predict via
    `model.predict_sequences`, and compute the `diff_mean` efficacy score
    per requested output type.

    Returns: {output_type_name: np.ndarray of shape (len(candidates),)}
    """
    variants: list[str] = []
    for cand in candidates:
        abs_pos_rel = ctx.start_rel + cand.position
        variants.append(_build_variant_sequence(ctx.ref_seq, abs_pos_rel, cand.length))

    n = len(variants)
    print(f"AlphaGenome: predicting {n} masked variants...")
    outputs = ctx.model.predict_sequences(
        intervals=[ctx.interval] * n,
        sequences=variants,
        requested_outputs=ctx.requested_outputs,
        ontology_terms=cfg.ontology_terms,
    )

    results: dict[str, np.ndarray] = {}
    for output_type in ctx.requested_outputs:
        variant_cols: list[np.ndarray] = []
        for i, out in enumerate(outputs):
            td = filter_td(out.get(output_type), cfg)
            if td is None:
                print(
                    f"Skipping variant {i} for {output_type.name}: no tracks after filtering"
                )
                continue
            variant_cols.append(td.values.mean(axis=1))

        if not variant_cols:
            print(f"Skipping {output_type.name}: no variant tracks")
            continue

        variant_array = np.stack(variant_cols, axis=1)

        ref_td = filter_td(ctx.ref_output.get(output_type), cfg)
        if ref_td is None:
            print(f"Skipping {output_type.name}: missing reference output")
            continue
        ref_mean = ref_td.values.mean(axis=1)[:, None]

        aso_scores = diff_mean(
            ref_mean,
            variant_array,
            ctx.exon_start_rel,
            ctx.exon_end_rel,
            ctx.gene_start_rel,
            ctx.gene_end_rel,
        )
        results[output_type.name] = np.asarray(aso_scores, dtype=np.float32)

    return results
