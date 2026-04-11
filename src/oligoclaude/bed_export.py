"""BED export for ASO predictions.

Ports the BED-writing logic from AlphaGenome_ASO/aso.ipynb (score_asos_and_export),
decoupled from AlphaGenome types so it works for any score source.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .aso_enum import AsoCandidate
from .config import OligoConfig


def _build_rows(
    selected: list[tuple[AsoCandidate, float]],
    chrom: str,
    strand: str,
    invert: bool,
) -> pd.DataFrame:
    """Build a BED DataFrame from a list of (candidate, score) pairs."""
    if not selected:
        return pd.DataFrame(
            columns=[
                "chrom",
                "chromStart",
                "chromEnd",
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
            ]
        )

    cands, scores = zip(*selected)
    scores_arr = np.asarray(scores, dtype=np.float32)

    norm_color = 1 + scores_arr if invert else 1 - scores_arr
    vals = np.clip(norm_color * 255, 1, 255).astype(int)
    if invert:
        itemRgb = [f"{int(v)},{int(v)},255" for v in vals]
        label_prefix = "bottom"
    else:
        itemRgb = [f"255,{int(v)},{int(v)}" for v in vals]
        label_prefix = "top"

    chrom_starts: list[int] = []
    chrom_ends: list[int] = []
    names: list[str] = []
    for i, (cand, s) in enumerate(zip(cands, scores_arr)):
        # The BED coordinates in the reference project are based on the scan
        # region; here we already receive absolute genomic coordinates via the
        # caller (which passes variant_interval_start separately).
        chrom_starts.append(cand.position)
        chrom_ends.append(cand.position + cand.length)
        names.append(f"{label_prefix}{i + 1}({s * 100:.2f}%)")

    bed_strand = "+" if strand == "-" else "-"

    df = pd.DataFrame(
        {
            "chrom": [chrom] * len(cands),
            "chromStart": chrom_starts,
            "chromEnd": chrom_ends,
            "name": names,
            "score": scores_arr,
            "strand": [bed_strand] * len(cands),
            "thickStart": chrom_starts,
            "thickEnd": chrom_ends,
            "itemRgb": itemRgb,
        }
    )
    return df


def write_bed(
    results_dir: Path,
    config_name: str,
    source: str,
    chrom: str,
    strand: str,
    candidates: list[AsoCandidate],
    scores: np.ndarray,
    variant_interval_start: int,
    samples_max: int = 20,
) -> tuple[Path, Path]:
    """Write `{config_name}_ASO_{source}.bed` and `{config_name}_ASO_{source}_full.bed`.

    Format mirrors AlphaGenome_ASO/aso.ipynb output: 9-col BED with
    `track name=... useScore=1` header. Positive scores get red itemRgb,
    negative scores get blue. Strand is inverted relative to the target
    (ASO binds the opposite strand).

    `candidates[i].position` is the 0-based offset within the scan region;
    absolute genomic coordinates are computed by adding `variant_interval_start`.
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    candidates_abs: list[AsoCandidate] = []
    for cand in candidates:
        candidates_abs.append(
            AsoCandidate(
                aso_id=cand.aso_id,
                aso_sequence_antisense=cand.aso_sequence_antisense,
                genomic_target_seq=cand.genomic_target_seq,
                position=variant_interval_start + cand.position,
                length=cand.length,
                measured=cand.measured,
                exon_label=cand.exon_label,
            )
        )

    pairs = list(zip(candidates_abs, scores.tolist()))
    pairs_pos = sorted(
        [p for p in pairs if p[1] > 0], key=lambda x: x[1], reverse=True
    )
    pairs_neg = sorted([p for p in pairs if p[1] < 0], key=lambda x: x[1])

    top = pairs_pos[:samples_max]
    bot = pairs_neg[:samples_max]

    def _safe_concat(parts: list[pd.DataFrame]) -> pd.DataFrame:
        non_empty = [p for p in parts if len(p)]
        if not non_empty:
            return _build_rows([], chrom, strand, invert=False)
        return pd.concat(non_empty, ignore_index=True)

    top_df = _build_rows(top, chrom, strand, invert=False)
    bot_df = _build_rows(bot, chrom, strand, invert=True)
    compact_df = _safe_concat([top_df, bot_df])

    full_top_df = _build_rows(pairs_pos, chrom, strand, invert=False)
    full_bot_df = _build_rows(pairs_neg, chrom, strand, invert=True)
    full_df = _safe_concat([full_top_df, full_bot_df])

    compact_path = results_dir / f"{config_name}_ASO_{source}.bed"
    full_path = results_dir / f"{config_name}_ASO_{source}_full.bed"

    with open(compact_path, "w") as f:
        f.write(
            f'track name={config_name}_ASO_{source} '
            f'description="ASO scores for {source}" visibility="pack" useScore=1\n'
        )
        compact_df.to_csv(f, sep="\t", index=False, header=False)

    with open(full_path, "w") as f:
        f.write(
            f'track name={config_name}_ASO_{source}_full '
            f'description="Full ASO scores for {source}" visibility="pack" useScore=1\n'
        )
        full_df.to_csv(f, sep="\t", index=False, header=False)

    return compact_path, full_path


def export_all(
    cfg: OligoConfig,
    chrom: str,
    variant_interval_start: int,
    candidates: list[AsoCandidate],
    all_scores: dict[str, np.ndarray],
    samples_max: int = 20,
) -> list[Path]:
    """Write a pair of BED files per score source."""
    written: list[Path] = []
    for source, arr in all_scores.items():
        if arr is None or len(arr) == 0:
            continue
        compact, full = write_bed(
            results_dir=cfg.results_dir,
            config_name=cfg.config_name,
            source=source,
            chrom=chrom,
            strand=cfg.strand,
            candidates=candidates,
            scores=np.asarray(arr, dtype=np.float32),
            variant_interval_start=variant_interval_start,
            samples_max=samples_max,
        )
        written.extend([compact, full])
    return written
