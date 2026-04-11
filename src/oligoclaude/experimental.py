"""Experimental data loading and reverse-complement matching."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .aso_enum import AsoCandidate
from .sequence_utils import reverse_complement


_MATCH_SUFFIX_RE = re.compile(r"_m\d+$")


def _strip_match_suffix(aso_id: str) -> str:
    return _MATCH_SUFFIX_RE.sub("", aso_id)


def load_experimental(csv_path: Path) -> pd.DataFrame:
    """Load and normalize an experimental CSV file.

    Strips whitespace from column names. Expected columns:
    - `ASO_ID` (str)
    - `ASO sequence` (str, antisense drug sequence)
    - `Measured (RT-PCR)` (float)
    - `Region (Exon)` (optional str, e.g. "3a", "3b")
    """
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()
    return df


def match_scores_to_experimental(
    exp_df: pd.DataFrame,
    candidates: list[AsoCandidate],
    scores: dict[str, np.ndarray],
    pred_length: Optional[int] = None,
) -> pd.DataFrame:
    """Match predicted scores to experimental ASOs via reverse-complement lookup.

    Port of `match_aso_score` from AlphaGenome_ASO/aso_comparison.ipynb cell-2.

    For each experimental row:
    1. Reverse-complement its ASO sequence to get the genomic target.
    2. Slide a window of length `pred_length` (default: candidates[0].length)
       across the target.
    3. For each window, look up matching candidate genomic_target_seq.
    4. Take the mean of all matched scores (one per score source).

    Returns the experimental DataFrame with added columns, one per score key.
    Unmatched rows have NaN for that source.
    """
    if pred_length is None and candidates:
        pred_length = candidates[0].length

    lookups: dict[str, dict[str, float]] = {}
    for source, arr in scores.items():
        lookup: dict[str, list[float]] = {}
        for cand, val in zip(candidates, arr):
            lookup.setdefault(cand.genomic_target_seq, []).append(float(val))
        lookups[source] = {k: float(np.mean(v)) for k, v in lookup.items()}

    out = exp_df.copy()
    for source in scores.keys():
        out[source] = np.nan

    for idx, row in out.iterrows():
        antisense = str(row.get("ASO sequence", "")).strip().upper()
        if not antisense:
            continue
        target = reverse_complement(antisense)
        for source, lookup in lookups.items():
            matches: list[float] = []
            for i in range(len(target) - pred_length + 1):
                sub = target[i : i + pred_length]
                if sub in lookup:
                    matches.append(lookup[sub])
            if matches:
                out.at[idx, source] = float(np.mean(matches))

    return out


def aggregate_experimental_candidates(
    candidates: list[AsoCandidate],
    scores: dict[str, np.ndarray],
    measured_col: str = "Measured (RT-PCR)",
    aso_id_col: str = "ASO_ID",
    exon_col: str = "Region (Exon)",
) -> pd.DataFrame:
    """Aggregate candidate-level scores up to one row per experimental ASO.

    When ASOs are enumerated from experimental data, a single experimental row
    may produce multiple candidates (if its reverse-complemented target appears
    at more than one site in the scan region). This function groups candidates
    by the stripped experimental ID (e.g. `a1468` from `a1468_m2`) and averages
    the scores so each experimental row gets exactly one summary point.

    Returns a DataFrame with columns: [aso_id_col, exon_col, measured_col, *score_keys].
    """
    rows: dict[str, dict] = {}
    for i, cand in enumerate(candidates):
        base = _strip_match_suffix(cand.aso_id)
        if base not in rows:
            rows[base] = {
                aso_id_col: base,
                exon_col: cand.exon_label,
                measured_col: cand.measured,
                "ASO sequence": cand.aso_sequence_antisense,
                "_count": 0,
                **{k: [] for k in scores.keys()},
            }
        rows[base]["_count"] += 1
        for k, arr in scores.items():
            v = float(arr[i])
            if not np.isnan(v):
                rows[base][k].append(v)

    final_rows = []
    for base, r in rows.items():
        out = {
            aso_id_col: r[aso_id_col],
            exon_col: r[exon_col],
            measured_col: r[measured_col],
            "ASO sequence": r["ASO sequence"],
        }
        for k in scores.keys():
            vals = r[k]
            out[k] = float(np.mean(vals)) if vals else np.nan
        final_rows.append(out)

    return pd.DataFrame(final_rows)
