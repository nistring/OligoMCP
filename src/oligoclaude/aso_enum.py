"""ASO candidate enumeration: sliding window or experimental-driven."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

from .sequence_utils import reverse_complement


@dataclass
class AsoCandidate:
    """A single ASO candidate with genomic context.

    - `aso_sequence_antisense` is the drug sequence (reverse complement of genomic).
    - `genomic_target_seq` is the sequence the ASO binds (the genome strand).
    - `position` is 0-based offset within the scan region (variant_interval).
    - `length` is the ASO length (may vary in experimental mode).
    """

    aso_id: str
    aso_sequence_antisense: str
    genomic_target_seq: str
    position: int
    length: int
    measured: Optional[float] = None
    exon_label: Optional[str] = None


def enumerate_sliding(
    ref_seq: str, start_rel: int, end_rel: int, aso_length: int, step: int = 1
) -> list[AsoCandidate]:
    """Sliding-window ASO enumeration across [start_rel, end_rel) within ref_seq.

    For each window start position, the genomic target is ref_seq[i:i+aso_length]
    and the antisense ASO drug sequence is its reverse complement. `position` is
    stored as the offset relative to start_rel so downstream BED coordinates use
    variant_interval.start as the anchor.
    """
    candidates: list[AsoCandidate] = []
    last = end_rel - aso_length
    for i in range(start_rel, last + 1, step):
        target = ref_seq[i : i + aso_length]
        antisense = reverse_complement(target)
        candidates.append(
            AsoCandidate(
                aso_id=f"win_{i - start_rel}",
                aso_sequence_antisense=antisense,
                genomic_target_seq=target,
                position=i - start_rel,
                length=aso_length,
            )
        )
    return candidates


def enumerate_from_experimental(
    exp_df: pd.DataFrame,
    ref_seq: str,
    variant_interval_start_rel: int,
    variant_interval_end_rel: int,
    aso_seq_col: str = "ASO sequence",
    aso_id_col: str = "ASO_ID",
    measured_col: Optional[str] = "Measured (RT-PCR)",
    exon_col: Optional[str] = "Region (Exon)",
) -> list[AsoCandidate]:
    """Build candidates from an experimental CSV.

    Each experimental ASO sequence (antisense) is reverse-complemented to get
    the genomic target, then searched within the scan region of `ref_seq`.
    Every occurrence produces one AsoCandidate (rare duplicates get suffixed IDs).

    Coordinates `variant_interval_*_rel` are offsets within the `ref_seq`
    (which covers the full resized interval). Candidates are only kept if
    they fall entirely inside the scan region.
    """
    candidates: list[AsoCandidate] = []
    for _, row in exp_df.iterrows():
        antisense = str(row[aso_seq_col]).strip().upper()
        if not antisense:
            continue
        target = reverse_complement(antisense)
        length = len(target)
        exp_id = str(row[aso_id_col]) if aso_id_col in row else f"exp_{_}"
        measured = float(row[measured_col]) if measured_col and measured_col in row and pd.notna(row[measured_col]) else None
        exon_label = str(row[exon_col]) if exon_col and exon_col in row and pd.notna(row[exon_col]) else None

        search_from = variant_interval_start_rel
        match_idx = 0
        while True:
            pos = ref_seq.find(target, search_from, variant_interval_end_rel)
            if pos == -1 or pos + length > variant_interval_end_rel:
                break
            match_idx += 1
            suffix = "" if match_idx == 1 else f"_m{match_idx}"
            candidates.append(
                AsoCandidate(
                    aso_id=f"{exp_id}{suffix}",
                    aso_sequence_antisense=antisense,
                    genomic_target_seq=target,
                    position=pos - variant_interval_start_rel,
                    length=length,
                    measured=measured,
                    exon_label=exon_label,
                )
            )
            search_from = pos + 1
    return candidates
