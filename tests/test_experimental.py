"""Tests for experimental matching and aggregation."""
import numpy as np
import pandas as pd

from oligoclaude.aso_enum import AsoCandidate
from oligoclaude.experimental import (
    aggregate_experimental_candidates,
    match_scores_to_experimental,
)
from oligoclaude.sequence_utils import reverse_complement


def test_aggregate_single_candidate_per_experiment():
    cands = [
        AsoCandidate(
            aso_id="a1",
            aso_sequence_antisense="ACGT",
            genomic_target_seq="ACGT",
            position=0,
            length=4,
            measured=2.5,
            exon_label="3a",
        ),
        AsoCandidate(
            aso_id="a2",
            aso_sequence_antisense="TTTT",
            genomic_target_seq="AAAA",
            position=4,
            length=4,
            measured=1.0,
            exon_label="3a",
        ),
    ]
    scores = {"SpliceAI": np.array([0.5, -0.2]), "AlphaGenome_RNA_SEQ": np.array([0.3, 0.1])}
    df = aggregate_experimental_candidates(cands, scores)
    assert len(df) == 2
    assert df.loc[df["ASO_ID"] == "a1", "SpliceAI"].iloc[0] == 0.5
    assert df.loc[df["ASO_ID"] == "a2", "AlphaGenome_RNA_SEQ"].iloc[0] == 0.1


def test_aggregate_multiple_matches_averaged():
    cands = [
        AsoCandidate("a1", "ACGT", "ACGT", 0, 4, measured=3.0, exon_label="3a"),
        AsoCandidate("a1_m2", "ACGT", "ACGT", 20, 4, measured=3.0, exon_label="3a"),
    ]
    scores = {"SpliceAI": np.array([0.4, 0.6])}
    df = aggregate_experimental_candidates(cands, scores)
    assert len(df) == 1
    assert df["SpliceAI"].iloc[0] == 0.5
    assert df["Measured (RT-PCR)"].iloc[0] == 3.0


def test_match_scores_sliding_window():
    # Sliding-window candidates covering a 4-mer window
    cands = []
    seq = "GCTTACAGAAAG"
    for i in range(len(seq) - 4 + 1):
        target = seq[i : i + 4]
        cands.append(
            AsoCandidate(
                aso_id=f"win_{i}",
                aso_sequence_antisense=reverse_complement(target),
                genomic_target_seq=target,
                position=i,
                length=4,
            )
        )
    scores = {"SpliceAI": np.linspace(0, 1, len(cands), dtype=np.float32)}

    # Experimental ASO is the antisense of a 6-mer target; should match 3 windows.
    target_6mer = "GCTTAC"
    exp_df = pd.DataFrame(
        {
            "ASO_ID": ["e1"],
            "ASO sequence": [reverse_complement(target_6mer)],
            "Measured (RT-PCR)": [1.0],
        }
    )
    matched = match_scores_to_experimental(
        exp_df, cands, scores, pred_length=4
    )
    assert not np.isnan(matched["SpliceAI"].iloc[0])
