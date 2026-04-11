"""Tests for aso_enum."""
import pandas as pd

from oligoclaude.aso_enum import (
    enumerate_from_experimental,
    enumerate_sliding,
)
from oligoclaude.sequence_utils import reverse_complement


def test_sliding_count_step_1():
    seq = "A" * 100
    cands = enumerate_sliding(seq, start_rel=10, end_rel=50, aso_length=18, step=1)
    # (50 - 18 + 1) - 10 = 23
    assert len(cands) == 23


def test_sliding_count_step_2():
    seq = "A" * 100
    cands = enumerate_sliding(seq, start_rel=0, end_rel=20, aso_length=5, step=2)
    # positions 0, 2, 4, ..., 14, 16 (16+5=21 > 20, so stop at 15?)
    # last valid i is 20-5=15, range(0, 16, 2) = 0,2,4,6,8,10,12,14
    assert len(cands) == 8


def test_sliding_positions_are_scan_relative():
    seq = "ACGT" * 25  # 100 bp
    cands = enumerate_sliding(seq, start_rel=4, end_rel=16, aso_length=4, step=1)
    assert cands[0].position == 0
    assert cands[0].genomic_target_seq == seq[4:8]
    assert cands[-1].position == 8
    assert cands[-1].genomic_target_seq == seq[12:16]


def test_sliding_antisense_is_reverse_complement():
    seq = "ACGTACGTACGTACGTAC"
    cands = enumerate_sliding(seq, 0, 18, aso_length=6, step=1)
    for c in cands:
        assert c.aso_sequence_antisense == reverse_complement(c.genomic_target_seq)


def test_experimental_enumeration_basic():
    ref = "AAAA" + "GCTTACAGAAAGGT" + "CCCC"
    antisense = reverse_complement("GCTTACAG")  # CTGTAAGC
    exp_df = pd.DataFrame(
        {
            "ASO_ID": ["a1"],
            "ASO sequence": [antisense],
            "Measured (RT-PCR)": [2.5],
            "Region (Exon)": ["3a"],
        }
    )
    cands = enumerate_from_experimental(
        exp_df,
        ref_seq=ref,
        variant_interval_start_rel=0,
        variant_interval_end_rel=len(ref),
    )
    assert len(cands) == 1
    assert cands[0].aso_id == "a1"
    assert cands[0].position == 4
    assert cands[0].length == 8
    assert cands[0].measured == 2.5
    assert cands[0].exon_label == "3a"
    assert cands[0].genomic_target_seq == "GCTTACAG"


def test_experimental_restricted_to_scan_region():
    ref = "GCTTACAG" + "CCCC" + "GCTTACAG"  # target appears twice
    antisense = reverse_complement("GCTTACAG")
    exp_df = pd.DataFrame(
        {
            "ASO_ID": ["a1"],
            "ASO sequence": [antisense],
            "Measured (RT-PCR)": [1.0],
        }
    )
    cands = enumerate_from_experimental(
        exp_df, ref_seq=ref, variant_interval_start_rel=0, variant_interval_end_rel=8
    )
    assert len(cands) == 1
    assert cands[0].position == 0

    cands2 = enumerate_from_experimental(
        exp_df,
        ref_seq=ref,
        variant_interval_start_rel=0,
        variant_interval_end_rel=len(ref),
    )
    assert len(cands2) == 2
    assert cands2[0].aso_id == "a1"
    assert cands2[1].aso_id == "a1_m2"
