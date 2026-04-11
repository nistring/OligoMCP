"""Tests for bed_export."""
from pathlib import Path

import numpy as np

from oligoclaude.aso_enum import AsoCandidate
from oligoclaude.bed_export import write_bed


def _mk_candidates(n: int, length: int = 18) -> list:
    cands = []
    for i in range(n):
        cands.append(
            AsoCandidate(
                aso_id=f"win_{i}",
                aso_sequence_antisense="N" * length,
                genomic_target_seq="N" * length,
                position=i,
                length=length,
            )
        )
    return cands


def test_write_bed_header_and_columns(tmp_path: Path):
    cands = _mk_candidates(5)
    scores = np.array([0.5, -0.3, 0.1, -0.2, 0.4], dtype=np.float32)

    compact, full = write_bed(
        results_dir=tmp_path,
        config_name="TEST",
        source="RNA_SEQ",
        chrom="chr1",
        strand="+",
        candidates=cands,
        scores=scores,
        variant_interval_start=1000,
        samples_max=20,
    )

    assert compact.name == "TEST_ASO_RNA_SEQ.bed"
    assert full.name == "TEST_ASO_RNA_SEQ_full.bed"

    lines = compact.read_text().strip().split("\n")
    assert lines[0].startswith("track name=TEST_ASO_RNA_SEQ ")
    assert "useScore=1" in lines[0]
    assert 'visibility="pack"' in lines[0]

    # 5 rows (3 positive + 2 negative)
    data_lines = lines[1:]
    assert len(data_lines) == 5


def test_write_bed_strand_inversion_plus(tmp_path: Path):
    cands = _mk_candidates(1)
    scores = np.array([0.5], dtype=np.float32)
    compact, _ = write_bed(
        results_dir=tmp_path,
        config_name="T",
        source="X",
        chrom="chr1",
        strand="+",
        candidates=cands,
        scores=scores,
        variant_interval_start=0,
    )
    data_line = compact.read_text().strip().split("\n")[1]
    cols = data_line.split("\t")
    assert cols[5] == "-"


def test_write_bed_strand_inversion_minus(tmp_path: Path):
    cands = _mk_candidates(1)
    scores = np.array([0.5], dtype=np.float32)
    compact, _ = write_bed(
        results_dir=tmp_path,
        config_name="T",
        source="X",
        chrom="chr1",
        strand="-",
        candidates=cands,
        scores=scores,
        variant_interval_start=0,
    )
    data_line = compact.read_text().strip().split("\n")[1]
    cols = data_line.split("\t")
    assert cols[5] == "+"


def test_write_bed_color_positive_is_red(tmp_path: Path):
    cands = _mk_candidates(1)
    scores = np.array([0.5], dtype=np.float32)
    compact, _ = write_bed(
        results_dir=tmp_path,
        config_name="T",
        source="X",
        chrom="chr1",
        strand="+",
        candidates=cands,
        scores=scores,
        variant_interval_start=0,
    )
    cols = compact.read_text().strip().split("\n")[1].split("\t")
    rgb = cols[8].split(",")
    assert rgb[0] == "255"


def test_write_bed_color_negative_is_blue(tmp_path: Path):
    cands = _mk_candidates(1)
    scores = np.array([-0.5], dtype=np.float32)
    compact, _ = write_bed(
        results_dir=tmp_path,
        config_name="T",
        source="X",
        chrom="chr1",
        strand="+",
        candidates=cands,
        scores=scores,
        variant_interval_start=0,
    )
    cols = compact.read_text().strip().split("\n")[1].split("\t")
    rgb = cols[8].split(",")
    assert rgb[2] == "255"


def test_write_bed_absolute_coordinates(tmp_path: Path):
    cands = _mk_candidates(1, length=18)
    scores = np.array([0.3], dtype=np.float32)
    compact, _ = write_bed(
        results_dir=tmp_path,
        config_name="T",
        source="X",
        chrom="chr3",
        strand="+",
        candidates=cands,
        scores=scores,
        variant_interval_start=9429626,
    )
    cols = compact.read_text().strip().split("\n")[1].split("\t")
    assert cols[0] == "chr3"
    assert int(cols[1]) == 9429626
    assert int(cols[2]) == 9429644
    assert int(cols[6]) == 9429626
    assert int(cols[7]) == 9429644
