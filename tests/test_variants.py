"""Tests for `oligomcp.variants` — parser, applier, and coord map."""
from __future__ import annotations

import pytest

from oligomcp.variants import (
    AppliedVariant,
    ChromMismatch,
    OutOfWindow,
    OverlappingVariants,
    ParsedVariant,
    RefMismatch,
    VariantCoordMap,
    VariantParseError,
    apply_variants_to_ref,
    pad_or_trim_to_length,
    parse_variant,
)


# ---- parsers ---------------------------------------------------------


def test_parse_vcf_basic():
    v = parse_variant("chr5:70070740:G:A", gene_symbol="SMN2", assembly="hg38")
    assert v.chrom == "chr5"
    assert v.position == 70070740
    assert v.ref == "G"
    assert v.alt == "A"
    assert v.source == "vcf"
    assert v.variant_id == "chr5_70070740_G_A"
    assert v.notation == "chr5:70070740:G:A"


def test_parse_vcf_with_dash_delimiter_and_bare_chrom():
    v = parse_variant("5-70070740-G-A", gene_symbol="SMN2", assembly="hg38")
    assert v.chrom == "chr5"
    assert v.position == 70070740
    assert v.source == "vcf"


def test_parse_explicit_dict():
    v = parse_variant(
        {"chrom": "chr5", "position": 70070700, "ref": "AAG", "alt": ""},
        gene_symbol="SMN2", assembly="hg38",
    )
    assert v.ref == "AAG"
    assert v.alt == ""
    assert v.source == "explicit"
    assert v.variant_id == "chr5_70070700_AAG_-"


def test_parse_explicit_respects_custom_id():
    v = parse_variant(
        {"id": "pathognomonic", "chrom": "chr5", "position": 1, "ref": "A", "alt": "G"},
        gene_symbol="SMN2", assembly="hg38",
    )
    assert v.variant_id == "pathognomonic"


def test_parse_empty_string_raises():
    with pytest.raises(VariantParseError):
        parse_variant("", gene_symbol="SMN2", assembly="hg38")


def test_parse_non_dna_raises():
    with pytest.raises(VariantParseError):
        parse_variant("chr5:100:Z:G", gene_symbol="SMN2", assembly="hg38")


def test_parse_refseq_accession_normalizes_to_chrN():
    v = parse_variant("NC_000005.10:70070740:G:A", gene_symbol="SMN2", assembly="hg38")
    assert v.chrom == "chr5"


# ---- apply_variants_to_ref -------------------------------------------


def _pv(position: int, ref: str, alt: str, chrom: str = "chr5", vid: str = "v") -> ParsedVariant:
    return ParsedVariant(
        chrom=chrom, position=position, ref=ref, alt=alt,
        variant_id=vid, notation=f"{chrom}:{position}:{ref or '-'}:{alt or '-'}",
        source="explicit",
    )


def test_apply_snv_unchanged_length():
    ref_seq = "ABCDEFGHIJ"  # anchor 100 → bases at 1-based 101..110
    v = _pv(103, "C", "T")   # 1-based 103 → 0-based offset 2 → 'C'
    assert ref_seq[2] == "C"
    patient, cmap = apply_variants_to_ref(ref_seq, [v], anchor_genomic=100, chrom="chr5")
    assert patient == "ABTDEFGHIJ"
    assert cmap.total_delta() == 0
    assert cmap.ref_to_patient(102) == 2         # 0-based genomic 102 → offset 2 (the 'T')
    assert cmap.patient_to_ref(2) == 102         # round-trip


def test_apply_insertion_and_deletion_multi():
    ref_seq = "ABCDEFGHIJ"
    v1 = _pv(103, "C", "xyz", vid="sub")   # substitution expanding 1→3
    v2 = _pv(107, "G", "", vid="del")      # pure deletion
    patient, cmap = apply_variants_to_ref(
        ref_seq, [v1, v2], anchor_genomic=100, chrom="chr5",
    )
    assert patient == "ABxyzDEFHIJ"
    assert cmap.total_delta() == 2 - 1  # +2 from insertion, -1 from deletion = +1
    # Deleted base maps to None.
    assert cmap.ref_to_patient(106) is None
    # Insertion-interior offsets have no reference counterpart.
    assert cmap.patient_to_ref(3) is None
    assert cmap.patient_to_ref(4) is None
    # Bases preserved before and after the variants round-trip correctly.
    assert cmap.patient_to_ref(0) == 100
    assert cmap.patient_to_ref(7) == 105


def test_apply_rejects_wrong_chromosome():
    ref_seq = "AAAAA"
    v = _pv(101, "A", "G", chrom="chr3")
    with pytest.raises(ChromMismatch):
        apply_variants_to_ref(ref_seq, [v], anchor_genomic=100, chrom="chr5")


def test_apply_rejects_out_of_window():
    ref_seq = "AAAAA"
    v = _pv(200, "A", "G")
    with pytest.raises(OutOfWindow):
        apply_variants_to_ref(ref_seq, [v], anchor_genomic=100, chrom="chr5")


def test_apply_rejects_ref_mismatch():
    ref_seq = "ABCDEFGHIJ"
    v = _pv(103, "X", "T")     # ref at 0-based 2 is 'C', not 'X'
    with pytest.raises(RefMismatch):
        apply_variants_to_ref(ref_seq, [v], anchor_genomic=100, chrom="chr5")


def test_apply_rejects_overlap():
    ref_seq = "ABCDEFGHIJ"
    v1 = _pv(102, "BC", "xx", vid="v1")     # covers 0-based 1..2
    v2 = _pv(103, "C", "y", vid="v2")       # covers 0-based 2..2
    with pytest.raises(OverlappingVariants):
        apply_variants_to_ref(ref_seq, [v1, v2], anchor_genomic=100, chrom="chr5")


def test_apply_pure_insertion_via_anchor_convention():
    # "insertion" via VCF anchor convention: ref='C', alt='Cxy' expands C → Cxy.
    ref_seq = "ABCDEFG"
    v = _pv(103, "C", "Cxy")
    patient, cmap = apply_variants_to_ref(ref_seq, [v], anchor_genomic=100, chrom="chr5")
    assert patient == "ABCxyDEFG"
    assert cmap.patient_to_ref(2) == 102   # the anchor 'C'
    assert cmap.patient_to_ref(3) is None  # inside insertion
    assert cmap.patient_to_ref(4) is None  # inside insertion
    assert cmap.patient_to_ref(5) == 103   # 'D' preserved


# ---- pad_or_trim_to_length -------------------------------------------


def test_pad_to_length_extends_with_downstream():
    seq = "ABCDE"  # length 5

    def fetcher(chrom, start, end):
        # Simulate that genomic 105..end is 'xyz...'
        assert chrom == "chr5"
        available = "XYZ"
        return available[: end - start]

    padded = pad_or_trim_to_length(
        seq, target=8, fetcher=fetcher,
        chrom="chr5", anchor_genomic=100, original_length=5,
    )
    assert padded == "ABCDEXYZ"


def test_pad_to_length_trims_when_too_long():
    seq = "ABCDEFGHIJ"

    def fetcher(*a, **kw):
        raise AssertionError("fetcher must not be called on truncation path")

    trimmed = pad_or_trim_to_length(
        seq, target=4, fetcher=fetcher,
        chrom="chr5", anchor_genomic=100, original_length=10,
    )
    assert trimmed == "ABCD"


def test_pad_to_length_passthrough_when_equal():
    seq = "ABCDE"
    out = pad_or_trim_to_length(
        seq, target=5, fetcher=lambda *a, **k: "",
        chrom="chr5", anchor_genomic=100, original_length=5,
    )
    assert out == "ABCDE"


# ---- VariantCoordMap round-trip ---------------------------------------


def test_coord_map_no_variants_is_identity():
    cmap = VariantCoordMap(anchor_genomic=100, chrom="chr5", applied=())
    assert cmap.ref_to_patient(100) == 0
    assert cmap.ref_to_patient(105) == 5
    assert cmap.patient_to_ref(0) == 100
    assert cmap.patient_to_ref(5) == 105


def test_coord_map_snv_round_trip():
    ref_seq = "AAAAA"
    v = _pv(102, "A", "C")
    _, cmap = apply_variants_to_ref(ref_seq, [v], anchor_genomic=100, chrom="chr5")
    for g in range(100, 105):
        p = cmap.ref_to_patient(g)
        assert p is not None
        assert cmap.patient_to_ref(p) == g


def test_coord_map_deletion_round_trip():
    ref_seq = "ABCDE"
    v = _pv(102, "B", "")   # delete 'B' at 0-based 1
    patient, cmap = apply_variants_to_ref(ref_seq, [v], anchor_genomic=100, chrom="chr5")
    assert patient == "ACDE"
    # Deleted genomic position returns None; others round-trip.
    assert cmap.ref_to_patient(101) is None
    for g in (100, 102, 103, 104):
        p = cmap.ref_to_patient(g)
        assert p is not None
        assert cmap.patient_to_ref(p) == g


def test_coord_map_insertion_round_trip():
    ref_seq = "ABCD"
    v = _pv(102, "B", "Bxx")   # insert xx after 'B'
    patient, cmap = apply_variants_to_ref(ref_seq, [v], anchor_genomic=100, chrom="chr5")
    assert patient == "ABxxCD"
    # Every reference position has a patient offset; insertion interior
    # has none on the reverse mapping.
    for g in range(100, 104):
        assert cmap.patient_to_ref(cmap.ref_to_patient(g)) == g
    assert cmap.patient_to_ref(2) is None  # first inserted x
    assert cmap.patient_to_ref(3) is None  # second inserted x
