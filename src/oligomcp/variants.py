"""Variant parsing + application: turn user-supplied mutations into a patient
baseline for ASO design.

Supported input forms (see `parse_variant`):

  - VCF-style:            ``chr5:70070740:C:T``  or  ``5-70070740-C-T``
  - HGVS genomic (g.):    ``chr5:g.70070740C>T``, ``NC_000005.10:g.70070740_70070742del``
  - HGVS coding (c.):     ``c.840C>T``   (gene_symbol + strand resolved via mygene)
  - rsID:                 ``rs1800112``  (dbSNP esummary; disk-cached)
  - ClinVar accession:    ``VCV000000001``  (ClinVar esummary; disk-cached)
  - Explicit dict:        ``{"chrom":"chr5","position":70070740,"ref":"C","alt":"T"}``

Every parser normalizes to the same `ParsedVariant` record — chrom as
``chrN``, position 1-based (VCF convention), ref/alt uppercase (``""`` for
pure indels) and always on the positive genomic strand regardless of the
source notation.

`apply_variants_to_ref` edits variants into ``ref_seq`` and returns a
`VariantCoordMap` for later coord translation:

  - ``ref_to_patient(gpos)`` — 0-based genomic position → offset in
    ``patient_seq``. Returns ``None`` if the base was deleted.
  - ``patient_to_ref(poff)`` — offset in ``patient_seq`` → 0-based
    genomic position. Returns ``None`` for offsets inside an insertion
    (no reference counterpart exists).

`pad_or_trim_to_length` enforces the fixed-width constraints of the
downstream predictors (AlphaGenome ``interval.width``, SpliceAI's
15 000 bp window) by extending from the downstream reference or
truncating the right edge after indels.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, replace
from typing import Callable, Iterable, Optional

from .core import COMP, reverse_complement


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class VariantError(RuntimeError):
    """Base class for all variant-related errors surfaced to the user."""


class VariantParseError(VariantError):
    """Raised when a notation cannot be parsed into a ParsedVariant."""


class VariantApplicationError(VariantError):
    """Raised while editing variants into ref_seq (ref mismatch, overlap, ...)."""


class ChromMismatch(VariantApplicationError):
    """The variant's chromosome doesn't match the loaded ref_seq window."""


class OutOfWindow(VariantApplicationError):
    """The variant's position falls outside the loaded ref_seq window."""


class RefMismatch(VariantApplicationError):
    """ref_seq bases at the variant position don't match variant.ref."""


class OverlappingVariants(VariantApplicationError):
    """Two variants cover overlapping positions in ref_seq."""


class IndelTooLarge(VariantApplicationError):
    """An applied indel would grow patient_seq past the available window."""


class ExonDeletedByVariant(VariantApplicationError):
    """The applied variants remove (part of) the target exon, leaving no
    design surface. The ASO pipeline has nothing to score in that case."""


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ParsedVariant:
    chrom: str          # normalized "chr5"
    position: int       # 1-based genomic (VCF convention)
    ref: str            # uppercase DNA on the + strand; "" for pure insertion
    alt: str            # uppercase DNA on the + strand; "" for pure deletion
    variant_id: str     # user-supplied or auto "chr5_70070740_C_T"
    notation: str       # original input string (or repr(dict)) for provenance
    source: str         # "vcf" | "hgvs_g" | "hgvs_c" | "rsid" | "clinvar" | "explicit"


@dataclass(frozen=True)
class AppliedVariant:
    variant: ParsedVariant
    ref_offset: int    # 0-based offset in the ORIGINAL ref_seq
    delta: int         # len(alt) - len(ref). +N = insertion, -N = deletion.


@dataclass(frozen=True)
class VariantCoordMap:
    """Maps between reference-genome coordinates and patient_seq offsets.

    Both inputs and outputs use 0-based genomic positions. Position
    conventions match `cfg.exon_intervals` (0-based) and `ag_ctx.interval`
    (half-open 0-based), so the same coord map plugs straight into the
    existing workflow wiring.
    """

    anchor_genomic: int
    chrom: str
    applied: tuple[AppliedVariant, ...]   # sorted by ref_offset ascending

    def ref_to_patient(self, genomic_pos: int) -> Optional[int]:
        """Return the 0-based patient_seq offset of `genomic_pos`, or None
        if the base at that genomic position was deleted by a variant.

        For substitutions / delins where the queried position lies inside
        the edited region, the return value points into the alt string.
        """
        ref_offset = genomic_pos - self.anchor_genomic
        cursor_patient = 0
        cursor_ref = 0
        for av in self.applied:
            seg_len = av.ref_offset - cursor_ref
            if ref_offset < cursor_ref + seg_len:
                return cursor_patient + (ref_offset - cursor_ref)
            cursor_patient += seg_len
            cursor_ref += seg_len
            ref_len = len(av.variant.ref)
            alt_len = len(av.variant.alt)
            if ref_offset < cursor_ref + ref_len:
                offset_in_ref = ref_offset - cursor_ref
                if offset_in_ref < alt_len:
                    return cursor_patient + offset_in_ref
                return None  # deleted by this variant
            cursor_patient += alt_len
            cursor_ref += ref_len
        return cursor_patient + (ref_offset - cursor_ref)

    def patient_to_ref(self, patient_offset: int) -> Optional[int]:
        """Return the 0-based genomic position corresponding to
        `patient_offset`, or None if that offset falls inside an insertion
        (no reference base exists for the extra alt content).
        """
        cursor_patient = 0
        cursor_ref = 0
        for av in self.applied:
            seg_len = av.ref_offset - cursor_ref
            if patient_offset < cursor_patient + seg_len:
                return self.anchor_genomic + cursor_ref + (patient_offset - cursor_patient)
            cursor_patient += seg_len
            cursor_ref += seg_len
            ref_len = len(av.variant.ref)
            alt_len = len(av.variant.alt)
            if patient_offset < cursor_patient + alt_len:
                offset_in_alt = patient_offset - cursor_patient
                if offset_in_alt < ref_len:
                    return self.anchor_genomic + cursor_ref + offset_in_alt
                return None  # inside an insertion
            cursor_patient += alt_len
            cursor_ref += ref_len
        return self.anchor_genomic + cursor_ref + (patient_offset - cursor_patient)

    def total_delta(self) -> int:
        return sum(av.delta for av in self.applied)


# ---------------------------------------------------------------------------
# Chromosome normalization
# ---------------------------------------------------------------------------


# RefSeq genomic accessions → UCSC-style chromosome names. 24 human
# chromosomes + mitochondrion. Hardcoded so the parser works offline.
_REFSEQ_TO_CHROM = {
    "NC_000001": "chr1", "NC_000002": "chr2", "NC_000003": "chr3",
    "NC_000004": "chr4", "NC_000005": "chr5", "NC_000006": "chr6",
    "NC_000007": "chr7", "NC_000008": "chr8", "NC_000009": "chr9",
    "NC_000010": "chr10", "NC_000011": "chr11", "NC_000012": "chr12",
    "NC_000013": "chr13", "NC_000014": "chr14", "NC_000015": "chr15",
    "NC_000016": "chr16", "NC_000017": "chr17", "NC_000018": "chr18",
    "NC_000019": "chr19", "NC_000020": "chr20", "NC_000021": "chr21",
    "NC_000022": "chr22", "NC_000023": "chrX", "NC_000024": "chrY",
    "NC_012920": "chrM",
}


def _normalize_chrom(raw: str) -> str:
    s = str(raw).strip()
    if not s:
        raise VariantParseError("Empty chromosome.")
    base = s.split(".", 1)[0]
    if base in _REFSEQ_TO_CHROM:
        return _REFSEQ_TO_CHROM[base]
    if s.lower().startswith("chr"):
        return "chr" + s[3:]
    return f"chr{s}"


def _uppercase_bases(s: Optional[str]) -> str:
    if s is None:
        return ""
    s = s.strip().upper()
    if s in {"-", ".", "DEL"}:
        return ""
    if not re.fullmatch(r"[ACGTN]*", s):
        raise VariantParseError(f"Non-DNA characters in base string {s!r}")
    return s


def _default_variant_id(chrom: str, position: int, ref: str, alt: str) -> str:
    return f"{chrom}_{position}_{ref or '-'}_{alt or '-'}"


# ---------------------------------------------------------------------------
# VCF-style parser
# ---------------------------------------------------------------------------


_VCF_RE = re.compile(
    r"^(?P<chrom>[\w.]+)[:\-](?P<pos>\d+)[:\-]"
    r"(?P<ref>[A-Za-z.\-]*)[:\-](?P<alt>[A-Za-z.\-]*)$"
)


def _parse_vcf(s: str, *, raw_notation: str) -> ParsedVariant:
    m = _VCF_RE.match(s)
    if not m:
        raise VariantParseError(f"Not VCF-shaped: {s!r}")
    chrom = _normalize_chrom(m.group("chrom"))
    position = int(m.group("pos"))
    ref = _uppercase_bases(m.group("ref"))
    alt = _uppercase_bases(m.group("alt"))
    return ParsedVariant(
        chrom=chrom, position=position, ref=ref, alt=alt,
        variant_id=_default_variant_id(chrom, position, ref, alt),
        notation=raw_notation, source="vcf",
    )


# ---------------------------------------------------------------------------
# HGVS-genomic parser
# ---------------------------------------------------------------------------


# Accepts an optional `<chrom>:` or `<NC accession>:` prefix, then `g.`,
# then one of the five operator shapes below. Reference bases embedded in
# `del`/`dup`/`delins` are optional — we fetch from UCSC when absent so
# `apply_variants_to_ref`'s ref-check remains meaningful.
_HGVS_G_RE = re.compile(
    r"^(?:(?P<chrom>[\w.]+):)?g\."
    r"(?P<start>\d+)(?:_(?P<end>\d+))?"
    r"(?P<op>"
        r"(?P<sub_ref>[ACGT])>(?P<sub_alt>[ACGT])"
        r"|delins(?P<delins_alt>[ACGT]+)"
        r"|del(?P<del_bases>[ACGT]*)"
        r"|ins(?P<ins_alt>[ACGT]+)"
        r"|dup(?P<dup_bases>[ACGT]*)"
    r")$",
    re.IGNORECASE,
)


FetchRef = Callable[[str, int, int], str]   # (chrom, start, end) -> DNA


def _parse_hgvs_g(
    s: str,
    *,
    raw_notation: str,
    fetch_ref: Optional[FetchRef],
    default_chrom: Optional[str] = None,
    source_tag: str = "hgvs_g",
) -> ParsedVariant:
    m = _HGVS_G_RE.match(s.strip())
    if not m:
        raise VariantParseError(f"Not HGVS-genomic: {s!r}")
    chrom_raw = m.group("chrom") or default_chrom
    if not chrom_raw:
        raise VariantParseError(
            f"HGVS-genomic {s!r} has no chromosome prefix and no default was available."
        )
    chrom = _normalize_chrom(chrom_raw)
    start = int(m.group("start"))
    end = int(m.group("end")) if m.group("end") else start

    def _fetch(a: int, b: int) -> str:
        if fetch_ref is None:
            raise VariantParseError(
                f"HGVS notation {s!r} needs reference bases but no fetcher was "
                "provided (parser was invoked offline)."
            )
        return fetch_ref(chrom, a, b).upper()

    if m.group("sub_ref"):
        ref = m.group("sub_ref").upper()
        alt = m.group("sub_alt").upper()
        position = start
    elif m.group("delins_alt"):
        ref = _fetch(start - 1, end)
        alt = m.group("delins_alt").upper()
        position = start
    elif m.group("op").lower().startswith("delins"):  # guarded above, but keep the branch explicit
        raise VariantParseError(f"delins branch not taken for {s!r}")
    elif m.group("op").lower().startswith("del"):
        given = (m.group("del_bases") or "").upper()
        ref = given if given else _fetch(start - 1, end)
        alt = ""
        position = start
    elif m.group("ins_alt"):
        # HGVS insertion: c/g.N_N+1insXXX means insert XXX between bases
        # N and N+1. Re-encode as a VCF-anchored substitution so the
        # applier doesn't need an insertion-specific path: ref=base[N],
        # alt=base[N]+XXX, position=N.
        anchor_base = _fetch(start - 1, start)
        ref = anchor_base
        alt = anchor_base + m.group("ins_alt").upper()
        position = start
    elif m.group("op").lower().startswith("dup"):
        given = (m.group("dup_bases") or "").upper()
        dup_seq = given if given else _fetch(start - 1, end)
        anchor = _fetch(end - 1, end)
        ref = anchor
        alt = anchor + dup_seq
        position = end
    else:
        raise VariantParseError(f"Unrecognized HGVS-genomic operator in {s!r}")

    return ParsedVariant(
        chrom=chrom, position=position, ref=ref, alt=alt,
        variant_id=_default_variant_id(chrom, position, ref, alt),
        notation=raw_notation, source=source_tag,
    )


# ---------------------------------------------------------------------------
# HGVS-coding parser
# ---------------------------------------------------------------------------


# Positive integer positions only; intron offsets (+/-) and UTR (*) are
# rejected with a clear message. Everything else mirrors the g. shapes.
_HGVS_C_RE = re.compile(
    r"^c\.(?P<start>-?\*?\d+(?:[+\-]\d+)?)(?:_(?P<end>-?\*?\d+(?:[+\-]\d+)?))?"
    r"(?P<op>"
        r"(?P<sub_ref>[ACGT])>(?P<sub_alt>[ACGT])"
        r"|delins(?P<delins_alt>[ACGT]+)"
        r"|del(?P<del_bases>[ACGT]*)"
        r"|ins(?P<ins_alt>[ACGT]+)"
        r"|dup(?P<dup_bases>[ACGT]*)"
    r")$",
    re.IGNORECASE,
)


def _cdot_to_genomic(
    cdot: int,
    *,
    exons: list[list[int]],
    cdsstart: int,
    cdsend: int,
    strand: str,
) -> int:
    """Map c.N (1-based, positive integer) to a 0-based genomic position.

    Exons are ``[start, end)`` half-open in ascending genomic order (the
    convention mygene.info returns). ``cdsstart`` / ``cdsend`` define the
    half-open CDS interval; c.1 is the first CDS base in transcription
    order (genomic ``cdsstart`` on ``+`` strand, genomic ``cdsend-1`` on
    ``-`` strand).
    """
    if cdot <= 0:
        raise VariantParseError(
            f"c.{cdot} is a UTR/negative coordinate; v1 supports positive integer c.N only."
        )

    # Restrict each exon to its CDS intersection, in transcription order.
    cds_segments: list[tuple[int, int]] = []
    for a, b in exons:
        seg_start = max(a, cdsstart)
        seg_end = min(b, cdsend)
        if seg_end > seg_start:
            cds_segments.append((seg_start, seg_end))
    if strand == "-":
        cds_segments.sort(key=lambda s: -s[0])
    else:
        cds_segments.sort(key=lambda s: s[0])

    remaining = cdot
    for seg_start, seg_end in cds_segments:
        length = seg_end - seg_start
        if remaining <= length:
            if strand == "+":
                return seg_start + (remaining - 1)
            return seg_end - remaining
        remaining -= length

    raise VariantParseError(f"c.{cdot} is beyond the coding region.")


def _parse_hgvs_c(
    s: str,
    *,
    raw_notation: str,
    gene_symbol: str,
    assembly: str,
) -> ParsedVariant:
    m = _HGVS_C_RE.match(s.strip())
    if not m:
        raise VariantParseError(f"Not HGVS-coding: {s!r}")
    start_raw = m.group("start")
    end_raw = m.group("end")
    if any(ch in start_raw for ch in "*+-") or (end_raw and any(ch in end_raw for ch in "*+-")):
        raise VariantParseError(
            f"c. notation with UTR (*) or intron offsets (+/-) is not supported "
            f"in v1: {s!r}. Use genomic form (chrN:g...) instead."
        )

    from .resources import canonical_transcript_exons, lookup_gene_info

    info = lookup_gene_info(gene_symbol, assembly)
    tx_id, exons, cdsstart, cdsend = canonical_transcript_exons(info)
    if cdsstart is None or cdsend is None:
        raise VariantParseError(
            f"Canonical transcript {tx_id} for {gene_symbol!r} has no CDS "
            "recorded in mygene.info; cannot resolve c. coordinates."
        )
    strand = info["strand"]
    start_cdot = int(start_raw)
    end_cdot = int(end_raw) if end_raw else start_cdot

    gstart = _cdot_to_genomic(start_cdot, exons=exons, cdsstart=cdsstart, cdsend=cdsend, strand=strand)
    gend = _cdot_to_genomic(end_cdot, exons=exons, cdsstart=cdsstart, cdsend=cdsend, strand=strand)
    gstart, gend = min(gstart, gend), max(gstart, gend)

    def _maybe_rc(b: str) -> str:
        return reverse_complement(b) if strand == "-" else b

    chrom = info["chrom"]
    fetcher = _default_fasta_fetcher(assembly)

    if m.group("sub_ref"):
        ref = m.group("sub_ref").upper().translate(COMP) if strand == "-" else m.group("sub_ref").upper()
        alt = m.group("sub_alt").upper().translate(COMP) if strand == "-" else m.group("sub_alt").upper()
        return ParsedVariant(
            chrom=chrom, position=gstart + 1, ref=ref, alt=alt,
            variant_id=_default_variant_id(chrom, gstart + 1, ref, alt),
            notation=raw_notation, source="hgvs_c",
        )

    # Indels: build an equivalent g. string and dispatch to the g. parser
    # so we don't duplicate the delins/ins/dup/del branches here. The only
    # subtlety is base complementation on - strand genes.
    if m.group("delins_alt"):
        alt = _maybe_rc(m.group("delins_alt").upper())
        synthetic = f"{chrom}:g.{gstart + 1}_{gend + 1}delins{alt}"
    elif m.group("op").lower().startswith("del"):
        given = (m.group("del_bases") or "").upper()
        alt_bases = _maybe_rc(given) if given else ""
        synthetic = (
            f"{chrom}:g.{gstart + 1}_{gend + 1}del{alt_bases}"
            if alt_bases else f"{chrom}:g.{gstart + 1}_{gend + 1}del"
        )
    elif m.group("ins_alt"):
        ins = _maybe_rc(m.group("ins_alt").upper())
        synthetic = f"{chrom}:g.{gstart + 1}_{gend + 1}ins{ins}"
    elif m.group("op").lower().startswith("dup"):
        given = (m.group("dup_bases") or "").upper()
        dup = _maybe_rc(given) if given else ""
        synthetic = (
            f"{chrom}:g.{gstart + 1}_{gend + 1}dup{dup}"
            if dup else f"{chrom}:g.{gstart + 1}_{gend + 1}dup"
        )
    else:
        raise VariantParseError(f"Unrecognized c. operator in {s!r}")

    intermediate = _parse_hgvs_g(
        synthetic,
        raw_notation=raw_notation,
        fetch_ref=fetcher,
        default_chrom=chrom,
        source_tag="hgvs_c",
    )
    return replace(intermediate, variant_id=_default_variant_id(
        intermediate.chrom, intermediate.position, intermediate.ref, intermediate.alt,
    ))


# ---------------------------------------------------------------------------
# rsID / ClinVar parsers
# ---------------------------------------------------------------------------


def _parse_rsid(s: str, *, raw_notation: str, assembly: str) -> ParsedVariant:
    from .resources import lookup_rsid_variant

    chrom, position, ref, alt = lookup_rsid_variant(s, assembly=assembly)
    chrom_n = _normalize_chrom(chrom)
    return ParsedVariant(
        chrom=chrom_n, position=position, ref=ref, alt=alt,
        variant_id=_default_variant_id(chrom_n, position, ref, alt),
        notation=raw_notation, source="rsid",
    )


def _parse_clinvar(s: str, *, raw_notation: str, assembly: str) -> ParsedVariant:
    from .resources import lookup_clinvar_variant

    chrom, position, ref, alt = lookup_clinvar_variant(s, assembly=assembly)
    chrom_n = _normalize_chrom(chrom)
    return ParsedVariant(
        chrom=chrom_n, position=position, ref=ref, alt=alt,
        variant_id=_default_variant_id(chrom_n, position, ref, alt),
        notation=raw_notation, source="clinvar",
    )


# ---------------------------------------------------------------------------
# Explicit dict parser
# ---------------------------------------------------------------------------


def _parse_explicit(d: dict) -> ParsedVariant:
    try:
        chrom = _normalize_chrom(d["chrom"])
        position = int(d["position"])
    except KeyError as e:
        raise VariantParseError(
            f"Explicit variant dict missing required field {e}: {d!r}"
        ) from e
    ref = _uppercase_bases(d.get("ref", ""))
    alt = _uppercase_bases(d.get("alt", ""))
    variant_id = str(
        d.get("id") or d.get("variant_id") or _default_variant_id(chrom, position, ref, alt)
    )
    return ParsedVariant(
        chrom=chrom, position=position, ref=ref, alt=alt,
        variant_id=variant_id, notation=json.dumps(d, sort_keys=True),
        source="explicit",
    )


# ---------------------------------------------------------------------------
# Fasta fetcher (UCSC-backed)
# ---------------------------------------------------------------------------


def _default_fasta_fetcher(assembly: str) -> FetchRef:
    from .resources import fetch_sequence_ucsc

    def _fetch(chrom: str, start: int, end: int) -> str:
        return fetch_sequence_ucsc(assembly, chrom, start, end)

    return _fetch


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


def parse_variant(
    item: object,
    *,
    gene_symbol: str,
    assembly: str,
) -> ParsedVariant:
    """Normalize any supported input form into a `ParsedVariant`.

    `gene_symbol` is only consulted for HGVS-coding (c.) notations;
    `assembly` is required for rsID / ClinVar lookups and for UCSC ref
    fetches used by some HGVS-g operators.

    Raises `VariantParseError` on unrecognized input.
    """
    if isinstance(item, dict):
        return _parse_explicit(item)
    if not isinstance(item, str):
        raise VariantParseError(
            f"Variant entry must be str or dict, got {type(item).__name__}: {item!r}"
        )

    s = item.strip()
    if not s:
        raise VariantParseError("Empty variant string.")

    if re.match(r"^rs\d+$", s, flags=re.IGNORECASE):
        return _parse_rsid(s, raw_notation=item, assembly=assembly)
    if re.match(r"^VCV\d+$", s, flags=re.IGNORECASE):
        return _parse_clinvar(s, raw_notation=item, assembly=assembly)
    if re.match(r"^c\.", s, flags=re.IGNORECASE):
        return _parse_hgvs_c(s, raw_notation=item, gene_symbol=gene_symbol, assembly=assembly)
    if re.search(r"(^|:)g\.", s, flags=re.IGNORECASE):
        return _parse_hgvs_g(
            s, raw_notation=item,
            fetch_ref=_default_fasta_fetcher(assembly),
        )
    return _parse_vcf(s, raw_notation=item)


# ---------------------------------------------------------------------------
# Sequence editing
# ---------------------------------------------------------------------------


def apply_variants_to_ref(
    ref_seq: str,
    variants: Iterable[ParsedVariant],
    *,
    anchor_genomic: int,
    chrom: str,
) -> tuple[str, VariantCoordMap]:
    """Edit each variant into ref_seq. Indels are allowed.

    Returns ``(patient_seq, coord_map)``. Raises `VariantApplicationError`
    subclasses on any validation failure (chrom mismatch, out-of-window,
    ref mismatch, overlap). Variants must not overlap; application order
    is deterministic (sorted by genomic position).
    """
    vlist = list(variants)
    applied: list[AppliedVariant] = []
    for v in vlist:
        if v.chrom != chrom:
            raise ChromMismatch(
                f"Variant {v.variant_id} is on {v.chrom} but ref_seq is on {chrom}."
            )
        ref_offset = v.position - 1 - anchor_genomic
        if ref_offset < 0 or ref_offset + len(v.ref) > len(ref_seq):
            raise OutOfWindow(
                f"Variant {v.variant_id} at {v.chrom}:{v.position} ({v.ref or '-'}>{v.alt or '-'}) "
                f"falls outside the loaded window "
                f"{chrom}:{anchor_genomic}-{anchor_genomic + len(ref_seq)}."
            )
        if v.ref:
            seen = ref_seq[ref_offset : ref_offset + len(v.ref)].upper()
            if seen != v.ref:
                raise RefMismatch(
                    f"Ref mismatch at {v.chrom}:{v.position}: variant says ref={v.ref!r}, "
                    f"ref_seq has {seen!r}. Wrong assembly, wrong strand, or a "
                    "position off by one."
                )
        applied.append(AppliedVariant(
            variant=v, ref_offset=ref_offset, delta=len(v.alt) - len(v.ref),
        ))

    applied.sort(key=lambda av: av.ref_offset)

    for a, b in zip(applied, applied[1:]):
        a_end = a.ref_offset + len(a.variant.ref)
        if b.ref_offset < a_end:
            raise OverlappingVariants(
                f"Variants {a.variant.variant_id} and {b.variant.variant_id} overlap "
                f"at ref_offset {b.ref_offset} (< {a_end})."
            )
        if a.ref_offset == b.ref_offset and (not a.variant.ref) and (not b.variant.ref):
            raise OverlappingVariants(
                f"Two pure insertions share anchor {a.ref_offset}: "
                f"{a.variant.variant_id} and {b.variant.variant_id}."
            )

    parts: list[str] = []
    cursor = 0
    for av in applied:
        if av.ref_offset > cursor:
            parts.append(ref_seq[cursor:av.ref_offset])
        parts.append(av.variant.alt)
        cursor = av.ref_offset + len(av.variant.ref)
    if cursor < len(ref_seq):
        parts.append(ref_seq[cursor:])
    patient_seq = "".join(parts)

    coord_map = VariantCoordMap(
        anchor_genomic=anchor_genomic,
        chrom=chrom,
        applied=tuple(applied),
    )
    return patient_seq, coord_map


# ---------------------------------------------------------------------------
# Fixed-length window shaping
# ---------------------------------------------------------------------------


def pad_or_trim_to_length(
    seq: str,
    *,
    target: int,
    fetcher: FetchRef,
    chrom: str,
    anchor_genomic: int,
    original_length: int,
) -> str:
    """Pad with downstream reference or trim the right edge to match ``target``.

    ``anchor_genomic + original_length`` is the right edge of the
    original pre-edit window — where padding is drawn from. For
    deletions (len(seq) < target) we fetch downstream reference bases to
    re-fill. For insertions (len(seq) > target) we truncate. Returned
    string always has exactly ``target`` characters.
    """
    current = len(seq)
    if current == target:
        return seq
    if current > target:
        return seq[:target]
    need = target - current
    pad_start = anchor_genomic + original_length
    pad_end = pad_start + need
    padding = fetcher(chrom, pad_start, pad_end)
    if len(padding) < need:
        raise IndelTooLarge(
            f"Not enough downstream reference to pad {current} → {target} bp "
            f"(need {need} bp from {chrom}:{pad_start}-{pad_end}, got {len(padding)}). "
            "Widen cfg.resize_width or split the variant set across multiple runs."
        )
    return seq + padding[:need].upper()


# ---------------------------------------------------------------------------
# Serialization helpers for run metadata
# ---------------------------------------------------------------------------


def applied_variants_to_records(coord_map: VariantCoordMap) -> list[dict]:
    """Flatten a VariantCoordMap into a JSON-serializable list of dicts."""
    return [
        {
            "variant_id": av.variant.variant_id,
            "notation": av.variant.notation,
            "source": av.variant.source,
            "chrom": av.variant.chrom,
            "position": av.variant.position,
            "ref": av.variant.ref,
            "alt": av.variant.alt,
            "ref_offset": av.ref_offset,
            "delta": av.delta,
        }
        for av in coord_map.applied
    ]
