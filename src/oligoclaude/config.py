"""Config loading and validation for OligoClaude."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


_DEFAULT_GTF_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/hg38/"
    "gencode.v46.annotation.gtf.gz.feather"
)


@dataclass
class OligoConfig:
    gene_symbol: str
    exon_intervals: tuple[int, int]
    strand: str
    assembly: str
    gtf_url: str
    results_dir: Path
    data_dir: Path
    ontology_terms: list[str]
    track_filter: str
    requested_outputs: list[str]
    ASO_length: int
    flank: tuple[int, int]

    fasta_path: Optional[Path] = None
    dna_api_key: Optional[str] = None
    aso_step: int = 1
    experimental_data: Optional[Path] = None
    target_mode: str = "exclude"
    spliceai_batch: int = 12
    spliceai_threads: Optional[int] = None
    resize_width: Optional[int] = None
    config_name: str = ""


def load_config(path: Path) -> OligoConfig:
    """Load a JSON config file into an OligoConfig dataclass.

    Paths in the JSON are resolved relative to the config file's directory
    so configs are portable.
    """
    path = Path(path).expanduser().resolve()
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    base_dir = path.parent

    def _resolve(p: str | None) -> Optional[Path]:
        if p is None:
            return None
        pp = Path(p).expanduser()
        return pp if pp.is_absolute() else (base_dir / pp).resolve()

    if raw.get("target_mode", "exclude") not in ("exclude", "include"):
        raise ValueError(
            f"target_mode must be 'exclude' or 'include', got {raw.get('target_mode')!r}"
        )

    exon = raw["exon_intervals"]
    flank = raw.get("flank", [200, 200])

    cfg = OligoConfig(
        gene_symbol=raw["gene_symbol"],
        exon_intervals=(int(exon[0]), int(exon[1])),
        strand=raw.get("strand", "+"),
        assembly=raw.get("assembly", "hg38"),
        fasta_path=_resolve(raw.get("fasta_path")),
        gtf_url=raw.get("gtf_url") or _DEFAULT_GTF_URL,
        results_dir=_resolve(raw.get("results_dir", "results")),
        data_dir=_resolve(raw.get("data_dir", "data")),
        dna_api_key=raw.get("dna_api_key") or None,
        ontology_terms=list(raw.get("ontology_terms", [])),
        track_filter=raw.get("track_filter", ""),
        requested_outputs=list(raw.get("requested_outputs", ["RNA_SEQ", "SPLICE_SITE_USAGE"])),
        ASO_length=int(raw.get("ASO_length", 18)),
        flank=(int(flank[0]), int(flank[1])),
        aso_step=int(raw.get("aso_step", 1)),
        experimental_data=_resolve(raw.get("experimental_data")),
        target_mode=raw.get("target_mode", "exclude"),
        spliceai_batch=int(raw.get("spliceai_batch", 12)),
        spliceai_threads=raw.get("spliceai_threads"),
        resize_width=raw.get("resize_width"),
        config_name=path.stem,
    )
    return cfg
