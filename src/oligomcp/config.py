"""Config loading and validation for OligoMCP."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


_DEFAULT_GTF_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/hg38/"
    "gencode.v46.annotation.gtf.gz.feather"
)

# Fields that materially change the ASO design output and whose defaults the
# user should consciously accept. Missing keys in a config JSON are surfaced
# as a `needs_info` status so Claude can ask the user before silently using
# the listed default. `strand`/`assembly` are excluded because they're almost
# always derivable (mygene) or the hg38 default is right.
_OPINIONATED_FIELDS: list[tuple[str, object, str]] = [
    ("ASO_length", 18, "ASO length in bp (typical: 15-22)."),
    ("aso_step", 1, "Sliding-window step; 1 = score every position, 5 = every 5th."),
    ("flank", [200, 200], "[upstream, downstream] flank around the exon in bp."),
    ("target_mode", "exclude",
     "'exclude' to skip the target exon, 'include' to retain it."),
]


def missing_opinionated_fields(raw: dict) -> list[dict]:
    """Return entries for opinionated fields absent from `raw` (the loaded JSON).

    Each entry is `{"name", "default", "description"}`. Use this before
    running `predict_aso_efficacy` to prompt the user to confirm defaults
    rather than silently applying them.
    """
    return [
        {"name": name, "default": default, "description": desc}
        for name, default, desc in _OPINIONATED_FIELDS
        if name not in raw
    ]


@dataclass
class OligoConfig:
    gene_symbol: str
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

    exon_intervals: Optional[tuple[int, int]] = None
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

    Input paths (`fasta_path`, `experimental_data`) are resolved relative
    to the config file's directory so input bundles shipped next to the
    JSON stay portable. Output paths (`results_dir`, `data_dir`) are
    resolved relative to the current working directory, so running the
    tool from the project root writes results under `./results/`
    regardless of where the config lives.
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

    def _resolve_cwd(p: str, default: str) -> Path:
        pp = Path(p or default).expanduser()
        return pp.resolve() if pp.is_absolute() else (Path.cwd() / pp).resolve()

    if raw.get("target_mode", "exclude") not in ("exclude", "include"):
        raise ValueError(
            f"target_mode must be 'exclude' or 'include', got {raw.get('target_mode')!r}"
        )

    raw_exon = raw.get("exon_intervals")
    exon_intervals: Optional[tuple[int, int]] = (
        (int(raw_exon[0]), int(raw_exon[1])) if raw_exon else None
    )
    flank = raw.get("flank", [200, 200])

    cfg = OligoConfig(
        gene_symbol=raw["gene_symbol"],
        exon_intervals=exon_intervals,
        strand=raw.get("strand", "+"),
        assembly=raw.get("assembly", "hg38"),
        fasta_path=_resolve(raw.get("fasta_path")),
        gtf_url=raw.get("gtf_url") or _DEFAULT_GTF_URL,
        results_dir=_resolve_cwd(raw.get("results_dir"), "results"),
        data_dir=_resolve_cwd(raw.get("data_dir"), "data"),
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
