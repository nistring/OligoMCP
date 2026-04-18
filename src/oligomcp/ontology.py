"""Build + load the AlphaGenome ontology-term snapshot.

Queries `DnaClient.output_metadata(Organism.HOMO_SAPIENS)`, iterates over
every output type that ships biosample-level ontology metadata, and
emits a deduplicated table where each row is a unique `ontology_curie`
annotated with:

  - biosample_name / _type / _life_stage / gtex_tissue — human labels.
  - available_outputs — comma-joined list of AlphaGenome output types
    that have at least one track for this CURIE.
  - track_filters — pipe-joined list of assay-descriptor substrings
    (e.g. "polyA plus RNA-seq", "total RNA-seq", "ATAC-seq",
    "Histone ChIP-seq H3K27ac") that actually appear in track names
    for this CURIE. These are the values to pass as the `track_filter`
    field in OligoConfig / `predict_aso_efficacy_inline` — mismatched
    filters silently drop every AlphaGenome track for the CURIE.
  - total_tracks — sum of track rows across output types.

Exposes:

  - `build_ontology_table(client)` — pure builder, returns a DataFrame.
  - `save_ontology_snapshot(data_dir)` — writes TSV + meta file.
  - `load_ontology_snapshot()` — reads the committed TSV into a list
    of dicts for the MCP `search_ontology_terms` tool.

Used by:
  - `oligomcp fetch-ontology-terms` (CLI subcommand).
  - `mcp_server.search_ontology_terms` (for substring search).
"""
from __future__ import annotations

import datetime as _dt
import re
from pathlib import Path
from typing import Optional

import pandas as pd

# Output attributes on alphagenome.data.output_metadata.OutputMetadata that
# ship biosample-level ontology metadata.
_OUTPUT_ATTRS = [
    "rna_seq",
    "cage",
    "procap",
    "atac",
    "dnase",
    "chip_histone",
    "chip_tf",
    "splice_sites",
    "splice_site_usage",
    "splice_junctions",
    "contact_maps",
]

# Track-name prefix strippers. Names look like:
#   "CL:0000100 total RNA-seq"
#   "usage_CL:0000047 polyA plus RNA-seq"      (splice_site_usage)
#   "CLO:0014078 ATAC-seq"
#   "CL:0000084 Histone ChIP-seq H3K27ac"
# After stripping the leading `(usage_|ref_|junc_)?<CURIE>\s+`, the tail
# is the "assay descriptor" that users pass as `track_filter`.
_CURIE_PREFIX_RE = re.compile(
    r"^(?:usage_|ref_|junc_|junction_)?"
    r"(?:CL|CLO|UBERON|EFO|NTR|GO)(?::|_)\d+"
    r"\s+"
)


def _assay_descriptor(track_name: str) -> str:
    """Return the 'track_filter' substring for an AlphaGenome track name."""
    stripped = _CURIE_PREFIX_RE.sub("", str(track_name), count=1).strip()
    return stripped or str(track_name).strip()


def _project_root() -> Path:
    # oligomcp/ontology.py → oligomcp/ → src/ → repo root
    return Path(__file__).resolve().parents[2]


def default_snapshot_dir() -> Path:
    """Where the committed TSV lives (repo-root `data/`)."""
    return _project_root() / "data"


def default_tsv_path() -> Path:
    return default_snapshot_dir() / "alphagenome_ontology_terms.tsv"


def default_meta_path() -> Path:
    return default_snapshot_dir() / "alphagenome_ontology_terms.meta.txt"


def build_ontology_table(client) -> pd.DataFrame:
    """Fetch ontology metadata from AlphaGenome; return one row per CURIE.

    Args:
        client: An `alphagenome.models.dna_client.DnaClient` instance.

    Returns:
        DataFrame with columns: ontology_curie, biosample_name,
        biosample_type, biosample_life_stage, gtex_tissue,
        available_outputs, track_filters, total_tracks.
    """
    from alphagenome.models.dna_client import Organism

    meta = client.output_metadata(Organism.HOMO_SAPIENS)

    frames: list[pd.DataFrame] = []
    for attr in _OUTPUT_ATTRS:
        df = getattr(meta, attr, None)
        if df is None or not hasattr(df, "columns"):
            continue
        if "ontology_curie" not in df.columns:
            continue
        keep = ["ontology_curie"]
        for optional in (
            "name",
            "biosample_name",
            "biosample_type",
            "biosample_life_stage",
            "gtex_tissue",
        ):
            if optional in df.columns:
                keep.append(optional)
        piece = df[keep].dropna(subset=["ontology_curie"]).copy()
        piece["__output_type"] = attr
        piece["__track_filter"] = (
            piece["name"].map(_assay_descriptor) if "name" in piece.columns else ""
        )
        frames.append(piece)

    if not frames:
        raise RuntimeError(
            "No ontology metadata found — is the AlphaGenome API key valid?"
        )

    cat = pd.concat(frames, ignore_index=True)

    def _merge_filters(series: pd.Series) -> str:
        uniq = sorted({f for f in series.dropna().astype(str) if f})
        return "|".join(uniq)

    grouped = (
        cat.groupby("ontology_curie", as_index=False)
        .agg(
            biosample_name=("biosample_name", "first")
            if "biosample_name" in cat.columns else ("ontology_curie", "first"),
            biosample_type=("biosample_type", "first")
            if "biosample_type" in cat.columns else ("ontology_curie", "first"),
            biosample_life_stage=("biosample_life_stage", "first")
            if "biosample_life_stage" in cat.columns else ("ontology_curie", "first"),
            gtex_tissue=("gtex_tissue", "first")
            if "gtex_tissue" in cat.columns else ("ontology_curie", "first"),
            available_outputs=("__output_type", lambda s: ",".join(sorted(set(s)))),
            track_filters=("__track_filter", _merge_filters),
            total_tracks=("ontology_curie", "size"),
        )
        .sort_values("ontology_curie")
        .reset_index(drop=True)
    )
    return grouped


def save_ontology_snapshot(
    data_dir: Optional[Path] = None,
    *,
    client=None,
) -> tuple[Path, Path]:
    """Regenerate the snapshot TSV + meta file and return their paths.

    Writes as the last step so an interrupted run never leaves partial
    output in place.
    """
    from alphagenome.models.dna_client import create

    from .resources import require_alphagenome_api_key

    data_dir = Path(data_dir) if data_dir else default_snapshot_dir()
    data_dir.mkdir(parents=True, exist_ok=True)

    if client is None:
        client = create(require_alphagenome_api_key())

    df = build_ontology_table(client)

    tsv_path = data_dir / "alphagenome_ontology_terms.tsv"
    meta_path = data_dir / "alphagenome_ontology_terms.meta.txt"
    # Write to a sibling `.partial` first so an interruption never clobbers
    # an existing good snapshot. Atomic rename on success.
    tmp_tsv = tsv_path.with_suffix(tsv_path.suffix + ".partial")
    df.to_csv(tmp_tsv, sep="\t", index=False)
    tmp_tsv.replace(tsv_path)

    try:
        import alphagenome
        ag_version = getattr(alphagenome, "__version__", "unknown")
    except Exception:
        ag_version = "unknown"

    outputs = sorted(
        {o for row in df["available_outputs"] for o in str(row).split(",") if o}
    )
    meta_lines = [
        f"Generated at:         {_dt.datetime.now(_dt.timezone.utc).isoformat()}",
        f"AlphaGenome version:  {ag_version}",
        f"Unique ontology CURIEs: {len(df)}",
        f"Outputs covered:      {outputs}",
        "",
        "Columns:",
        "  ontology_curie         — pass this in `ontology_terms` list.",
        "  biosample_*            — human-readable labels.",
        "  available_outputs      — AlphaGenome output types with tracks for this CURIE.",
        "  track_filters          — assay descriptor substrings to pass as `track_filter`;",
        "                            pipe-separated. If `polyA plus RNA-seq` isn't listed",
        "                            for a CURIE, the default inline track_filter drops every",
        "                            AlphaGenome track for that term — use `total RNA-seq`",
        "                            or leave track_filter empty.",
        "  total_tracks           — row count across all output types.",
        "",
        "Empty `ontology_terms` = average across every track AlphaGenome returns.",
    ]
    meta_path.write_text("\n".join(meta_lines) + "\n", encoding="utf-8")

    return tsv_path, meta_path


_SNAPSHOT_CACHE: Optional[list[dict]] = None


def load_ontology_snapshot(path: Optional[Path] = None) -> list[dict]:
    """Read the committed snapshot into a list of dicts (cached)."""
    global _SNAPSHOT_CACHE
    if _SNAPSHOT_CACHE is not None and path is None:
        return _SNAPSHOT_CACHE
    p = Path(path) if path else default_tsv_path()
    if not p.exists():
        result: list[dict] = []
    else:
        df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        result = df.to_dict(orient="records")
    if path is None:
        _SNAPSHOT_CACHE = result
    return result


def clear_snapshot_cache() -> None:
    """Drop the in-memory cache (so the next `load_` re-reads the TSV)."""
    global _SNAPSHOT_CACHE
    _SNAPSHOT_CACHE = None
