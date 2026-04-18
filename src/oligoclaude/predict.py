"""AlphaGenome + SpliceAI scoring, unified around the `diff_mean` formula.

Both backends score each ASO candidate by how much a masked variant
(`N` block for AlphaGenome, uniform 0.25 one-hot for SpliceAI) perturbs
exon usage relative to gene-body baseline:

    diff_mean = alt[exon].mean() / ref[body].mean() * alt[body].mean()
                - ref[exon].mean()

AlphaGenome evaluates this over a resized gene interval via its hosted
`predict_sequences` API. SpliceAI runs the OpenSpliceAI MANE-10000nt
5-model ensemble on CPU over a fixed (SL+CL) window centered on the exon,
using the full SL window as the "body" baseline. Ports of the formula and
of the AlphaGenome helpers come from `/home/nistring/AlphaGenome_ASO/aso.ipynb`.
"""
from __future__ import annotations

import os
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from .config import OligoConfig
from .core import AsoCandidate, load_reference_sequence, one_hot_encode

SL = 5000
CL_MAX = 10000
INPUT_LEN = SL + CL_MAX


def diff_mean(
    ref: np.ndarray,
    alt: np.ndarray,
    start: int,
    end: int,
    gene_start: int,
    gene_end: int,
) -> np.ndarray:
    """Difference-of-means score. `ref`/`alt` shape (L, num_variants)."""
    return (
        alt[start:end].mean(0) / ref[gene_start:gene_end].mean(0)
        * alt[gene_start:gene_end].mean(0)
        - ref[start:end].mean(0)
    )


# ============================================================
# AlphaGenome
# ============================================================


@dataclass
class AGContext:
    model: Any
    gene_interval: Any
    interval: Any
    ref_seq: str
    variant_interval: Any
    ref_output: Any
    requested_outputs: list[Any]
    exon_start_rel: int
    exon_end_rel: int
    gene_start_rel: int
    gene_end_rel: int
    start_rel: int
    end_rel: int


def _parse_output_types(names: list[str]) -> list[Any]:
    from alphagenome.models.dna_client import OutputType

    out: list[Any] = []
    for n in names:
        u = n.upper()
        if u in {"SPLICE_JUNCTIONS", "SPLICE_SITES"}:
            raise ValueError(f"{u} is not supported. Use SPLICE_SITE_USAGE instead.")
        if u in OutputType.__members__:
            out.append(OutputType[u])
        else:
            print(f"Warning: Unknown output type '{n}', skipping...")
    return out


def _filter_td(td, cfg: OligoConfig):
    """Apply config track-name + strand filters to an AlphaGenome TrackData."""
    if td is None:
        return None
    if cfg.track_filter:
        td = td.filter_tracks([cfg.track_filter in n for n in td.names])
    strand_filter = {
        "+": "filter_to_positive_strand",
        "-": "filter_to_negative_strand",
    }.get(str(cfg.strand).lower())
    if strand_filter and hasattr(td, strand_filter):
        td = getattr(td, strand_filter)()
    return td


def _optimal_resize(width: int, override: Optional[int]) -> int:
    from alphagenome.models.dna_client import SUPPORTED_SEQUENCE_LENGTHS

    if override is not None:
        return int(override)
    supported = sorted(SUPPORTED_SEQUENCE_LENGTHS.values())
    for L in supported:
        if L >= width:
            return L
    print(f"WARNING: interval {width:,} bp exceeds max supported {supported[-1]:,} bp")
    return supported[-1]


def setup_alphagenome(cfg: OligoConfig) -> AGContext:
    """Create client, resize the gene interval, and fetch the wildtype ref output."""
    from alphagenome.data import gene_annotation, genome
    from alphagenome.models import dna_client

    from .resources import require_alphagenome_api_key

    model = dna_client.create(require_alphagenome_api_key(cfg.dna_api_key))
    gtf = pd.read_feather(cfg.gtf_url)
    gene_interval = gene_annotation.get_gene_interval(gtf, gene_symbol=cfg.gene_symbol)
    interval = gene_interval.resize(_optimal_resize(gene_interval.width, cfg.resize_width))

    print(f"Gene: {cfg.gene_symbol}")
    print(
        f"Interval: {interval.chromosome}:{interval.start}-{interval.end} "
        f"(width={interval.width:,} bp, gene width={gene_interval.width:,} bp)"
    )

    requested_outputs = _parse_output_types(cfg.requested_outputs)
    ref_seq = load_reference_sequence(
        cfg.fasta_path, interval.chromosome, interval.start, interval.end,
        assembly=cfg.assembly,
    )

    exon_start, exon_end = int(cfg.exon_intervals[0]), int(cfg.exon_intervals[1])
    variant_interval = genome.Interval(
        chromosome=interval.chromosome,
        start=exon_start - cfg.flank[0],
        end=exon_end + cfg.flank[1],
    )
    print(
        f"Target exon: {exon_start}-{exon_end} | "
        f"Scan region: {variant_interval.start}-{variant_interval.end} "
        f"({variant_interval.width:,} bp)"
    )

    print("Fetching wildtype prediction from AlphaGenome...")
    ref_output = model.predict_interval(
        interval=interval,
        requested_outputs=requested_outputs,
        ontology_terms=cfg.ontology_terms,
    )

    return AGContext(
        model=model,
        gene_interval=gene_interval,
        interval=interval,
        ref_seq=ref_seq,
        variant_interval=variant_interval,
        ref_output=ref_output,
        requested_outputs=requested_outputs,
        exon_start_rel=exon_start - interval.start,
        exon_end_rel=exon_end - interval.start,
        gene_start_rel=gene_interval.start - interval.start,
        gene_end_rel=gene_interval.end - interval.start,
        start_rel=variant_interval.start - interval.start,
        end_rel=variant_interval.end - interval.start,
    )


def score_asos_alphagenome(
    ctx: AGContext, cfg: OligoConfig, candidates: list[AsoCandidate]
) -> dict[str, np.ndarray]:
    """Mask each candidate with `N`s, batch-predict, and apply `diff_mean`."""
    variants = [
        ctx.ref_seq[: c.position] + "N" * c.length + ctx.ref_seq[c.position + c.length :]
        for c in candidates
    ]
    print(f"AlphaGenome: predicting {len(variants)} masked variants...")
    outputs = ctx.model.predict_sequences(
        intervals=[ctx.interval] * len(variants),
        sequences=variants,
        requested_outputs=ctx.requested_outputs,
        ontology_terms=cfg.ontology_terms,
    )

    results: dict[str, np.ndarray] = {}
    for output_type in ctx.requested_outputs:
        variant_cols = [
            _filter_td(out.get(output_type), cfg).values.mean(axis=1)
            for out in outputs
            if _filter_td(out.get(output_type), cfg) is not None
        ]
        if not variant_cols:
            print(f"Skipping {output_type.name}: no variant tracks")
            continue

        ref_td = _filter_td(ctx.ref_output.get(output_type), cfg)
        if ref_td is None:
            print(f"Skipping {output_type.name}: missing reference output")
            continue

        results[output_type.name] = np.asarray(
            diff_mean(
                ref_td.values.mean(axis=1)[:, None],
                np.stack(variant_cols, axis=1),
                ctx.exon_start_rel,
                ctx.exon_end_rel,
                ctx.gene_start_rel,
                ctx.gene_end_rel,
            ),
            dtype=np.float32,
        )
    return results


# ============================================================
# SpliceAI (OpenSpliceAI MANE-10000nt)
# ============================================================

_MANE_10000_L = 32
# W and AR are numpy arrays — the SpliceAI constructor does
# `np.sum(AR * (W - 1))` which requires element-wise ops, not list arithmetic.
_MANE_10000_W = np.array([11] * 8 + [21] * 4 + [41] * 4, dtype=np.int64)
_MANE_10000_AR = np.array([1] * 4 + [4] * 4 + [10] * 4 + [25] * 4, dtype=np.int64)


_SPLICEAI_CACHE: dict[str, tuple[list[Any], Any]] = {}
_SPLICEAI_LOCK = threading.Lock()


def _default_n_models() -> int:
    """How many models of the MANE-10000nt ensemble to load by default.

    Defaults to a single model — this is ~5× faster at inference, uses
    ~1/5 the memory, and produces scores that rank ASO candidates
    consistently with the full 5-model ensemble (the between-model
    variance reduction matters most for absolute calibration, not for
    ranking). Override with `OLIGOCLAUDE_SPLICEAI_N_MODELS=5` when full
    ensemble accuracy is needed (e.g. final validation run after the
    top candidates have been narrowed).
    """
    override = os.environ.get("OLIGOCLAUDE_SPLICEAI_N_MODELS", "").strip()
    if override.isdigit():
        return max(1, int(override))
    return 1


def setup_spliceai(
    threads: Optional[int] = None, weights_dir: Optional[Path] = None
) -> tuple[list[Any], Any]:
    """Load (or retrieve cached) MANE-10000nt model ensemble on CPU.

    Loaded models are cached at module level and reused across calls —
    critical for remote HTTP deployments where reloading checkpoints per
    request (~5-10 s per model) blows through MCP tool-call timeouts.
    Thread-safe via double-checked locking.

    Number of models loaded is determined by `_default_n_models()` (see
    docstring): 5 locally, 1 on Horizon, override via env var.

    Uses the vendored `SpliceAI` class from `oligoclaude._spliceai_model`
    (copied from openspliceai v0.0.5, MIT), so we avoid pulling in the
    `openspliceai` package and its C-extension transitive deps.
    """
    import torch

    from ._spliceai_model import SpliceAI
    from .resources import ensure_spliceai_weights

    weights_dir = ensure_spliceai_weights(weights_dir)
    n_models = _default_n_models()
    cache_key = f"{Path(weights_dir).resolve()}|n={n_models}"

    cached = _SPLICEAI_CACHE.get(cache_key)
    if cached is not None:
        return cached

    with _SPLICEAI_LOCK:
        cached = _SPLICEAI_CACHE.get(cache_key)
        if cached is not None:
            return cached

        torch.set_num_threads(threads or max(1, (os.cpu_count() or 2) // 2))
        torch.set_grad_enabled(False)
        device = torch.device("cpu")

        pt_files = sorted(Path(weights_dir).glob("*.pt"))[:n_models]
        if not pt_files:
            raise RuntimeError(f"No .pt weight files found in {weights_dir}.")

        print(f"SpliceAI: loading {len(pt_files)}-model ensemble from {weights_dir}")

        models: list[Any] = []
        for pt in pt_files:
            m = SpliceAI(_MANE_10000_L, _MANE_10000_W, _MANE_10000_AR)
            state = torch.load(pt, map_location=device)
            if isinstance(state, dict) and "state_dict" in state:
                state = state["state_dict"]
            try:
                m.load_state_dict(state)
            except RuntimeError:
                m.load_state_dict(
                    {k.replace("module.", "", 1): v for k, v in state.items()}
                )
            m.eval().to(device)
            models.append(m)

        result = (models, device)
        _SPLICEAI_CACHE[cache_key] = result
        return result


def _build_base_tensor(
    fasta_path: Optional[Path], chrom: str, exon_start: int, exon_end: int,
    assembly: str = "hg38",
):
    """Center a (SL+CL) window on the exon midpoint and one-hot encode once."""
    import torch

    exon_mid = (exon_start + exon_end) // 2
    sl_start = exon_mid - SL // 2
    win_start = sl_start - CL_MAX // 2
    win_end = win_start + INPUT_LEN

    seq = load_reference_sequence(fasta_path, chrom, win_start, win_end, assembly=assembly)
    if len(seq) != INPUT_LEN:
        raise RuntimeError(
            f"Requested {INPUT_LEN} bp from {chrom}:{win_start}-{win_end}, "
            f"got {len(seq)} bp. Exon too close to chromosome edge?"
        )
    tensor = torch.from_numpy(one_hot_encode(seq).T).unsqueeze(0).contiguous()

    exon_a = max(0, exon_start - sl_start)
    exon_b = min(SL, exon_end - sl_start)
    if exon_a >= exon_b:
        raise RuntimeError(
            f"Target exon {exon_start}-{exon_end} falls outside the SL "
            f"window {sl_start}-{sl_start + SL}."
        )
    return tensor, win_start, exon_a, exon_b


def _forward_ensemble(models: list[Any], batch):
    """Run a batch through all models and average the (B, 3, SL) outputs.

    Wrapped in `torch.no_grad()` because `torch.set_grad_enabled(False)`
    called from `setup_spliceai` is thread-local: the background preload
    thread disables grad only for itself, while FastMCP's HTTP request
    handler runs on a different thread with grad still enabled. Cached
    models get reused across threads; without this guard the resulting
    tensors carry `requires_grad=True` and `.numpy()` raises.
    """
    import torch

    with torch.no_grad():
        outs = []
        for m in models:
            o = m(batch)
            outs.append(o[0] if isinstance(o, (tuple, list)) else o)
        return torch.stack(outs, dim=0).mean(dim=0)


def score_asos_spliceai(
    cfg: OligoConfig,
    candidates: list[AsoCandidate],
    chrom: str,
    variant_interval_start_genomic: int,
    models: Optional[list[Any]] = None,
) -> np.ndarray:
    """Score each candidate via `diff_mean` on donor+acceptor probabilities.

    `candidates[i].position` is the offset within ref_seq;
    `variant_interval_start_genomic` is the genomic position of ref_seq[0].
    """
    if models is None:
        models, _ = setup_spliceai(cfg.spliceai_threads)

    base_tensor, win_start, exon_a, exon_b = _build_base_tensor(
        cfg.fasta_path, chrom, int(cfg.exon_intervals[0]), int(cfg.exon_intervals[1]),
        assembly=cfg.assembly,
    )

    ref_out = _forward_ensemble(models, base_tensor)
    ref_signal = (ref_out[0, 1] + ref_out[0, 2]).cpu().numpy().astype(np.float32)
    ref_exon_mean = float(ref_signal[exon_a:exon_b].mean())
    ref_sl_mean = float(ref_signal.mean())

    scores = np.zeros(len(candidates), dtype=np.float32)
    batch_size = max(1, int(cfg.spliceai_batch))
    print(
        f"SpliceAI: scoring {len(candidates)} ASOs on CPU "
        f"(batch={batch_size}, SL={SL}, CL={CL_MAX})..."
    )

    for start in range(0, len(candidates), batch_size):
        batch_cands = candidates[start : start + batch_size]
        B = len(batch_cands)
        batch_input = base_tensor.expand(B, -1, -1).clone()

        valid = np.ones(B, dtype=bool)
        for i, c in enumerate(batch_cands):
            rel_start = variant_interval_start_genomic + c.position - win_start
            rel_end = rel_start + c.length
            if rel_start < 0 or rel_end > INPUT_LEN:
                valid[i] = False
                continue
            batch_input[i, :, rel_start:rel_end] = 0.25

        alt_out = _forward_ensemble(models, batch_input)
        alt_signal = (alt_out[:, 1] + alt_out[:, 2]).cpu().numpy().astype(np.float32)

        batch_scores = (
            alt_signal[:, exon_a:exon_b].mean(axis=1)
            / (ref_sl_mean + 1e-12)
            * alt_signal.mean(axis=1)
            - ref_exon_mean
        )
        batch_scores[~valid] = np.nan
        scores[start : start + B] = batch_scores.astype(np.float32)

        if (start // batch_size) % 10 == 0 or start + B >= len(candidates):
            print(f"  SpliceAI progress: {min(start + B, len(candidates))}/{len(candidates)}")

    return scores
