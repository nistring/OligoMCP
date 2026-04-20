"""AlphaGenome + SpliceAI scoring, unified around the `diff_mean` formula.

Both backends score each ASO candidate by how much a masked variant
(`N` block for AlphaGenome, uniform 0.25 one-hot for SpliceAI) perturbs
exon usage relative to a gene-body baseline:

    diff_mean = alt[target].mean() / ref[body].mean() * alt[body].mean()
                - ref[target].mean()

where `target` means different things for different output tracks:

  - AlphaGenome RNA_SEQ: the whole target exon interval.
  - AlphaGenome SPLICE_SITE_USAGE: same interval, but the track is
    already peaked at the donor/acceptor so the mean captures mostly
    junction signal.
  - SpliceAI: only the two splice junctions (acceptor at exon_start,
    donor at exon_end-1). SpliceAI's softmax output is ~0 everywhere
    except at canonical splice sites, so averaging over the full exon
    would dilute the real signal with ~200+ near-zero bases.

AlphaGenome evaluates this over a resized gene interval via its hosted
`predict_sequences` API. SpliceAI runs the OpenSpliceAI MANE-10000nt
model (a single model by default; set `OLIGOMCP_SPLICEAI_N_MODELS=5` to
use the full 5-seed ensemble) on CPU over a fixed (SL+CL) window
centered on the exon, using the full SL window as the "body" baseline.
Ports of the formula and of the AlphaGenome helpers come from
`/home/nistring/AlphaGenome_ASO/aso.ipynb`.
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


_FRAC_EPS = 1e-6     # stabilizer for (raw / |ref[target]|) when baseline is ~0
_BODY_EPS = 1e-12    # stabilizer for alt[body] / ref[body]


def diff_mean_frac(
    ref: np.ndarray,
    alt: np.ndarray,
    target_sel,
    body_sel,
) -> np.ndarray:
    """Body-corrected fractional change of `alt` vs `ref` at `target`.

        α    = alt[body].mean / (ref[body].mean + BODY_EPS)
        raw  = α · alt[target].mean − ref[target].mean
        frac = raw / (|ref[target].mean| + FRAC_EPS)

    `α` is a global amplitude correction: meaningful for coverage-like
    tracks whose body mean shifts under masking, ≈ 1 for sparse tracks
    where body is mostly zero.

    `target_sel` / `body_sel` may be a `slice` or an int index array.
    `ref` shape (L,) or (L, 1); `alt` shape (L,) or (L, num_variants).
    Returned array shape matches the 2nd axis of `alt` (scalar or
    num_variants,).

    Returns the dimensionless fractional-of-baseline score. Negative =
    ASO reduces the reference signal (skip-promoting); positive =
    amplifies it (inclusion-enhancing).
    """
    ref_target = ref[target_sel].mean(0)
    alt_target = alt[target_sel].mean(0)
    ref_body = ref[body_sel].mean(0)
    alt_body = alt[body_sel].mean(0)
    alpha = alt_body / (ref_body + _BODY_EPS)
    raw = alpha * alt_target - ref_target
    return raw / (np.abs(ref_target) + _FRAC_EPS)


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
    """Apply config track-name + strand filters to an AlphaGenome TrackData.

    For a + strand gene we keep tracks that are NOT on the positive strand
    (i.e. `-` strand + unstranded), and vice versa. The ASO's antisense is
    opposite to the gene's mRNA sense, so the tracks that actually report
    on the ASO's target are on the opposite strand (plus any unstranded
    library, which pools both).
    """
    if td is None:
        return None
    if cfg.track_filter:
        td = td.filter_tracks([cfg.track_filter in n for n in td.names])
    strand_filter = {
        "+": "filter_to_nonpositive_strand",
        "-": "filter_to_nonnegative_strand",
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
    max_workers = max(1, int(getattr(cfg, "alphagenome_workers", 16) or 16))
    print(
        f"AlphaGenome: predicting {len(variants)} masked variants "
        f"(max_workers={max_workers})..."
    )
    outputs = ctx.model.predict_sequences(
        intervals=[ctx.interval] * len(variants),
        sequences=variants,
        requested_outputs=ctx.requested_outputs,
        ontology_terms=cfg.ontology_terms,
        max_workers=max_workers,
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

        ref_arr = ref_td.values.mean(axis=1)[:, None]
        alt_arr = np.stack(variant_cols, axis=1)

        # SPLICE_SITE_USAGE is peaked at the donor/acceptor — averaging
        # over the full exon dilutes the real signal with ~200 near-zero
        # bases, same reason SpliceAI now scores only at junctions.
        if output_type.name == "SPLICE_SITE_USAGE":
            target_sel = np.array(
                [ctx.exon_start_rel, ctx.exon_end_rel - 1], dtype=int
            )
        else:
            target_sel = slice(ctx.exon_start_rel, ctx.exon_end_rel)
        body_sel = slice(ctx.gene_start_rel, ctx.gene_end_rel)
        frac = diff_mean_frac(ref_arr, alt_arr, target_sel, body_sel)
        results[output_type.name] = np.asarray(frac, dtype=np.float32)
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
    ranking). Override with `OLIGOMCP_SPLICEAI_N_MODELS=5` when full
    ensemble accuracy is needed (e.g. final validation run after the
    top candidates have been narrowed).
    """
    override = os.environ.get("OLIGOMCP_SPLICEAI_N_MODELS", "").strip()
    if override.isdigit():
        return max(1, int(override))
    return 1


def _resolve_spliceai_device() -> Any:
    """Pick the torch device for SpliceAI.

    Default: CUDA if `torch.cuda.is_available()`, otherwise CPU.
    Override via `OLIGOMCP_SPLICEAI_DEVICE=cpu|cuda` — useful when the
    installed CUDA build mismatches the host driver (torch returns
    `is_available()=False` but we want an explicit error) or when
    forcing CPU for deterministic benchmarks.
    """
    import torch

    override = os.environ.get("OLIGOMCP_SPLICEAI_DEVICE", "").strip().lower()
    if override == "cpu":
        return torch.device("cpu")
    if override == "cuda":
        if not torch.cuda.is_available():
            print(
                "WARNING: OLIGOMCP_SPLICEAI_DEVICE=cuda but CUDA is not available; "
                "falling back to CPU."
            )
            return torch.device("cpu")
        return torch.device("cuda")
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def setup_spliceai(
    threads: Optional[int] = None, weights_dir: Optional[Path] = None
) -> tuple[list[Any], Any]:
    """Load (or retrieve cached) MANE-10000nt SpliceAI models.

    Auto-selects CUDA when available (override via
    `OLIGOMCP_SPLICEAI_DEVICE=cpu`). Loaded models are cached at module
    level and reused across every subsequent call, saving ~5-10 s per
    model on each request that would otherwise re-read the checkpoints.
    Thread-safe via double-checked locking (the MCP server may dispatch
    tool calls from multiple threads).

    Number of models loaded is determined by `_default_n_models()` —
    one model by default, set `OLIGOMCP_SPLICEAI_N_MODELS=5` for the
    full ensemble.

    Uses the vendored `SpliceAI` class from `oligomcp._spliceai_model`
    (copied from openspliceai v0.0.5, MIT), so we avoid pulling in the
    `openspliceai` package and its C-extension transitive deps.
    """
    import torch

    from ._spliceai_model import SpliceAI
    from .resources import ensure_spliceai_weights

    weights_dir = ensure_spliceai_weights(weights_dir)
    n_models = _default_n_models()
    device = _resolve_spliceai_device()
    # Device is part of the cache key so a CPU-preloaded model is not
    # accidentally reused when a later request forces CUDA (or vice versa).
    cache_key = f"{Path(weights_dir).resolve()}|n={n_models}|dev={device}"

    cached = _SPLICEAI_CACHE.get(cache_key)
    if cached is not None:
        return cached

    with _SPLICEAI_LOCK:
        cached = _SPLICEAI_CACHE.get(cache_key)
        if cached is not None:
            return cached

        if device.type == "cpu":
            torch.set_num_threads(threads or max(1, (os.cpu_count() or 2) // 2))
        torch.set_grad_enabled(False)

        pt_files = sorted(Path(weights_dir).glob("*.pt"))[:n_models]
        if not pt_files:
            raise RuntimeError(f"No .pt weight files found in {weights_dir}.")

        print(
            f"SpliceAI: loading {len(pt_files)}-model ensemble on {device} "
            f"from {weights_dir}"
        )

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
    applied_variants: Optional[list] = None,
):
    """Center a (SL+CL) window on the exon midpoint and one-hot encode once.

    When ``applied_variants`` is given, variants whose genomic positions
    fall inside the 15 kb window are edited into the loaded WT sequence
    before one-hot encoding. Indels change the patient window length;
    we pad with downstream reference or trim the right edge to keep
    exactly INPUT_LEN bp so the rest of the SpliceAI machinery still
    sees a 1×4×15000 tensor. ``exon_a`` / ``exon_b`` are returned in
    patient-window offsets so the scoring indexes the correct
    (post-edit) junction bases.
    """
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

    exon_a = max(0, exon_start - sl_start)
    exon_b = min(SL, exon_end - sl_start)
    local_coord_map = None

    if applied_variants:
        from .variants import apply_variants_to_ref, pad_or_trim_to_length

        in_window = [
            v for v in applied_variants
            if v.chrom == chrom
            and win_start <= (v.position - 1) < win_end
            and (v.position - 1 + len(v.ref)) <= win_end
        ]
        if in_window:
            seq, local_coord_map = apply_variants_to_ref(
                seq, in_window, anchor_genomic=win_start, chrom=chrom,
            )
            seq = pad_or_trim_to_length(
                seq, target=INPUT_LEN,
                fetcher=lambda c, a, b: load_reference_sequence(
                    fasta_path, c, a, b, assembly=assembly,
                ),
                chrom=chrom,
                anchor_genomic=win_start,
                original_length=INPUT_LEN,
            )
            # Remap exon boundaries into patient-window coords. If either
            # junction was deleted by a variant, scoring becomes
            # meaningless — fall through to the WT exon_a/exon_b and let
            # downstream code surface NaN scores rather than crash.
            mapped_a = local_coord_map.ref_to_patient(win_start + exon_a)
            mapped_b = local_coord_map.ref_to_patient(win_start + exon_b)
            if mapped_a is not None and mapped_b is not None:
                exon_a = max(0, min(SL, mapped_a))
                exon_b = max(0, min(SL, mapped_b))

    tensor = torch.from_numpy(one_hot_encode(seq).T).unsqueeze(0).contiguous()

    if exon_a >= exon_b:
        raise RuntimeError(
            f"Target exon {exon_start}-{exon_end} falls outside the SL "
            f"window {sl_start}-{sl_start + SL}."
        )
    return tensor, win_start, exon_a, exon_b, local_coord_map


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


def _auto_spliceai_batch(device: Any) -> int:
    """Pick a default batch size when the user leaves `spliceai_batch` at 0.

    Larger batches amortize fixed per-call overhead (Python, CUDA launch,
    BLAS setup), so the sweet spot depends on hardware: GPUs can carry
    hundreds of 15-kbase windows per forward pass, CPUs are bottlenecked
    by BLAS and memory bandwidth long before that.
    """
    return 256 if getattr(device, "type", "cpu") == "cuda" else 64


def score_asos_spliceai(
    cfg: OligoConfig,
    candidates: list[AsoCandidate],
    chrom: str,
    variant_interval_start_genomic: int,
    models: Optional[list[Any]] = None,
    device: Any = None,
    applied_variants: Optional[list] = None,
    coord_map: Optional[object] = None,
) -> np.ndarray:
    """Score each candidate via `diff_mean_frac` on donor+acceptor probabilities.

    `candidates[i].position` is the offset within ref_seq;
    `variant_interval_start_genomic` is the genomic position of ref_seq[0].

    When ``applied_variants`` is provided, the 15 kb SpliceAI window is
    built from the patient baseline (WT window + variants edited in);
    ``coord_map`` is used to map each candidate's ASO target position
    from patient_seq coords (relative to the AlphaGenome interval) into
    SpliceAI-window coords, so the ASO mask lands on the correct patient
    bases.

    Returns a dimensionless per-candidate array. Negative = ASO reduces
    junction probability (skip-promoting).
    """
    if models is None:
        models, device = setup_spliceai(cfg.spliceai_threads)
    if device is None:
        # The loaded model carries the authoritative device; infer from it
        # so stale cached callers still work after the GPU path landed.
        device = next(models[0].parameters()).device

    base_tensor, win_start, exon_a, exon_b, local_coord_map = _build_base_tensor(
        cfg.fasta_path, chrom, int(cfg.exon_intervals[0]), int(cfg.exon_intervals[1]),
        assembly=cfg.assembly,
        applied_variants=applied_variants,
    )
    base_tensor = base_tensor.to(device)

    ref_out = _forward_ensemble(models, base_tensor)
    # SpliceAI softmax output: channel 1 = acceptor, channel 2 = donor.
    ref_signal = (ref_out[0, 1] + ref_out[0, 2]).cpu().numpy().astype(np.float32)
    # Score only the two junctions of the target exon: acceptor at
    # exon_start (first exonic base), donor at exon_end-1 (last exonic
    # base). Averaging over the 200+ bp exon would dilute the peaked
    # softmax signal with hundreds of ~0 bases.
    junction_idx = np.array([exon_a, exon_b - 1], dtype=int)
    # Shape (L, 1) so `diff_mean_frac` broadcasts cleanly with alt (L, B).
    ref_arr = ref_signal[:, None]

    scores = np.zeros(len(candidates), dtype=np.float32)
    requested_batch = int(cfg.spliceai_batch)
    batch_size = requested_batch if requested_batch > 0 else _auto_spliceai_batch(device)
    print(
        f"SpliceAI: scoring {len(candidates)} ASOs on {device} "
        f"(batch={batch_size}, SL={SL}, CL={CL_MAX})..."
    )

    for start in range(0, len(candidates), batch_size):
        batch_cands = candidates[start : start + batch_size]
        B = len(batch_cands)
        batch_input = base_tensor.expand(B, -1, -1).clone()

        valid = np.ones(B, dtype=bool)
        for i, c in enumerate(batch_cands):
            # c.position is a 0-based offset within ref_seq. For the
            # default (no variants) path that equals a reference offset;
            # for the patient-baseline path it's a patient offset inside
            # the AG interval, which we first translate to genomic via
            # the AG coord_map, then into SpliceAI-window coords via
            # local_coord_map (if variants reached the window).
            if coord_map is not None:
                gpos = coord_map.patient_to_ref(c.position)
                if gpos is None:
                    valid[i] = False
                    continue
            else:
                gpos = variant_interval_start_genomic + c.position

            if local_coord_map is not None:
                rel_start = local_coord_map.ref_to_patient(gpos)
                if rel_start is None:
                    valid[i] = False
                    continue
            else:
                rel_start = gpos - win_start
            rel_end = rel_start + c.length
            if rel_start < 0 or rel_end > INPUT_LEN:
                valid[i] = False
                continue
            batch_input[i, :, rel_start:rel_end] = 0.25

        alt_out = _forward_ensemble(models, batch_input)
        # alt_arr shape (L, B) so the same diff_mean_frac broadcasting pattern
        # used for AlphaGenome applies.
        alt_arr = (alt_out[:, 1] + alt_out[:, 2]).cpu().numpy().astype(np.float32).T
        frac = np.asarray(
            diff_mean_frac(ref_arr, alt_arr, junction_idx, slice(0, SL)),
            dtype=np.float32,
        )
        frac[~valid] = np.nan
        scores[start : start + B] = frac

        if (start // batch_size) % 10 == 0 or start + B >= len(candidates):
            print(f"  SpliceAI progress: {min(start + B, len(candidates))}/{len(candidates)}")

    return scores
