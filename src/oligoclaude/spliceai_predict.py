"""CPU-efficient SpliceAI prediction via OpenSpliceAI.

The SpliceAI efficacy score follows the same `diff_mean` formula as the
AlphaGenome SPLICE_SITE_USAGE track (see alphagenome_predict.diff_mean),
applied to per-position donor+acceptor probabilities produced by the
OSAI_MANE-10000nt 5-model ensemble.

Performance strategy:
- One-hot encode the reference window (SL + CL_max) exactly once.
- Keep a preallocated torch tensor and clone+slice-assign per ASO variant.
- Batch B=12 candidates per ensemble forward pass.
- `torch.no_grad()`, `.eval()`, `torch.set_num_threads()` tuning.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np

from .aso_enum import AsoCandidate
from .config import OligoConfig
from .sequence_utils import load_reference_sequence, one_hot_encode

SL = 5000
CL_MAX = 10000
INPUT_LEN = SL + CL_MAX


@dataclass
class SpliceAIContext:
    models: list[Any]
    base_input: Any
    win_start_genomic: int
    exon_start_in_sl: int
    exon_end_in_sl: int


def setup_spliceai(
    threads: Optional[int] = None,
) -> tuple[list[Any], Any]:
    """Load the OpenSpliceAI MANE-10000nt 5-model ensemble on CPU.

    Returns (models, device). The models are in eval mode with gradients
    disabled. Threads default to half of available CPU cores.
    """
    import torch

    torch.set_num_threads(threads or max(1, (os.cpu_count() or 2) // 2))
    torch.set_grad_enabled(False)
    device = torch.device("cpu")

    models: list[Any] = []
    last_err: Optional[Exception] = None

    try:
        from openspliceai.predict.utils import (
            load_pytorch_models,
            initialize_constants,
        )

        consts = initialize_constants(flanking_size=CL_MAX)
        models_obj, _ = load_pytorch_models(None, device, consts.get("SL", SL), CL_MAX)
        if isinstance(models_obj, (list, tuple)):
            models = list(models_obj)
        else:
            models = [models_obj]
    except Exception as e:
        last_err = e
        try:
            from openspliceai.predict import load_models  # type: ignore

            models = list(load_models(flanking_size=CL_MAX, device=device))
        except Exception as e2:
            raise RuntimeError(
                "Could not load OpenSpliceAI models. Install with "
                "`pip install openspliceai` and ensure MANE-10000nt weights "
                f"are downloaded. Last errors: {last_err!r}, {e2!r}"
            ) from e2

    for m in models:
        m.eval()
        m.to(device)

    return models, device


def _build_base_tensor(
    fasta_path: Path, chrom: str, exon_start: int, exon_end: int
) -> tuple[Any, int, int, int]:
    """Center an (SL + CL_MAX) window on the exon midpoint and one-hot encode.

    Returns:
        (tensor, win_start_genomic, exon_start_in_sl, exon_end_in_sl)

        tensor: torch.Tensor shape (1, 4, INPUT_LEN) float32
        win_start_genomic: genomic position of input[0]
        exon_start_in_sl: exon start position within the SL output region (0..SL)
        exon_end_in_sl: exon end position within the SL output region
    """
    import torch

    exon_mid = (exon_start + exon_end) // 2
    sl_start = exon_mid - SL // 2
    sl_end = sl_start + SL
    win_start = sl_start - CL_MAX // 2
    win_end = win_start + INPUT_LEN

    seq = load_reference_sequence(fasta_path, chrom, win_start, win_end)
    if len(seq) != INPUT_LEN:
        raise RuntimeError(
            f"Requested {INPUT_LEN} bp from {chrom}:{win_start}-{win_end}, "
            f"got {len(seq)} bp. Exon too close to chromosome edge?"
        )

    arr = one_hot_encode(seq)
    tensor = torch.from_numpy(arr.T).unsqueeze(0).contiguous()

    exon_start_in_sl = max(0, exon_start - sl_start)
    exon_end_in_sl = min(SL, exon_end - sl_start)
    if exon_start_in_sl >= exon_end_in_sl:
        raise RuntimeError(
            f"Target exon {exon_start}-{exon_end} does not fall within the "
            f"SL window {sl_start}-{sl_end}. Increase SL or check coordinates."
        )

    return tensor, win_start, exon_start_in_sl, exon_end_in_sl


def _forward_ensemble(models: list[Any], batch):
    """Run a batch through all ensemble models and average the outputs.

    Input:  (B, 4, INPUT_LEN)
    Output: (B, 3, SL) — [none/background, acceptor, donor] probabilities
    """
    import torch

    outs = []
    for m in models:
        o = m(batch)
        if isinstance(o, (tuple, list)):
            o = o[0]
        outs.append(o)
    stacked = torch.stack(outs, dim=0)
    return stacked.mean(dim=0)


def score_asos_spliceai(
    cfg: OligoConfig,
    candidates: list[AsoCandidate],
    chrom: str,
    variant_interval_start_genomic: int,
    models: Optional[list[Any]] = None,
) -> np.ndarray:
    """Score each ASO candidate by diff_mean on donor+acceptor probabilities.

    The score uses the same `diff_mean` formula as AlphaGenome:
        alt[exon].mean() / ref[SL].mean() * alt[SL].mean() - ref[exon].mean()
    where `ref`/`alt` are (SL,) per-position donor+acceptor probability sums
    for the reference and variant respectively.

    `candidates[i].position` is the 0-based offset within the scan region
    (variant_interval), so absolute genomic position =
        variant_interval_start_genomic + candidates[i].position.
    """
    import torch

    if models is None:
        models, _ = setup_spliceai(cfg.spliceai_threads)

    exon_start = int(cfg.exon_intervals[0])
    exon_end = int(cfg.exon_intervals[1])
    base_tensor, win_start, exon_a, exon_b = _build_base_tensor(
        cfg.fasta_path, chrom, exon_start, exon_end
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

    for batch_start in range(0, len(candidates), batch_size):
        batch_cands = candidates[batch_start : batch_start + batch_size]
        B = len(batch_cands)
        batch_input = base_tensor.expand(B, -1, -1).clone()

        valid_mask = np.ones(B, dtype=bool)
        for i, cand in enumerate(batch_cands):
            cand_genomic_start = variant_interval_start_genomic + cand.position
            rel_start = cand_genomic_start - win_start
            rel_end = rel_start + cand.length
            if rel_start < 0 or rel_end > INPUT_LEN:
                valid_mask[i] = False
                continue
            batch_input[i, :, rel_start:rel_end] = 0.25

        alt_out = _forward_ensemble(models, batch_input)
        alt_signal = (alt_out[:, 1] + alt_out[:, 2]).cpu().numpy().astype(np.float32)

        alt_exon_mean = alt_signal[:, exon_a:exon_b].mean(axis=1)
        alt_sl_mean = alt_signal.mean(axis=1)

        batch_scores = (
            alt_exon_mean / (ref_sl_mean + 1e-12) * alt_sl_mean - ref_exon_mean
        )
        batch_scores[~valid_mask] = np.nan
        scores[batch_start : batch_start + B] = batch_scores.astype(np.float32)

        if (batch_start // batch_size) % 10 == 0 or batch_start + B >= len(candidates):
            done = min(batch_start + B, len(candidates))
            print(f"  SpliceAI progress: {done}/{len(candidates)}")

    return scores
