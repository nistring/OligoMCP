"""Sequence utilities: reverse complement, FASTA I/O, one-hot encoding."""
from __future__ import annotations

from pathlib import Path

import numpy as np

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def load_reference_sequence(fasta_path: Path, chrom: str, start: int, end: int) -> str:
    """Load a genomic subsequence from an indexed FASTA file using pyfaidx.

    Coordinates are 0-based, end-exclusive (matching Python slicing and BED).
    """
    from pyfaidx import Fasta

    fa = Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True)
    if chrom not in fa:
        raise KeyError(
            f"Chromosome {chrom!r} not found in {fasta_path}. "
            f"Available: {list(fa.keys())[:5]}..."
        )
    return str(fa[chrom][start:end])


_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def one_hot_encode(seq: str) -> np.ndarray:
    """One-hot encode a DNA sequence.

    Returns an (L, 4) float32 array. A/C/G/T map to standard one-hot vectors;
    N and any unknown base maps to [0.25, 0.25, 0.25, 0.25] (uniform).
    """
    L = len(seq)
    arr = np.zeros((L, 4), dtype=np.float32)
    for i, b in enumerate(seq.upper()):
        idx = _BASE_IDX.get(b)
        if idx is None:
            arr[i] = 0.25
        else:
            arr[i, idx] = 1.0
    return arr
