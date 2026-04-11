"""SpliceAI (OpenSpliceAI MANE-10000nt) weight cache.

OpenSpliceAI publishes the 5-model MANE-10000nt ensemble as standalone
PyTorch state-dict files on the Johns Hopkins CCB FTP mirror:
    ftp://ftp.ccb.jhu.edu/pub/data/OpenSpliceAI/OSAI-MANE/10000nt/

We download each .pt file to ~/.oligoclaude/spliceai/mane_10000nt/ so
Windows users (who cannot install the full openspliceai package because
of its pysam dependency) can still load the weights. If your network
blocks FTP, grab the files manually from the URL above and drop them in
the cache dir — `ensure_spliceai_weights()` skips already-present files.
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path
from typing import Optional

CACHE_DIR = Path.home() / ".oligoclaude" / "spliceai" / "mane_10000nt"

_MIRROR_BASES = [
    # JHU CCB FTP — the canonical publication location for the MANE weights.
    # Verified working 2025; the previously-tried HTTPS mirrors at
    # ccb.jhu.edu/openspliceai/data/... and ccb.jhu.edu/pub/data/... 404.
    "ftp://ftp.ccb.jhu.edu/pub/data/OpenSpliceAI/OSAI-MANE/10000nt",
]
_WEIGHT_FILES = [
    "model_10000nt_rs10.pt",
    "model_10000nt_rs11.pt",
    "model_10000nt_rs12.pt",
    "model_10000nt_rs13.pt",
    "model_10000nt_rs14.pt",
]


def default_weights_dir() -> Path:
    return CACHE_DIR


def _download_one(url: str, dst: Path, *, verbose: bool) -> bool:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "oligoclaude"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            total = resp.headers.get("Content-Length")
            total_int = int(total) if total and total.isdigit() else None
            tmp = dst.with_suffix(dst.suffix + ".partial")
            written = 0
            last_pct = -1
            with open(tmp, "wb") as f:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    f.write(chunk)
                    written += len(chunk)
                    if verbose and total_int:
                        pct = int(written * 100 / total_int)
                        if pct != last_pct:
                            sys.stderr.write(f"\r    {dst.name}: {pct}%")
                            sys.stderr.flush()
                            last_pct = pct
            if verbose and total_int:
                sys.stderr.write("\n")
            tmp.replace(dst)
        return True
    except Exception as e:
        if verbose:
            print(f"    failed {url}: {e}")
        if dst.with_suffix(dst.suffix + ".partial").exists():
            try:
                dst.with_suffix(dst.suffix + ".partial").unlink()
            except OSError:
                pass
        return False


def ensure_spliceai_weights(
    cache_dir: Optional[Path] = None, *, verbose: bool = True
) -> Path:
    """Download the 5-model MANE-10000nt ensemble if not already cached.

    Returns the directory containing the .pt files.
    """
    cache_dir = Path(cache_dir) if cache_dir else CACHE_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)

    missing = [f for f in _WEIGHT_FILES if not (cache_dir / f).exists()]
    if not missing:
        if verbose:
            print(f"Using cached SpliceAI weights: {cache_dir}")
        return cache_dir

    if verbose:
        print(f"Downloading {len(missing)} SpliceAI weight file(s) to {cache_dir}")
    for fname in missing:
        dst = cache_dir / fname
        downloaded = False
        for base in _MIRROR_BASES:
            url = f"{base}/{fname}"
            if verbose:
                print(f"  {fname} ← {url}")
            if _download_one(url, dst, verbose=verbose):
                downloaded = True
                break
        if not downloaded:
            raise RuntimeError(
                f"Could not download {fname}. Tried: "
                + ", ".join(f"{b}/{fname}" for b in _MIRROR_BASES)
                + "\nCheck your internet connection or download the file manually "
                f"and place it in {cache_dir}."
            )
    return cache_dir
