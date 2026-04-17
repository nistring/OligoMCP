"""External resources: API credentials, GRCh38 FASTA, SpliceAI weights.

All three are cached under ~/.oligoclaude/. The credentials file is written
mode 0600; the genome and weight caches are plain files. A single
`_download_with_progress` helper handles streaming + progress for both
downloads (with optional gzip decode for the genome).
"""
from __future__ import annotations

import gzip
import json
import os
import stat
import sys
import urllib.parse
import urllib.request
import warnings
from pathlib import Path
from typing import Optional

BASE_DIR = Path.home() / ".oligoclaude"
CRED_PATH = BASE_DIR / "credentials.json"
GENOME_DIR = BASE_DIR / "genomes"
HG38_FILENAME = "GRCh38.primary_assembly.genome.fa"
HG38_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_46/GRCh38.primary_assembly.genome.fa.gz"
)
SPLICEAI_DIR = BASE_DIR / "spliceai" / "mane_10000nt"
_SPLICEAI_MIRRORS = [
    "ftp://ftp.ccb.jhu.edu/pub/data/OpenSpliceAI/OSAI-MANE/10000nt",
]
_SPLICEAI_FILES = [f"model_10000nt_rs{i}.pt" for i in range(10, 15)]

ENV_VAR = "ALPHAGENOME_API_KEY"
_PLACEHOLDERS = {"", "REPLACE_WITH_YOUR_KEY", "YOUR_ALPHAGENOME_API_KEY"}


# ---------- credentials ----------

def _read_creds() -> dict:
    if not CRED_PATH.exists():
        return {}
    try:
        return json.loads(CRED_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def get_alphagenome_api_key(cfg_value: Optional[str] = None) -> Optional[str]:
    """Resolve the AlphaGenome key from env, file, or legacy config field."""
    env_key = os.environ.get(ENV_VAR, "").strip()
    if env_key and env_key not in _PLACEHOLDERS:
        return env_key

    file_key = _read_creds().get("alphagenome_api_key", "")
    if isinstance(file_key, str) and file_key.strip() and file_key.strip() not in _PLACEHOLDERS:
        return file_key.strip()

    if cfg_value and cfg_value.strip() and cfg_value.strip() not in _PLACEHOLDERS:
        warnings.warn(
            "dna_api_key is embedded in the config file. Move it to the "
            f"{ENV_VAR} environment variable or run `oligoclaude set-api-key` "
            "(saves to ~/.oligoclaude/credentials.json, mode 0600).",
            DeprecationWarning,
            stacklevel=2,
        )
        return cfg_value.strip()

    return None


def require_alphagenome_api_key(cfg_value: Optional[str] = None) -> str:
    key = get_alphagenome_api_key(cfg_value)
    if not key:
        raise RuntimeError(
            "No AlphaGenome API key found. Set one of:\n"
            f"  - {ENV_VAR} environment variable\n"
            f"  - Run: oligoclaude set-api-key <KEY>\n"
            "Or pass --skip-alphagenome to run without it."
        )
    return key


def save_alphagenome_api_key(key: str) -> Path:
    key = key.strip()
    if not key or key in _PLACEHOLDERS:
        raise ValueError("Refusing to save an empty or placeholder key.")
    BASE_DIR.mkdir(parents=True, exist_ok=True)
    data = _read_creds()
    data["alphagenome_api_key"] = key
    CRED_PATH.write_text(json.dumps(data, indent=2), encoding="utf-8")
    try:
        os.chmod(CRED_PATH, stat.S_IRUSR | stat.S_IWUSR)
    except OSError:
        pass
    return CRED_PATH


def clear_alphagenome_api_key() -> bool:
    data = _read_creds()
    if "alphagenome_api_key" not in data:
        return False
    del data["alphagenome_api_key"]
    if data:
        CRED_PATH.write_text(json.dumps(data, indent=2), encoding="utf-8")
    else:
        try:
            CRED_PATH.unlink()
        except FileNotFoundError:
            pass
    return True


# ---------- shared download ----------

def _download_with_progress(
    url: str, dst: Path, *, verbose: bool = True, decode_gzip: bool = False
) -> None:
    """Stream `url` into `dst`, optionally gunzipping on the fly.

    Writes to `<dst>.partial` then atomically renames. Prints percentage
    progress to stderr when Content-Length is known.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_suffix(dst.suffix + ".partial")
    req = urllib.request.Request(url, headers={"User-Agent": "oligoclaude"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        total = resp.headers.get("Content-Length") if hasattr(resp, "headers") else None
        total_int = int(total) if total and str(total).isdigit() else None
        src = gzip.GzipFile(fileobj=resp) if decode_gzip else resp
        written = 0
        last_pct = -1
        label = dst.name
        with open(tmp, "wb") as out:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
                written += len(chunk)
                if verbose and total_int and not decode_gzip:
                    pct = int(written * 100 / total_int)
                    if pct != last_pct:
                        sys.stderr.write(f"\r  {label}: {pct}%")
                        sys.stderr.flush()
                        last_pct = pct
        if verbose and total_int and not decode_gzip:
            sys.stderr.write("\n")
    tmp.replace(dst)


# ---------- GRCh38 FASTA ----------

def default_hg38_path() -> Path:
    return GENOME_DIR / HG38_FILENAME


def ensure_hg38_fasta(
    cache_dir: Optional[Path] = None, *, verbose: bool = True
) -> Path:
    """Download and stream-decompress GRCh38 if not already cached."""
    cache_dir = Path(cache_dir) if cache_dir else GENOME_DIR
    dst = cache_dir / HG38_FILENAME
    if dst.exists() and dst.stat().st_size > 2 * 1024**3:
        if verbose:
            print(f"Using cached GRCh38 FASTA: {dst}")
        return dst
    if verbose:
        print(f"Downloading GRCh38 primary assembly FASTA from {HG38_URL}")
        print(f"  → {dst} (streaming gunzip, ~900 MB download, ~3 GB on disk)")
    _download_with_progress(HG38_URL, dst, verbose=verbose, decode_gzip=True)
    if verbose:
        print(f"  Decompressed: {dst.stat().st_size / 1024**3:.2f} GB")
    return dst


def resolve_fasta_path(
    cfg_value: Optional[Path], *, verbose: bool = True
) -> Optional[Path]:
    """Resolve a FASTA path, returning None when no local file is available.

    When None is returned, callers should use `fetch_sequence_ucsc()` to
    retrieve only the needed region on-the-fly instead of downloading
    the entire ~3 GB genome.
    """
    if cfg_value and Path(cfg_value).exists():
        return Path(cfg_value)
    cached = default_hg38_path()
    if cached.exists() and cached.stat().st_size > 2 * 1024**3:
        if verbose and cfg_value:
            print(f"fasta_path {cfg_value} not found — using cached {cached}")
        return cached
    if verbose:
        print(
            "No local FASTA — genomic sequences will be fetched on-the-fly "
            "from the UCSC API. Run `oligoclaude fetch-genome` to cache "
            "the full genome locally for offline use."
        )
    return None


# ---------- UCSC online sequence fetch ----------

def lookup_gene_chromosome(gene_symbol: str, assembly: str = "hg38") -> str:
    """Look up a gene's chromosome via the mygene.info REST API.

    Returns the UCSC-style chromosome string (e.g. `chr3`). Works without
    AlphaGenome or a local GTF, so it unblocks `--skip-alphagenome` runs
    driven by only a gene symbol.
    """
    pos_field = "genomic_pos_hg19" if assembly.lower() in {"hg19", "grch37"} else "genomic_pos"
    params = urllib.parse.urlencode({
        "q": f"symbol:{gene_symbol}",
        "species": "human",
        "size": 1,
    })
    hits_url = f"https://mygene.info/v3/query?{params}"
    req = urllib.request.Request(hits_url, headers={"User-Agent": "oligoclaude/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            hits_data = json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        raise RuntimeError(
            f"Failed to look up gene {gene_symbol!r} on mygene.info: {e}"
        ) from e

    hits = hits_data.get("hits") or []
    if not hits:
        raise RuntimeError(f"mygene.info returned no hits for gene symbol {gene_symbol!r}")
    gene_id = hits[0].get("_id")
    if not gene_id:
        raise RuntimeError(f"mygene.info hit missing _id for {gene_symbol!r}: {hits[0]}")

    detail_url = f"https://mygene.info/v3/gene/{gene_id}"
    req = urllib.request.Request(detail_url, headers={"User-Agent": "oligoclaude/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            detail = json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        raise RuntimeError(
            f"Failed to fetch gene detail for {gene_symbol!r} (id={gene_id}): {e}"
        ) from e

    pos = detail.get(pos_field)
    if isinstance(pos, list):
        pos = pos[0] if pos else None
    if not pos or "chr" not in pos:
        raise RuntimeError(
            f"mygene.info has no {pos_field} for {gene_symbol!r} (id={gene_id})"
        )
    chrom = str(pos["chr"])
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def fetch_sequence_ucsc(assembly: str, chrom: str, start: int, end: int) -> str:
    """Fetch a genomic region from the UCSC REST API (0-based, end-exclusive).

    Returns the uppercase DNA string. Much faster than downloading the
    full ~3 GB genome FASTA — only retrieves the exact region needed.
    """
    params = urllib.parse.urlencode({
        "genome": assembly,
        "chrom": chrom,
        "start": start,
        "end": end,
    })
    url = f"https://api.genome.ucsc.edu/getData/sequence?{params}"
    req = urllib.request.Request(url, headers={"User-Agent": "oligoclaude/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        raise RuntimeError(
            f"Failed to fetch sequence from UCSC API for "
            f"{chrom}:{start}-{end} ({assembly}): {e}\n"
            "If you are offline, run `oligoclaude fetch-genome` to download "
            "the full GRCh38 FASTA for local use."
        ) from e

    seq = data.get("dna", "")
    if not seq:
        raise RuntimeError(
            f"UCSC API returned no sequence for {chrom}:{start}-{end} ({assembly}). "
            f"Response: {data}"
        )
    return seq.upper()


# ---------- SpliceAI weights ----------

def default_spliceai_weights_dir() -> Path:
    return SPLICEAI_DIR


def ensure_spliceai_weights(
    cache_dir: Optional[Path] = None, *, verbose: bool = True
) -> Path:
    """Download the MANE-10000nt 5-model ensemble if not already cached."""
    cache_dir = Path(cache_dir) if cache_dir else SPLICEAI_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)
    missing = [f for f in _SPLICEAI_FILES if not (cache_dir / f).exists()]
    if not missing:
        if verbose:
            print(f"Using cached SpliceAI weights: {cache_dir}")
        return cache_dir

    if verbose:
        print(f"Downloading {len(missing)} SpliceAI weight file(s) to {cache_dir}")
    for fname in missing:
        dst = cache_dir / fname
        last_err: Optional[Exception] = None
        for base in _SPLICEAI_MIRRORS:
            url = f"{base}/{fname}"
            if verbose:
                print(f"  {fname} ← {url}")
            try:
                _download_with_progress(url, dst, verbose=verbose)
                last_err = None
                break
            except Exception as e:
                last_err = e
                if verbose:
                    print(f"    failed: {e}")
        if last_err is not None:
            raise RuntimeError(
                f"Could not download {fname} from any mirror. Last error: {last_err!r}. "
                f"Download manually from "
                f"ftp://ftp.ccb.jhu.edu/pub/data/OpenSpliceAI/OSAI-MANE/10000nt/ "
                f"and place in {cache_dir}."
            )
    return cache_dir
