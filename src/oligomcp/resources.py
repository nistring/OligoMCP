"""External resources: API credentials, GRCh38 FASTA, SpliceAI weights.

All three are cached under ~/.oligomcp/. The credentials file is written
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

BASE_DIR = Path.home() / ".oligomcp"
CRED_PATH = BASE_DIR / "credentials.json"
GENOME_DIR = BASE_DIR / "genomes"
VARIANT_CACHE_DIR = BASE_DIR / "variant_cache"
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
            f"{ENV_VAR} environment variable or run `oligomcp set-api-key` "
            "(saves to ~/.oligomcp/credentials.json, mode 0600).",
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
            f"  - Run: oligomcp set-api-key <KEY>\n"
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
    req = urllib.request.Request(url, headers={"User-Agent": "oligomcp"})
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
            "from the UCSC API. Run `oligomcp fetch-genome` to cache "
            "the full genome locally for offline use."
        )
    return None


# ---------- UCSC online sequence fetch ----------

def _mygene_detail(gene_symbol: str) -> dict:
    """Resolve a gene symbol to its full mygene.info record."""
    params = urllib.parse.urlencode({
        "q": f"symbol:{gene_symbol}",
        "species": "human",
        "size": 1,
    })
    hits_url = f"https://mygene.info/v3/query?{params}"
    req = urllib.request.Request(hits_url, headers={"User-Agent": "oligomcp/1.0"})
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
    req = urllib.request.Request(detail_url, headers={"User-Agent": "oligomcp/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        raise RuntimeError(
            f"Failed to fetch gene detail for {gene_symbol!r} (id={gene_id}): {e}"
        ) from e


def lookup_gene_info(gene_symbol: str, assembly: str = "hg38") -> dict:
    """Resolve a gene symbol to chromosome, strand, and transcript exons.

    Returns a dict with keys: `chrom` (UCSC style, e.g. `chr3`),
    `strand` (`+` or `-`), `gene_start`, `gene_end`, `transcripts` (list of
    dicts with `transcript`, `cdsstart`, `cdsend`, `exons: [[start, end], ...]`).

    Uses `exons` for hg38 and `exons_hg19` for hg19. Works without
    AlphaGenome, so it unblocks NL-driven runs that specify only a gene.
    """
    detail = _mygene_detail(gene_symbol)
    pos_field = "genomic_pos_hg19" if assembly.lower() in {"hg19", "grch37"} else "genomic_pos"
    exons_field = "exons_hg19" if assembly.lower() in {"hg19", "grch37"} else "exons"

    pos = detail.get(pos_field)
    if isinstance(pos, list):
        pos = pos[0] if pos else None
    if not pos or "chr" not in pos:
        raise RuntimeError(f"mygene.info has no {pos_field} for {gene_symbol!r}")

    chrom = str(pos["chr"])
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    strand = "+" if int(pos.get("strand", 1)) >= 0 else "-"

    raw_tx = detail.get(exons_field) or []
    transcripts: list[dict] = []
    for tx in raw_tx:
        exons = tx.get("position") or []
        if not exons:
            continue
        transcripts.append({
            "transcript": tx.get("transcript"),
            "cdsstart": int(tx["cdsstart"]) if tx.get("cdsstart") is not None else None,
            "cdsend": int(tx["cdsend"]) if tx.get("cdsend") is not None else None,
            "exons": [[int(a), int(b)] for a, b in exons],
        })

    return {
        "chrom": chrom,
        "strand": strand,
        "gene_start": int(pos.get("start", 0)),
        "gene_end": int(pos.get("end", 0)),
        "transcripts": transcripts,
    }


def lookup_gene_chromosome(gene_symbol: str, assembly: str = "hg38") -> str:
    """Look up just the chromosome for a gene symbol (convenience wrapper)."""
    return lookup_gene_info(gene_symbol, assembly)["chrom"]


def canonical_transcript_exons(gene_info: dict) -> tuple[str, list[list[int]], Optional[int], Optional[int]]:
    """Return (transcript_id, exons, cdsstart, cdsend) for the canonical-ish tx.

    Picks the transcript with the most exons as a proxy for the canonical one.
    Does not auto-choose a target exon — callers are responsible for asking
    the user to select one explicitly.
    """
    transcripts = gene_info.get("transcripts") or []
    if not transcripts:
        raise RuntimeError("mygene.info returned no transcripts for this gene.")
    tx = max(transcripts, key=lambda t: len(t["exons"]))
    return (
        tx.get("transcript") or "unknown",
        tx["exons"],
        tx.get("cdsstart"),
        tx.get("cdsend"),
    )


# ---------- variant lookup (dbSNP, ClinVar) ----------

# RefSeq chromosome accessions (e.g. "NC_000005") → UCSC chrom names.
# Also used by `oligomcp.variants._normalize_chrom`; kept here because
# dbSNP / ClinVar responses report chromosomes in RefSeq form.
_REFSEQ_CHROM_TABLE = {
    **{f"NC_{i:06d}": f"chr{i}" for i in range(1, 23)},
    "NC_000023": "chrX",
    "NC_000024": "chrY",
    "NC_012920": "chrM",
}


def _refseq_to_chrom(accession_base: str) -> Optional[str]:
    """Return a UCSC chrom name for a RefSeq accession base ('NC_000005')."""
    return _REFSEQ_CHROM_TABLE.get(accession_base)


def _variant_cache_path(key: str) -> Path:
    """Return an on-disk cache file for a variant lookup key."""
    VARIANT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    safe = "".join(c if (c.isalnum() or c in "._-") else "_" for c in key)
    return VARIANT_CACHE_DIR / f"{safe}.json"


def _assembly_aliases(assembly: str) -> tuple[str, ...]:
    """Accepted aliases for matching NCBI's assembly names.

    NCBI records report assembly as e.g. "GRCh38.p14" or "GRCh37.p13";
    our config speaks "hg38"/"hg19"/"grch38"/"grch37". Map to a tuple of
    case-insensitive prefixes for the caller to match against.
    """
    a = assembly.lower()
    if a in {"hg38", "grch38"}:
        return ("grch38",)
    if a in {"hg19", "grch37"}:
        return ("grch37",)
    return (a,)


def _fetch_json(url: str, *, timeout: int = 30) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "oligomcp/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def lookup_rsid_variant(rsid: str, *, assembly: str = "hg38") -> tuple[str, int, str, str]:
    """Resolve an rsID to ``(chrom, position, ref, alt)`` for the given assembly.

    Uses NCBI dbSNP's refsnp JSON (via eutils), which reports every
    `placement_with_allele` entry including assembly and allele-string.
    Responses are cached to disk to sidestep NCBI's ~3-req/s throttle
    during iteration.

    When a variant has multiple alt alleles, the first non-reference
    allele is returned and a warning is emitted; use an explicit dict
    input for disambiguation.
    """
    s = rsid.strip()
    if not s.lower().startswith("rs"):
        raise ValueError(f"rsID must start with 'rs': {rsid!r}")
    rs_num = s[2:]
    cache_path = _variant_cache_path(f"rsid_{rs_num}_{assembly.lower()}")
    if cache_path.exists():
        cached = json.loads(cache_path.read_text(encoding="utf-8"))
    else:
        url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rs_num}"
        cached = _fetch_json(url)
        cache_path.write_text(json.dumps(cached), encoding="utf-8")

    aliases = _assembly_aliases(assembly)
    placements = cached.get("primary_snapshot_data", {}).get("placements_with_allele", [])
    if not placements:
        raise RuntimeError(f"dbSNP returned no placements for {rsid!r}.")

    # Pick the placement whose assembly name matches our requested one.
    chosen = None
    for pl in placements:
        for alt_asm in (pl.get("placement_annot", {}).get("seq_id_traits_by_assembly") or []):
            asm_name = str(alt_asm.get("assembly_name", "")).lower()
            if any(asm_name.startswith(a) for a in aliases):
                chosen = pl
                break
        if chosen is not None:
            break
    if chosen is None:
        raise RuntimeError(
            f"dbSNP has no {assembly} placement for {rsid!r}. Available: "
            + ", ".join(sorted({
                str(a.get("assembly_name"))
                for pl in placements
                for a in (pl.get("placement_annot", {}).get("seq_id_traits_by_assembly") or [])
            }))
        )

    alleles = chosen.get("alleles", [])
    if not alleles:
        raise RuntimeError(f"dbSNP placement for {rsid!r} carries no alleles.")

    # Each allele carries a `spdi` record: {seq_id, position, deleted_sequence, inserted_sequence}.
    ref_alleles = [
        a for a in alleles
        if a.get("allele", {}).get("spdi", {}).get("deleted_sequence")
           == a.get("allele", {}).get("spdi", {}).get("inserted_sequence")
    ]
    alt_alleles = [a for a in alleles if a not in ref_alleles]
    if not alt_alleles:
        raise RuntimeError(f"dbSNP {rsid!r} carries only reference alleles.")

    chosen_alt = alt_alleles[0]
    if len(alt_alleles) > 1:
        import warnings as _warnings
        others = [a.get("allele", {}).get("spdi", {}).get("inserted_sequence") for a in alt_alleles[1:]]
        _warnings.warn(
            f"dbSNP {rsid!r} has multiple alt alleles; using the first "
            f"({chosen_alt.get('allele', {}).get('spdi', {}).get('inserted_sequence')!r}). "
            f"Others: {others!r}. Use an explicit variant dict to disambiguate.",
            stacklevel=2,
        )

    spdi = chosen_alt.get("allele", {}).get("spdi") or {}
    seq_id = str(spdi.get("seq_id", ""))
    # seq_id is a RefSeq accession like "NC_000005.10"; map to chrN.
    base = seq_id.split(".", 1)[0]
    chrom = _refseq_to_chrom(base)
    if chrom is None:
        raise RuntimeError(f"Unrecognized RefSeq accession in dbSNP placement: {seq_id!r}")

    # SPDI position is 0-based, VCF convention is 1-based.
    position = int(spdi.get("position", 0)) + 1
    ref = str(spdi.get("deleted_sequence") or "").upper()
    alt = str(spdi.get("inserted_sequence") or "").upper()
    return chrom, position, ref, alt


def lookup_clinvar_variant(accession: str, *, assembly: str = "hg38") -> tuple[str, int, str, str]:
    """Resolve a ClinVar VCV accession to ``(chrom, position, ref, alt)``.

    Backed by NCBI's esummary JSON endpoint; the response contains the
    VCF-style variation record under ``variation_set[0].variation_loc``.
    Disk-cached like rsID lookups.
    """
    s = accession.strip()
    if not s.upper().startswith("VCV"):
        raise ValueError(f"ClinVar accession must start with 'VCV': {accession!r}")
    vcv_num = s[3:].lstrip("0") or "0"
    cache_path = _variant_cache_path(f"clinvar_{vcv_num}_{assembly.lower()}")
    if cache_path.exists():
        cached = json.loads(cache_path.read_text(encoding="utf-8"))
    else:
        url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=clinvar&id={vcv_num}&retmode=json"
        )
        cached = _fetch_json(url)
        cache_path.write_text(json.dumps(cached), encoding="utf-8")

    result = cached.get("result") or {}
    record = result.get(str(vcv_num)) or {}
    var_set = record.get("variation_set") or []
    if not var_set:
        raise RuntimeError(f"ClinVar {accession!r} has no variation_set.")
    locs = var_set[0].get("variation_loc") or []
    aliases = _assembly_aliases(assembly)
    chosen = None
    for loc in locs:
        name = str(loc.get("assembly_name", "")).lower()
        if any(name.startswith(a) for a in aliases):
            chosen = loc
            break
    if chosen is None:
        raise RuntimeError(
            f"ClinVar {accession!r} has no {assembly} placement. Available: "
            + ", ".join({str(l.get("assembly_name")) for l in locs})
        )
    chrom_raw = str(chosen.get("chr", "")).strip()
    if not chrom_raw:
        raise RuntimeError(f"ClinVar {accession!r} placement is missing chromosome.")
    chrom = chrom_raw if chrom_raw.startswith("chr") else f"chr{chrom_raw}"
    position = int(chosen.get("start") or chosen.get("display_start") or 0)
    ref = str(chosen.get("ref") or chosen.get("reference_allele_vcf") or "").upper()
    alt = str(chosen.get("alt") or chosen.get("alternate_allele_vcf") or "").upper()
    if not ref and not alt:
        raise RuntimeError(f"ClinVar {accession!r} placement has no ref/alt alleles.")
    return chrom, position, ref, alt


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
    req = urllib.request.Request(url, headers={"User-Agent": "oligomcp/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        raise RuntimeError(
            f"Failed to fetch sequence from UCSC API for "
            f"{chrom}:{start}-{end} ({assembly}): {e}\n"
            "If you are offline, run `oligomcp fetch-genome` to download "
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

_BUNDLED_WEIGHTS_DIR = Path(__file__).parent / "_spliceai_weights"


def default_spliceai_weights_dir() -> Path:
    return SPLICEAI_DIR


def _bundled_weights_complete() -> bool:
    """True when all 5 MANE-10000nt `.pt` files ship inside the package."""
    return _BUNDLED_WEIGHTS_DIR.is_dir() and all(
        (_BUNDLED_WEIGHTS_DIR / f).exists() for f in _SPLICEAI_FILES
    )


def ensure_spliceai_weights(
    cache_dir: Optional[Path] = None, *, verbose: bool = True
) -> Path:
    """Resolve the MANE-10000nt SpliceAI weights location.

    Resolution order:
      1. Explicit `cache_dir` if complete.
      2. Bundled copy inside the installed package
         (`oligomcp/_spliceai_weights/`) — the default, so SpliceAI runs
         offline with zero configuration.
      3. User cache at `~/.oligomcp/spliceai/mane_10000nt/`, downloading
         any missing files from the upstream FTP mirror as a fallback.

    All 5 MANE-10000nt `.pt` files ship with the package (~14 MB total);
    whether `setup_spliceai()` loads 1 of them (default) or all 5 is
    controlled separately by `OLIGOMCP_SPLICEAI_N_MODELS`.
    """
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        missing = [f for f in _SPLICEAI_FILES if not (cache_dir / f).exists()]
        if not missing:
            if verbose:
                print(f"Using cached SpliceAI weights: {cache_dir}")
            return cache_dir
    else:
        if _bundled_weights_complete():
            if verbose:
                print(f"Using bundled SpliceAI weights: {_BUNDLED_WEIGHTS_DIR}")
            return _BUNDLED_WEIGHTS_DIR
        cache_dir = SPLICEAI_DIR
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
