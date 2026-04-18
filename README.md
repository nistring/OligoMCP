# OligoClaude

Predict antisense oligonucleotide (ASO) efficacy using **AlphaGenome** and **SpliceAI**, compare to experimental RT-PCR data, and visualize results on the **UCSC Genome Browser** — all from a single JSON config.

Exposed four ways:
1. CLI (`oligoclaude run --config …`)
2. Python library (`from oligoclaude import run_workflow`)
3. **Local MCP server** so Claude Desktop / Claude Code can call it as a connector
4. **Hosted HTTP MCP server** on FastMCP Cloud / Prefect Horizon, so anyone
   can use it from a Claude client without cloning the repo (see the
   "Remote HTTP deploy" section below).

## Install

```bash
git clone <this-repo>
cd OligoClaude
pip install -e .
```

Works on Linux, macOS, and Windows with no extra steps. OligoClaude vendors
the pure-torch SpliceAI class from OpenSpliceAI (see
`src/oligoclaude/_spliceai_model.py` for MIT attribution), so it doesn't
depend on the `openspliceai` package or its C-extension transitive deps
(`mappy`, `pysam`). Only the pretrained weights are fetched at runtime.

For unit tests: `pip install -e .[dev]`

### One-time setup

Run once after installing:

```bash
oligoclaude init
```

This downloads the SpliceAI MANE-10000nt weights (~3 MB × 5) to
`~/.oligoclaude/spliceai/` and prompts (hidden) for your AlphaGenome API key
(saved to `~/.oligoclaude/credentials.json`, mode 0600). Both steps are
idempotent — re-running `init` skips what's already configured.

Non-interactive / CI / Docker:

```bash
ALPHAGENOME_API_KEY=sk-... oligoclaude init -y           # use env var, skip prompt
oligoclaude init --skip-api-key                           # just download weights
oligoclaude init --skip-spliceai                          # AlphaGenome-only install
```

**No genome download required.** OligoClaude fetches only the needed genomic
region on-the-fly from the UCSC REST API (typically <1 second for 500 kb vs
downloading the full ~3 GB genome). If you prefer offline use or run many
analyses, you can optionally cache the full genome locally:

```bash
oligoclaude fetch-genome                        # optional: downloads GRCh38 FASTA (~3 GB)
```

When a local FASTA is available (via `fetch-genome` or a `fasta_path` in the
config), OligoClaude uses it automatically via pyfaidx. Otherwise it
transparently falls back to the UCSC API.

Individual-step commands are also available if you prefer: `oligoclaude
set-api-key`, `oligoclaude fetch-spliceai-weights`, `oligoclaude fetch-genome`.
You can also set the API key via the `ALPHAGENOME_API_KEY` environment
variable instead of the credentials file.

## Config

Create a JSON config (see `config/SETD5_e1.json` for a full example):

```json
{
  "gene_symbol": "SETD5",
  "assembly": "hg38",
  "gtf_url": "https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather",
  "results_dir": "results",
  "data_dir": "data",
  "track_filter": "polyA plus RNA-seq",
  "ontology_terms": ["CL:0000127", "CL:0002319"],
  "requested_outputs": ["RNA_SEQ", "SPLICE_SITE_USAGE"],
  "exon_intervals": [9429826, 9430051],
  "flank": [200, 200],
  "strand": "+",
  "ASO_length": 18,
  "aso_step": 1,
  "experimental_data": "ASOseq_incl.csv",
  "target_mode": "exclude",
  "spliceai_batch": 12
}
```

Required fields: `gene_symbol`, `assembly`, `gtf_url`, `exon_intervals`, `strand`, `ASO_length`, `flank`, `ontology_terms`, `requested_outputs`.

**Do NOT put `dna_api_key` in the config file.** Use `oligoclaude set-api-key`
or the `ALPHAGENOME_API_KEY` environment variable. A legacy `dna_api_key`
field is still honored for backward compatibility, but triggers a
`DeprecationWarning` because it gets committed to version control.

**`fasta_path` is optional.** If omitted, OligoClaude fetches only the needed
region from the UCSC REST API on-the-fly. If a local genome is cached
(`oligoclaude fetch-genome`), it is used automatically via pyfaidx.

Optional fields:
- `experimental_data`: path to a CSV with columns `ASO_ID`, `ASO sequence`, `Measured (RT-PCR)`, optionally `Region (Exon)`. If present, OligoClaude scores those ASOs directly and produces a correlation plot; otherwise it runs a sliding window across `[exon_start - flank[0], exon_end + flank[1]]`.
- `aso_step`: sliding window step size (default 1).
- `target_mode`: `"exclude"` (default) or `"include"`. Only affects interpretation, not raw scores.
- `spliceai_batch`, `spliceai_threads`: SpliceAI CPU batching knobs.

Relative paths inside the JSON are resolved against the config file's directory.

## CLI usage

```bash
# Full workflow (AlphaGenome + SpliceAI + correlation plot + UCSC browser)
oligoclaude run --config config/SETD5_e1.json -v

# Skip SpliceAI (faster; AlphaGenome only)
oligoclaude run --config config/SETD5_e1.json --skip-spliceai

# Skip AlphaGenome (no API key needed)
oligoclaude run --config config/SETD5_e1.json --skip-alphagenome

# Suppress auto-opening the UCSC Genome Browser
oligoclaude run --config config/SETD5_e1.json --no-browser
```

Outputs land in `<results_dir>/ASO/`:
- `<config_name>_ASO_scores.csv` — raw scores per ASO, all sources
- `<config_name>_ASO_<source>.bed` / `_full.bed` — UCSC custom track files (predicted)
- `<config_name>_ASO_Measured.bed` — experimental RT-PCR track (if experimental data provided)
- `<config_name>_correlation.png` — per-exon + combined correlation plot (if experimental data)
- `<config_name>_experimental_matched.csv` — aggregated predictions vs measured

## UCSC Genome Browser

By default, `oligoclaude run` **automatically opens the UCSC Genome Browser**
in your default browser with all compact BED tracks (predicted + experimental)
pre-loaded. The view is centered on the target exon with the configured flank.

The predicted tracks use a **red/blue gradient** (red = exon weakening,
blue = exon strengthening) and the experimental track uses a **green gradient**
(darker green = higher measured RT-PCR inclusion).

If you prefer to upload manually, pass `--no-browser` and use the printed
file paths:

1. Go to `https://genome.ucsc.edu/cgi-bin/hgCustom?db=hg38`
2. Click **Choose File** and pick one of the generated `*.bed` files
3. Click **Submit**

## Use from Claude (MCP connector)

OligoClaude exposes three MCP tools that Claude Code and Claude Desktop can
call as a first-class connector:

- **`list_gene_exons(gene_symbol, assembly)`** — discover exons of a gene's
  canonical transcript via mygene.info. Use this first when the user only
  supplies a gene symbol so you can show them the exons and ask which to
  target. No AlphaGenome key required.
- **`predict_aso_efficacy_inline(gene_symbol, exon_intervals, ...)`** —
  scoring API that takes config fields as tool arguments and returns CSV
  and BED content inline (no files on disk for the caller). Designed for
  **remote HTTP deployments** (FastMCP Cloud / Prefect Horizon) but works
  locally too.
- **`predict_aso_efficacy(config_path, ...)`** — the original file-based
  tool. Requires a JSON config on the server's filesystem, so only useful
  for **local stdio** deployments.

### Local stdio: register in Claude Code / Desktop

```bash
pip install -e /path/to/OligoClaude
claude mcp add oligoclaude -- oligoclaude-mcp
```

Or Claude Desktop `~/.config/Claude/claude_desktop_config.json` (Linux),
`~/Library/Application Support/Claude/claude_desktop_config.json` (macOS), or
`%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "oligoclaude": {
      "command": "oligoclaude-mcp",
      "env": { "ALPHAGENOME_API_KEY": "your-alphagenome-key-here" }
    }
  }
}
```

Alternatively, omit the `env` block and run `oligoclaude set-api-key` once —
the server reads `~/.oligoclaude/credentials.json` on startup.

Restart Claude Desktop. Then in any conversation you can say:

> *"Analyze ASO candidates for the SETD5 gene"*

Claude will call `list_gene_exons("SETD5")`, present the exons, ask you
which to target, and then run `predict_aso_efficacy_inline` (or the
file-based tool if you point it at a config).

### Remote HTTP deploy on Prefect Horizon (FastMCP Cloud)

Anyone can use the hosted server without cloning the repo. Deploy once:

1. Fork this repo (or leave it as-is if you own it) — the deploy needs a
   GitHub repo to pull from.
2. Sign in at [horizon.prefect.io](https://horizon.prefect.io) with your
   GitHub account.
3. Create a new server, select the repo, and use entrypoint
   **`server.py:mcp`**. Horizon auto-installs the dependencies from
   `pyproject.toml`.
4. Click deploy. Horizon publishes the server at
   `https://<your-name>.fastmcp.app/mcp` in ~60 seconds and
   auto-redeploys on every push to `main`.

Horizon also supports PR-preview URLs out of the box. No AlphaGenome key
is baked into the deployment — end users pass their own key as the
`alphagenome_api_key` argument to each `predict_aso_efficacy_inline`
call (in-memory only, not persisted server-side). SpliceAI weights are
downloaded into a per-container cache on first use (~14 MB).

To register the hosted URL in Claude Code:

```bash
claude mcp add --transport http oligoclaude https://<your-name>.fastmcp.app/mcp
```

Remote-safety knobs:
- Request cap: `predict_aso_efficacy_inline` is capped at 300 ASO
  candidates per call so CPU-bound SpliceAI inference stays under
  ~60 seconds. The tool returns a `too_many_candidates` status with a
  suggested `aso_step` when you exceed it.
- Private install: Horizon offers an "Authentication" toggle so only
  your org can connect. Leave it off for a public demo.

### Tool signatures

```python
list_gene_exons(gene_symbol: str, assembly: str = "hg38") -> dict

predict_aso_efficacy_inline(
    gene_symbol: str,
    exon_intervals: list[int],        # [start, end]
    assembly: str = "hg38",
    strand: str = "+",
    aso_length: int = 18,
    aso_step: int = 5,                # larger = fewer candidates
    flank: list[int] = [200, 200],
    target_mode: str = "exclude",
    skip_alphagenome: bool = True,
    skip_spliceai: bool = False,
    alphagenome_api_key: str | None = None,   # per-request, in-memory only
    samples_max: int = 20,
    requested_outputs: list[str] | None = None,
    ontology_terms: list[str] | None = None,
    track_filter: str = "polyA plus RNA-seq",
) -> dict  # status, n_candidates, scores, top_candidates, bed_tracks, ucsc_url

predict_aso_efficacy(
    config_path: str,
    skip_alphagenome: bool = False,
    skip_spliceai: bool = False,
    samples_max: int = 20,
) -> dict  # local-stdio only; returns server-side file paths
```

## How scoring works

**AlphaGenome** — for each ASO position, the target nucleotides are replaced with `N`s and the resulting sequence is predicted via `model.predict_sequences`. Per-output-type efficacy follows the `diff_mean` formula from AlphaGenome_ASO/aso.ipynb:

```
score = alt[exon].mean() / ref[gene_body].mean() * alt[gene_body].mean() - ref[exon].mean()
```

**SpliceAI** — the reference sequence (`SL=5000` centered on the exon, plus `CL_max=10000` flanking context) is one-hot encoded once; per-ASO variants are built by `clone() + slice-assign` to set `[0.25, 0.25, 0.25, 0.25]` at the ASO position. The 5-model OpenSpliceAI MANE-10000nt ensemble runs on CPU with configurable batching. The efficacy score uses the same `diff_mean` formula applied to per-position donor+acceptor probability sums.

Both scores are positive when the ASO weakens exon usage (exon skipping) and negative when it strengthens it.

## Package structure

```
server.py         Top-level entrypoint for Horizon / FastMCP Cloud
                  (re-exports the `mcp` object via absolute imports)
src/oligoclaude/
  config.py       JSON config loading + OligoConfig dataclass
  core.py         Sequence utils, ASO enumeration, experimental-data helpers
  predict.py      AlphaGenome + SpliceAI scoring (diff_mean formula)
  output.py       BED export, UCSC URL + browser auto-open, correlation plots
  resources.py    API credentials, UCSC sequence fetch, SpliceAI weight cache,
                  mygene.info gene / exon lookup
  workflow.py     End-to-end pipeline orchestration
  cli.py          CLI entry point (oligoclaude, oligoclaude-mcp)
  mcp_server.py   FastMCP server with three tools (list_gene_exons,
                  predict_aso_efficacy_inline, predict_aso_efficacy)
```

## Library usage

```python
from oligoclaude import run_workflow

result = run_workflow(
    "config/SETD5_e1.json",
    verbose=True,
    open_browser=False,    # suppress UCSC auto-open
)
print(result.scores_csv)
print(result.stats)
print(result.ucsc_url)     # URL to open manually if needed
```

## License

See `LICENSE`.
