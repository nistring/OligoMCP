# OligoMCP

Predict antisense oligonucleotide (ASO) efficacy with **AlphaGenome** and
**SpliceAI**, compare to experimental RT-PCR data when available, and
drive the whole pipeline from natural language via a Claude **MCP server**
or from the CLI.

| Mode | Command | When to use |
|---|---|---|
| **MCP** | `claude mcp add oligomcp -- oligomcp-mcp` | Interactive, natural-language flows in Claude Code / Desktop / Cursor / … |
| **CLI** | `oligomcp run --config foo.json` | Scripting, batch runs, reproducible experiments |
| **Library** | `from oligomcp import run_workflow` | Embed in notebooks / pipelines |

## Install

```bash
git clone <this-repo>
cd <repo>
pip install -e .
oligomcp init            # one-time: SpliceAI weights + AlphaGenome key prompt
```

Works on Linux, macOS, and Windows. SpliceAI weights (~14 MB, MANE-10000nt)
are bundled with the package — no `openspliceai` / `pysam` / `mappy` C
dependencies. Genomic sequence is fetched on demand from the UCSC REST API,
so there's no multi-GB genome download required for typical runs.

## Quick start — CLI

```bash
# 1. One-time setup (prompts for your AlphaGenome key, or skip and add later)
oligomcp init

# 2. Run on an example config (scores 34 ASOs across SETD5 exon 11)
oligomcp run --config config/SETD5_e1.json -v

# 3. Drop AlphaGenome if you only want SpliceAI (no key needed)
oligomcp run --config config/SETD5_e1.json --skip-alphagenome -v
```

Outputs land in `results/<gene_symbol>/` — one CSV of per-ASO scores, a
pair of BED files per score source (compact top-20 + full), and a
correlation PNG if `experimental_data` is provided.

Minimum config — just `gene_symbol` and `exon_intervals`:

```json
{ "gene_symbol": "SETD5", "exon_intervals": [9429826, 9430051] }
```

Full schema (every optional field):

```json
{
  "gene_symbol":      "SETD5",
  "exon_intervals":   [9429826, 9430051],
  "assembly":         "hg38",
  "strand":           "+",
  "ASO_length":       18,
  "aso_step":         1,
  "flank":            [200, 200],
  "target_mode":      "exclude",
  "requested_outputs":["RNA_SEQ", "SPLICE_SITE_USAGE"],
  "ontology_terms":   ["CL:0000127"],
  "experimental_data":"ASOseq_incl.csv"
}
```

See `config/SETD5_e1.json` and `config/SMN2.json` for real examples.

## Quick start — MCP

Register the stdio server with your MCP-capable client, then prompt in
natural language. Example:

> *Analyze ASO candidates targeting SETD5 exon 11 (chr3:9429826-9430051).*

Claude calls `list_gene_exons("SETD5")` first to confirm the exon, then
`predict_aso_efficacy` / `predict_aso_efficacy_inline` and surfaces
per-ASO scores, BED file paths, and correlation stats inline.

### Claude Code (CLI)

```bash
claude mcp add oligomcp -- oligomcp-mcp
```

Restart the session; the three tools appear in the toolbelt. If
`oligomcp-mcp` isn't on the shell PATH (common with venvs / conda envs),
pass the absolute path: `claude mcp add oligomcp -- /path/to/oligomcp-mcp`.

### Claude Desktop

Edit the desktop config:

| OS | Path |
|----|------|
| Linux | `~/.config/Claude/claude_desktop_config.json` |
| macOS | `~/Library/Application Support/Claude/claude_desktop_config.json` |
| Windows | `%APPDATA%\Claude\claude_desktop_config.json` |

```json
{
  "mcpServers": {
    "oligomcp": {
      "command": "oligomcp-mcp",
      "env": { "ALPHAGENOME_API_KEY": "your-key-here" }
    }
  }
}
```

Omit the `env` block if you already ran `oligomcp init` or `oligomcp
set-api-key` — the server reads `~/.oligomcp/credentials.json` on
startup. If the `command` isn't found, use the absolute path (find it
with `which oligomcp-mcp`). Restart Claude Desktop.

### Other MCP clients (Cursor, Cline, Continue, Zed, Goose, ChatGPT Desktop, …)

Any client that can launch a **local stdio MCP server** works. The
block is nearly identical across clients — a `command` and optional
`env`:

```json
{
  "oligomcp": {
    "command": "oligomcp-mcp",
    "env": { "ALPHAGENOME_API_KEY": "your-key-here" }
  }
}
```

Per-client paths:

- **Cursor** — `~/.cursor/mcp.json` (project: `.cursor/mcp.json`)
- **Cline** (VS Code) — settings panel → MCP servers
- **Continue** (VS Code / JetBrains) — `~/.continue/config.json` under `mcpServers`
- **Zed** — `settings.json` under `context_servers`
- **Goose** — `goose configure` → Add extension → Command-line extension
- **ChatGPT Desktop** — Settings → Connectors → Add custom connector

## MCP tools

- `list_gene_exons(gene_symbol, assembly="hg38")` — canonical-transcript
  exons via mygene.info, annotated CDS / UTR / first / last. No
  AlphaGenome key needed.
- `predict_aso_efficacy(config_path, ...)` — file-based scoring.
  Returns `status="needs_info"` if design-critical fields
  (`ASO_length`, `aso_step`, `flank`, `target_mode`) are missing so
  Claude can ask the user before silently defaulting. Pass
  `confirm_defaults=True` to accept defaults.
- `predict_aso_efficacy_inline(gene_symbol, exon_intervals, ...)` —
  arg-driven scoring (no config file needed), capped at 300 ASO
  candidates per request. Returns CSV + BED content inline.

## Environment variables

| Variable | Default | Purpose |
|---|---|---|
| `ALPHAGENOME_API_KEY` | *(unset)* | AlphaGenome key; alternative to `oligomcp set-api-key` |
| `OLIGOMCP_SPLICEAI_N_MODELS` | `1` | Set to `5` for full-ensemble SpliceAI (slower, marginal calibration gain) |
| `OLIGOMCP_PRELOAD_SPLICEAI` | `1` | Set to `0` to skip background model preload (faster startup, slower first call) |

## How scoring works

For each ASO position, the target nucleotides are masked (N for
AlphaGenome, uniform 0.25 one-hot for SpliceAI). The perturbed sequence
is re-scored and compared to reference via `diff_mean`:

```
score = alt[target].mean() / ref[body].mean() * alt[body].mean()
        - ref[target].mean()
```

Three score sources are emitted, each with a different `target` slice:

- **AlphaGenome RNA_SEQ** — `target` = whole exon interval (coverage
  track; signal is distributed).
- **AlphaGenome SPLICE_SITE_USAGE** — `target` = just the two splice
  junctions (acceptor at `exon_start`, donor at `exon_end − 1`), since
  the track is peaked there.
- **SpliceAI** — same junction-only target. SpliceAI's softmax output
  is ~0 everywhere except at canonical splice sites; averaging over the
  full exon would dilute the signal with hundreds of near-zero bases.

Positive scores = strengthens inclusion; negative = promotes skipping.
For `target_mode="exclude"` you want the most-negative scores.

## Outputs

Everything lands in `<results_dir>/<gene_symbol>/`:

- `<name>_ASO_scores.csv` — per-ASO scores, all sources
- `<name>_ASO_<source>.bed` / `_full.bed` — UCSC custom tracks (red
  gradient = inclusion-enhancing, blue = skip-promoting, saturation
  normalized within each track)
- `<name>_ASO_Measured.bed` — RT-PCR track (green gradient), if
  `experimental_data` is provided
- `<name>_correlation.png` — per-exon regression plot with
  percentile-rank x-axis so all three tracks overlay cleanly
- `<name>_experimental_matched.csv` — aggregated predictions ↔ measured

To view in UCSC Genome Browser: open
<https://genome.ucsc.edu/cgi-bin/hgCustom>, upload any of the BED files
(or paste the BED text). The same files work directly in IGV, JBrowse,
and most other genome viewers.

## Package layout

```
src/oligomcp/
  config.py               OligoConfig dataclass + JSON loader
  core.py                 Sequence utils, ASO enumeration, experimental helpers
  predict.py              AlphaGenome + SpliceAI scoring (diff_mean)
  output.py               BED export + per-exon correlation plot
  resources.py            Credentials, UCSC sequence fetch, mygene.info, weights
  workflow.py             End-to-end pipeline orchestration
  cli.py                  CLI entrypoint (`oligomcp run`, `init`, …)
  mcp_server.py           FastMCP server with 3 tools
  _spliceai_model.py      Vendored SpliceAI class (openspliceai v0.0.5, MIT)
  _spliceai_weights/*.pt  Bundled MANE-10000nt weights (~14 MB)
```

## License

See `LICENSE`.
