# OligoMCP

Predict antisense oligonucleotide (ASO) efficacy with **AlphaGenome** and
**SpliceAI**, compare to experimental RT-PCR data, and visualize on the
**UCSC Genome Browser** — exposed as a Claude **MCP server** so you can
drive the whole pipeline from natural language in Claude Code / Desktop.

Three usage modes (all run the same workflow):

| Mode | Command | When to use |
|---|---|---|
| CLI | `oligomcp run --config foo.json` | Scripting, batch runs |
| Library | `from oligomcp import run_workflow` | Embed in notebooks / pipelines |
| **MCP** | `claude mcp add oligomcp -- oligomcp-mcp` | Interactive, natural-language flows |

## Install

```bash
git clone <this-repo>
cd <repo>
pip install -e .
oligomcp init            # one-time: downloads/caches what's needed, prompts for AlphaGenome key
```

Works on Linux, macOS, and Windows. The vendored pure-torch SpliceAI class
(bundled weights ship inside the package) means no openspliceai / pysam /
mappy C-extension dependencies. `oligomcp init` downloads the GRCh38 FASTA
on demand (optional — the tool falls back to the UCSC REST API for small
queries) and stores your AlphaGenome API key at
`~/.oligomcp/credentials.json` (mode 0600).

## Quick start (MCP)

After installing, register OligoMCP with whichever MCP-capable client you
use. The server speaks the standard **stdio** transport, so any host that
can launch a subprocess and speak MCP over stdin/stdout works — you're
configuring different clients, not different servers.

Then in any connected chat, try:

> *Analyze ASO candidates targeting SETD5 exon 11.*

Claude / the model will call `list_gene_exons("SETD5")`, ask you which
exon to target, then run `predict_aso_efficacy_inline` and return scored
candidates, a compact BED track, and a pre-loaded UCSC Genome Browser URL.

### Claude Code (CLI)

```bash
claude mcp add oligomcp -- oligomcp-mcp
```

Restart the session; the three tools appear in the toolbelt.

### Claude Desktop

Edit the desktop config file:

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

Omit the `env` block if you ran `oligomcp set-api-key` — the server reads
`~/.oligomcp/credentials.json` on startup. If `oligomcp-mcp` isn't on your
shell PATH (common with venvs / conda envs), replace `"command"` with the
absolute path, e.g. `"/home/you/.conda/envs/myenv/bin/oligomcp-mcp"`. Find
it with `which oligomcp-mcp` (macOS/Linux) or `where oligomcp-mcp`
(Windows). Restart Claude Desktop after saving.

### ChatGPT Desktop (OpenAI)

Recent ChatGPT Desktop versions support custom MCP connectors. Open
**Settings → Connectors → Add custom connector** and fill in:

- **Name:** `oligomcp`
- **Command:** absolute path to `oligomcp-mcp` (use `which oligomcp-mcp`)
- **Environment variables:** add `ALPHAGENOME_API_KEY` if you want
  AlphaGenome scoring without passing a key per request

The UI evolves; if those exact fields aren't present, look for an option
to add a local command / stdio MCP server and reference the `oligomcp-mcp`
binary. Consult OpenAI's current connector docs if in doubt.

### Other MCP-capable clients (Cursor, Cline, Continue, Zed, Goose, …)

Any client that supports launching a **local stdio MCP server** can run
OligoMCP. The configuration usually takes the same shape across
ecosystems — a JSON block with a `command` and optional `env` fields:

```json
{
  "oligomcp": {
    "command": "oligomcp-mcp",
    "env": { "ALPHAGENOME_API_KEY": "your-key-here" }
  }
}
```

For example:

- **Cursor** — `~/.cursor/mcp.json` (project scope: `.cursor/mcp.json`)
- **Cline** (VS Code) — settings panel → MCP servers
- **Continue** (VS Code / JetBrains) — `~/.continue/config.json` under `mcpServers`
- **Zed editor** — `settings.json` under `context_servers`
- **Goose** (Block) — `goose configure` and select "Add extension" → "Command-line extension"

Environment variables the server reads on startup:

| Variable | Default | Purpose |
|---|---|---|
| `ALPHAGENOME_API_KEY` | *(unset)* | AlphaGenome key; alternatively saved via `oligomcp set-api-key` or passed per-request |
| `OLIGOMCP_SPLICEAI_N_MODELS` | `1` | Set to `5` for full-ensemble SpliceAI scoring (slower, slightly better calibration) |
| `OLIGOMCP_PRELOAD_SPLICEAI` | `1` | Set to `0` to disable the background model preload (faster startup, slower first SpliceAI call) |

## Tools exposed by the MCP server

- **`list_gene_exons(gene_symbol, assembly="hg38")`** — resolves the
  canonical transcript via mygene.info, returns all exons with CDS/UTR
  annotations. No AlphaGenome key needed.
- **`predict_aso_efficacy_inline(gene_symbol, exon_intervals, ...)`** —
  full scoring pipeline driven by tool arguments; returns CSV + BED text
  + UCSC URL inline. Accepts `alphagenome_api_key` per request.
- **`predict_aso_efficacy(config_path, ...)`** — same workflow but driven
  from a JSON config on disk. Useful for reproducible local runs.

## Config schema (for CLI / file-based use)

```json
{
  "gene_symbol": "SETD5",
  "exon_intervals": [9440455, 9440698],
  "assembly": "hg38",
  "strand": "+",
  "ASO_length": 18,
  "aso_step": 5,
  "flank": [200, 200],
  "target_mode": "exclude",
  "requested_outputs": ["RNA_SEQ", "SPLICE_SITE_USAGE"],
  "ontology_terms": ["CL:0000127"],
  "experimental_data": "rt_pcr.csv"
}
```

Only `gene_symbol` and `exon_intervals` are strictly required; everything
else defaults sensibly. API keys are resolved from (in order):
`alphagenome_api_key` tool argument → `$ALPHAGENOME_API_KEY` env var →
`~/.oligomcp/credentials.json` → legacy `dna_api_key` field in the config
(discouraged). See `config/SETD5_e1.json` for a full example.

## How scoring works

For each ASO position, the target nucleotides are masked (N for
AlphaGenome, uniform 0.25 one-hot for SpliceAI) and the modified sequence
is re-scored by the model. The per-candidate efficacy follows the
`diff_mean` formula:

```
score = alt[exon].mean() / ref[body].mean() * alt[body].mean()
        - ref[exon].mean()
```

AlphaGenome evaluates this over a resized gene interval via
`model.predict_sequences`. SpliceAI runs the vendored MANE-10000nt
single-model (full 5-model ensemble available via
`OLIGOMCP_SPLICEAI_N_MODELS=5`) on a `(SL=5000) + (CL=10000)` window
centered on the exon.

Positive scores = ASO strengthens exon inclusion; negative = weakens /
promotes skipping. For `target_mode="exclude"` you want the most-negative
scores.

## Outputs

Files land in `<results_dir>/ASO/`:

- `<name>_ASO_scores.csv` — per-ASO scores, all sources
- `<name>_ASO_<source>.bed` / `_full.bed` — UCSC custom tracks (red gradient
  = weakening, blue = strengthening)
- `<name>_ASO_Measured.bed` — experimental RT-PCR track, if provided (green)
- `<name>_correlation.png` — per-exon regression plot, if experimental data
- `<name>_experimental_matched.csv` — aggregated predictions vs measured

The CLI auto-opens the UCSC Browser with compact tracks pre-loaded (pass
`--no-browser` to suppress). The MCP's inline tool returns the URL in its
response payload.

## Package layout

```
server.py                 Top-level FastMCP entrypoint (for remote deploys)
src/oligomcp/
  config.py               Config loading + OligoConfig dataclass
  core.py                 Sequence utils, ASO enumeration, experimental helpers
  predict.py              AlphaGenome + SpliceAI scoring (diff_mean)
  output.py               BED export, UCSC URL + browser auto-open, plotting
  resources.py            Credentials, UCSC sequence fetch, mygene.info, weights
  workflow.py             End-to-end pipeline orchestration
  cli.py                  CLI entry (oligomcp, oligomcp-mcp)
  mcp_server.py           FastMCP server with 3 tools
  _spliceai_model.py      Vendored SpliceAI class (openspliceai v0.0.5, MIT)
  _spliceai_weights/*.pt  Bundled MANE-10000nt weights (~14 MB)
```

## License

See `LICENSE`.
