# OligoMCP

Predict antisense oligonucleotide (ASO) efficacy with **AlphaGenome** +
**SpliceAI**, compare to experimental RT-PCR data, and drive the
whole pipeline from natural language in Claude (MCP) or the CLI.

## Install

```bash
git clone <this-repo>
cd <repo>
pip install -e .
oligomcp init            # prompts for AlphaGenome key, bundles SpliceAI weights
```

Linux / macOS / Windows. SpliceAI weights ship with the package;
genomic sequence is fetched on-demand from UCSC — no genome download.

## Use it from Claude (MCP)

**Register once** — Claude Code:

```bash
claude mcp add oligomcp -- oligomcp-mcp
```

Other clients (Claude Desktop, Cursor, Cline, Continue, Zed, Goose,
ChatGPT Desktop) — add this block to their MCP config file:

```json
{ "mcpServers": { "oligomcp": {
  "command": "oligomcp-mcp",
  "env": { "ALPHAGENOME_API_KEY": "your-key-here" }
} } }
```

The `env` line is optional if you already ran `oligomcp init`. If
`oligomcp-mcp` isn't on PATH, replace with the full path from
`which oligomcp-mcp`.

**Then prompt naturally:**

> **"Design ASOs to skip SETD5 exon 11. Show the top 10 by SpliceAI."**

> **"Score ASOs for SMN2 exon 7, restricted to motor-neuron tracks."**
> (Claude calls `search_ontology_terms("motor neuron")` → `CL:0000100`,
> passes it as `ontology_terms`.)

> **"Run `config/SETD5_e1.json` and correlate with the experimental CSV."**

> **"Just SpliceAI — I don't have an AlphaGenome key."**
> (Claude adds `skip_alphagenome=true`.)

> **"Which ontology terms cover liver RNA-seq?"**

Claude chains `list_gene_exons` → `search_ontology_terms` →
`predict_aso_efficacy{_inline}` as needed and returns per-ASO scores,
BED file paths, and correlation stats.

## Use it from the CLI

```bash
# Full AlphaGenome + SpliceAI run from a config
oligomcp run --config config/SETD5_e1.json -v

# SpliceAI-only (no AlphaGenome key needed)
oligomcp run --config config/SETD5_e1.json --skip-alphagenome
```

Minimum config — just the gene + exon:

```json
{ "gene_symbol": "SETD5", "exon_intervals": [9429826, 9430051] }
```

Everything else defaults. See `config/SETD5_e1.json` for a realistic
example (custom flank, step, ontology terms, experimental CSV).

**Picking `ontology_terms`**: all 704 AlphaGenome CURIEs are
pre-snapshotted at `data/alphagenome_ontology_terms.tsv` — `grep` it
or call `search_ontology_terms` from the MCP. Empty list = average
across all tracks.

## Outputs

Land under `results/<gene_symbol>/`. All BED9 — usable in UCSC Genome
Browser (upload at <https://genome.ucsc.edu/cgi-bin/hgCustom>), IGV,
JBrowse, etc.

| File | Contents |
|---|---|
| `<name>_ASO_scores.csv` | Per-ASO scores across all sources |
| `<name>_ASO_<source>.bed` / `_full.bed` | Custom tracks — red = inclusion-enhancing, blue = skip-promoting |
| `<name>_ASO_Measured.bed` | Green-gradient RT-PCR track (if `experimental_data` set) |
| `<name>_correlation.png` | Per-exon regression with percentile-rank x-axis |
| `<name>_experimental_matched.csv` | Aggregated predictions ↔ measured |

**Sign convention**: positive score = strengthens exon inclusion;
negative = promotes skipping. For `target_mode="exclude"`, pick the
most-negative hits.

## MCP tools (4)

| Tool | Purpose |
|---|---|
| `list_gene_exons(gene_symbol, assembly="hg38")` | Exons of the canonical transcript. No API key. |
| `search_ontology_terms(query, output_type=None, limit=20)` | Substring search over the 704-CURIE ontology snapshot. No API key. |
| `predict_aso_efficacy(config_path, ...)` | File-based scoring. Returns `needs_info` if design fields missing. |
| `predict_aso_efficacy_inline(gene_symbol, exon_intervals, ...)` | Arg-driven scoring, up to 300 candidates/request. |

## Environment variables (optional)

| Variable | Default | Purpose |
|---|---|---|
| `ALPHAGENOME_API_KEY` | *(unset)* | AlphaGenome key (or use `oligomcp init` / `set-api-key`) |
| `OLIGOMCP_SPLICEAI_N_MODELS` | `1` | Set to `5` for full-ensemble SpliceAI |
| `OLIGOMCP_PRELOAD_SPLICEAI` | `1` | Set to `0` to skip background model preload |

## License

See `LICENSE`.
