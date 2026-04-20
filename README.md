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
pre-snapshotted at `data/alphagenome_ontology_terms.tsv` (one row per
CURIE, including each term's available `track_filters` like `polyA
plus RNA-seq`) — `grep` it or call `search_ontology_terms` from the
MCP. Empty list = average across all tracks. Regenerate the snapshot
with `oligomcp fetch-ontology-terms` (needs AlphaGenome key).

## Patient baseline (apply variants before ASO design)

Add a `variants` list to your config to **edit mutations into
`ref_seq` before ASO enumeration** — the pipeline then designs against
the patient's actual mRNA context, not wild-type. Useful when a known
SNP would disrupt ASO affinity, or when you need to overlap a
disease-causing variant.

```json
{
  "gene_symbol": "SMN2",
  "exon_intervals": [70070640, 70070751],
  "variants": [
    "chr5:70070740:G:A",
    {"id": "demo_3bp_del", "chrom": "chr5", "position": 70070700,
     "ref": "AAG", "alt": ""}
  ]
}
```

The list is heterogeneous — every entry below can appear in the same
`variants` array, and all are applied together to form one combined
patient baseline (overlaps are rejected):

| Form | Example |
|---|---|
| VCF-style | `"chr5:70070740:G:A"` or `"5-70070740-G-A"` |
| HGVS genomic | `"chr5:g.70070740G>A"`, `"chr5:g.70070700_70070702del"`, `"chr5:g.70070700_70070701insAT"` |
| HGVS coding | `"c.840C>T"` (resolved against the gene's canonical transcript + strand) |
| rsID | `"rs1800112"` (dbSNP lookup, disk-cached under `~/.oligomcp/variant_cache/`) |
| ClinVar | `"VCV000000001"` |
| Explicit dict | `{"chrom":"chr5","position":70070700,"ref":"AAG","alt":"","id":"my_del"}` |

**Indels are supported** — AlphaGenome's interval is kept at its fixed
width by padding with downstream reference (for deletions) or trimming
the right edge (for insertions); SpliceAI's 15 kb window is rebuilt
the same way. The WT AlphaGenome baseline is **re-fetched on the
patient sequence** (`predict_sequence`) so `diff_mean_frac` compares
patient-WT vs patient-ASO — otherwise every ASO score would absorb the
variant's own splicing effect.

Outputs:
- `<name>_applied_variants.json` — run-reproducibility metadata (every
  variant, its genomic position, ref→alt, Δ length, offset inside the
  loaded window).
- BED records are mapped back to reference-genome coordinates via
  `VariantCoordMap.patient_to_ref`. ASOs whose binding site lands
  inside an inserted region have no reference position and are kept in
  the CSV but skipped in BED.

From the MCP tools, pass the same list via the `variants` argument to
`predict_aso_efficacy_inline`, or just run
`predict_aso_efficacy(config_path)` against a JSON config that
includes the field. Parser / ref-mismatch / overlap failures surface
as `status="variant_error"` with the offending record.

## Outputs

Land under `results/<gene_symbol>/<YYYYMMDD_HHMMSS>/` — each run gets
its own timestamped subdirectory, so reruns never overwrite prior
output. All BED9, usable in UCSC Genome Browser (upload at
<https://genome.ucsc.edu/cgi-bin/hgCustom>), IGV, JBrowse, etc.

| File | Contents |
|---|---|
| `<name>_scores.csv` | Per-ASO scores across all sources |
| `<name>_<source>.bed` | Custom tracks (one file per source, every candidate) — red = inclusion-enhancing, blue = skip-promoting |
| `<name>_Measured.bed` | Green-gradient RT-PCR track (if `experimental_data` set) |
| `<name>_correlation.png` | Per-exon regression with percentile-rank x-axis |
| `<name>_experimental_matched.csv` | Aggregated predictions ↔ measured |

**Sign convention**: positive score = strengthens exon inclusion;
negative = promotes skipping. For `target_mode="exclude"`, pick the
most-negative hits.

## MCP tools (4)

| Tool | Purpose |
|---|---|
| `list_gene_exons(gene_symbol, assembly="hg38")` | Exons of the canonical transcript. No API key. |
| `search_ontology_terms(query, output_type=None, track_filter=None, limit=20)` | Substring search over the 704-CURIE ontology snapshot, including each CURIE's available `track_filters` (e.g. "polyA plus RNA-seq") so you can confirm a filter exists before passing it. No API key. |
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
