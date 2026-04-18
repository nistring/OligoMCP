"""CLI entrypoint for OligoMCP."""
from __future__ import annotations

import argparse
import getpass
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="oligomcp",
        description="Predict ASO efficacy via AlphaGenome and SpliceAI.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="Run the full ASO prediction workflow.")
    r.add_argument("--config", type=Path, required=True)
    r.add_argument("--skip-alphagenome", action="store_true")
    r.add_argument("--skip-spliceai", action="store_true")
    r.add_argument(
        "--samples-max",
        type=int,
        default=20,
        help="Top/bottom ASOs to emit in compact BED (does not affect correlation).",
    )
    r.add_argument("--verbose", "-v", action="store_true")

    i = sub.add_parser(
        "init",
        help="One-time setup: download SpliceAI weights and save AlphaGenome key.",
    )
    i.add_argument(
        "--skip-spliceai",
        action="store_true",
        help="Do not download SpliceAI weights (use for AlphaGenome-only installs).",
    )
    i.add_argument(
        "--skip-api-key",
        action="store_true",
        help="Do not prompt for the AlphaGenome API key (set it later via set-api-key).",
    )
    i.add_argument(
        "--yes", "-y", action="store_true",
        help="Accept defaults non-interactively (useful for CI / Docker).",
    )

    k = sub.add_parser("set-api-key", help="Save AlphaGenome API key (mode 0600).")
    k.add_argument("key", nargs="?", help="Key value; prompted if omitted.")

    sub.add_parser("clear-api-key", help="Remove the saved AlphaGenome API key.")

    for name, helptxt in [
        ("fetch-genome", "Download GRCh38 FASTA to ~/.oligomcp/genomes/."),
        ("fetch-spliceai-weights", "Download the MANE-10000nt ensemble."),
    ]:
        f = sub.add_parser(name, help=helptxt)
        f.add_argument("--cache-dir", type=Path, default=None)

    o = sub.add_parser(
        "fetch-ontology-terms",
        help="Regenerate data/alphagenome_ontology_terms.tsv from the live AlphaGenome API.",
    )
    o.add_argument(
        "--data-dir", type=Path, default=None,
        help="Output directory (default: repo-root `data/`).",
    )

    return p


def _cmd_run(args) -> int:
    from .workflow import run_workflow

    result = run_workflow(
        config_path=args.config,
        skip_alphagenome=args.skip_alphagenome,
        skip_spliceai=args.skip_spliceai,
        samples_max=args.samples_max,
        verbose=args.verbose,
    )
    print()
    print(f"Candidates scored: {result.n_candidates}")
    print(f"Scores CSV:        {result.scores_csv}")
    if result.correlation_plot:
        print(f"Correlation plot:  {result.correlation_plot}")
    print("BED files:")
    for p in result.bed_files:
        print(f"  - {p}")
    return 0


def _cmd_set_api_key(args) -> int:
    from .resources import save_alphagenome_api_key

    key = args.key or getpass.getpass("AlphaGenome API key (hidden): ")
    path = save_alphagenome_api_key(key)
    print(f"Saved AlphaGenome API key to {path} (mode 0600).")
    return 0


def _cmd_clear_api_key(args) -> int:
    from .resources import CRED_PATH, clear_alphagenome_api_key

    if clear_alphagenome_api_key():
        print(f"Removed AlphaGenome API key from {CRED_PATH}.")
    else:
        print("No saved AlphaGenome API key to remove.")
    return 0


def _cmd_fetch_genome(args) -> int:
    from .resources import ensure_hg38_fasta

    print(f"GRCh38 FASTA ready at: {ensure_hg38_fasta(args.cache_dir, verbose=True)}")
    return 0


def _cmd_fetch_spliceai(args) -> int:
    from .resources import ensure_spliceai_weights

    print(f"SpliceAI weights ready at: {ensure_spliceai_weights(args.cache_dir, verbose=True)}")
    return 0


def _cmd_fetch_ontology(args) -> int:
    from .ontology import save_ontology_snapshot

    print("Fetching AlphaGenome ontology metadata …")
    tsv, meta = save_ontology_snapshot(args.data_dir)
    print(f"Wrote {tsv}")
    print(f"Wrote {meta}")
    return 0


def _cmd_init(args) -> int:
    """One-time setup: download SpliceAI weights + save AlphaGenome API key."""
    from .resources import (
        ENV_VAR,
        ensure_spliceai_weights,
        get_alphagenome_api_key,
        save_alphagenome_api_key,
    )

    print("=== OligoMCP one-time setup ===")

    if args.skip_spliceai:
        print("\n[1/2] Skipping SpliceAI weight check (--skip-spliceai).")
    else:
        print("\n[1/2] Ensuring SpliceAI MANE-10000nt weights are available...")
        try:
            path = ensure_spliceai_weights(verbose=True)
            print(f"  ✓ Weights ready at: {path}")
        except Exception as e:
            print(f"  ! SpliceAI weight resolution failed: {e}")
            print("    You can retry later with `oligomcp fetch-spliceai-weights`.")

    print("\n[2/2] AlphaGenome API key")
    existing = get_alphagenome_api_key()
    if existing:
        print("  ✓ An AlphaGenome API key is already configured. Skipping.")
    elif args.skip_api_key:
        print("  Skipping (--skip-api-key). Save one later with `oligomcp set-api-key`.")
    elif args.yes or not sys.stdin.isatty():
        print(
            f"  Non-interactive mode — skipping prompt.\n"
            f"  Set ${ENV_VAR} or run `oligomcp set-api-key <KEY>` before running."
        )
    else:
        print(
            "  Get a free key at https://deepmind.google.com/science/alphagenome\n"
            "  (Leave blank to skip — you can set it later with `oligomcp set-api-key`.)"
        )
        try:
            key = getpass.getpass("  AlphaGenome API key (hidden, or Enter to skip): ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n  Skipped.")
            key = ""
        if key:
            try:
                path = save_alphagenome_api_key(key)
                print(f"  ✓ Saved to {path} (mode 0600).")
            except ValueError as e:
                print(f"  ! {e}")
        else:
            print("  Skipped. Save one later with `oligomcp set-api-key`.")

    print(
        "\nSetup complete. Genomic sequences are fetched on-the-fly from the UCSC API —\n"
        "no multi-GB genome download is required. Try:\n"
        "  oligomcp run --config config/SETD5_e1.json -v"
    )
    return 0


_HANDLERS = {
    "run": _cmd_run,
    "init": _cmd_init,
    "set-api-key": _cmd_set_api_key,
    "clear-api-key": _cmd_clear_api_key,
    "fetch-genome": _cmd_fetch_genome,
    "fetch-spliceai-weights": _cmd_fetch_spliceai,
    "fetch-ontology-terms": _cmd_fetch_ontology,
}


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    handler = _HANDLERS.get(args.cmd)
    if handler is None:
        parser.print_help()
        return 1
    return handler(args)


if __name__ == "__main__":
    sys.exit(main())
