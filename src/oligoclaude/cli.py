"""CLI entrypoint for OligoClaude."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .workflow import run_workflow


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="oligoclaude",
        description="Predict ASO efficacy via AlphaGenome and SpliceAI, and "
        "compare to experimental data.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser(
        "run", help="Run the full ASO prediction workflow for one config."
    )
    p_run.add_argument("--config", type=Path, required=True, help="Path to JSON config.")
    p_run.add_argument(
        "--skip-alphagenome", action="store_true", help="Skip AlphaGenome scoring."
    )
    p_run.add_argument(
        "--skip-spliceai", action="store_true", help="Skip SpliceAI scoring."
    )
    p_run.add_argument(
        "--samples-max",
        type=int,
        default=20,
        help="Top/bottom ASOs to emit in compact BED.",
    )
    p_run.add_argument("--verbose", "-v", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "run":
        result = run_workflow(
            config_path=args.config,
            skip_alphagenome=args.skip_alphagenome,
            skip_spliceai=args.skip_spliceai,
            samples_max=args.samples_max,
            verbose=args.verbose,
        )
        print()
        print(f"Candidates scored: {result.n_candidates}")
        print(f"Scores CSV: {result.scores_csv}")
        if result.correlation_plot:
            print(f"Correlation plot: {result.correlation_plot}")
        print()
        print(result.ucsc_instructions)
        return 0

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
