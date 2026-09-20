"""Run Project 16 authentic cross-reference projection and explanation."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from variant_liftover_explorer.authentic_analysis import authentic_analysis
from variant_liftover_explorer.liftover import LiftoverError
from variant_liftover_explorer.source_review import ReviewError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--project15-comparison-dir", type=Path, required=True)
    parser.add_argument("--project15-orthology-dir", type=Path, required=True)
    parser.add_argument("--project15-annotation-dir", type=Path, required=True)
    parser.add_argument("--connection-review-dir", type=Path, required=True)
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--baseline-accession", default="GCF_000005845.2")
    parser.add_argument("--sample-id", default="SRR13921545")
    parser.add_argument("--alignment-preset", default="asm20")
    parser.add_argument("--minimum-alignment-length", type=int, default=10_000)
    parser.add_argument("--minimum-alignment-identity", type=float, default=0.85)
    parser.add_argument("--minimum-alignment-mapq", type=int, default=20)
    parser.add_argument("--ambiguity-score-fraction", type=float, default=0.95)
    parser.add_argument("--minimum-comparable-coverage", type=float, default=0.70)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest = authentic_analysis(args)
    except (
        OSError,
        ValueError,
        json.JSONDecodeError,
        subprocess.SubprocessError,
        LiftoverError,
        ReviewError,
    ) as error:
        print(f"authentic analysis failed: {error}", file=sys.stderr)
        return 1
    print(
        "authentic analysis published: "
        f"{manifest['counts']['directional_projection_attempts']} projection attempts"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
