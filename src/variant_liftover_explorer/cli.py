"""Command-line interface for exact cross-reference variant liftover."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .workflow import WorkflowConfig, WorkflowError, run_workflow


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="variant-liftover-explorer",
        description=(
            "Align six bacterial references, project normalized variants reciprocally, "
            "and explain Project 15 gene-burden differences."
        ),
    )
    parser.add_argument("--reference-manifest", required=True, type=Path)
    parser.add_argument("--project15-comparison-dir", required=True, type=Path)
    parser.add_argument("--project15-orthology-dir", required=True, type=Path)
    parser.add_argument("--project15-annotation-dir", required=True, type=Path)
    parser.add_argument("--baseline-accession", default="GCF_000005845.2")
    parser.add_argument("--sample-id", default="SRR13921545")
    parser.add_argument("--alignment-preset", default="asm20")
    parser.add_argument("--minimum-alignment-length", type=int, default=10_000)
    parser.add_argument("--minimum-alignment-identity", type=float, default=0.85)
    parser.add_argument("--minimum-alignment-mapq", type=int, default=20)
    parser.add_argument("--ambiguity-score-fraction", type=float, default=0.95)
    parser.add_argument("--minimum-comparable-coverage", type=float, default=0.70)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        output = run_workflow(
            WorkflowConfig(
                reference_manifest=args.reference_manifest,
                project15_comparison_dir=args.project15_comparison_dir,
                project15_orthology_dir=args.project15_orthology_dir,
                project15_annotation_dir=args.project15_annotation_dir,
                output_dir=args.output_dir,
                baseline_accession=args.baseline_accession,
                sample_id=args.sample_id,
                alignment_preset=args.alignment_preset,
                minimum_alignment_length=args.minimum_alignment_length,
                minimum_alignment_identity=args.minimum_alignment_identity,
                minimum_alignment_mapq=args.minimum_alignment_mapq,
                ambiguity_score_fraction=args.ambiguity_score_fraction,
                minimum_comparable_coverage=args.minimum_comparable_coverage,
                threads=args.threads,
            )
        )
    except WorkflowError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(f"Completed cross-reference variant liftover: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
