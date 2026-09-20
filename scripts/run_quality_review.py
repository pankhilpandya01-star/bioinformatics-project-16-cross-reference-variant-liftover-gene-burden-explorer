"""Publish the Project 16 quality-review record from a passing JUnit report."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from variant_liftover_explorer.liftover import LiftoverError
from variant_liftover_explorer.quality_review import publish_quality_review
from variant_liftover_explorer.source_review import ReviewError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--junit-xml", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--project15-comparison-dir", type=Path, required=True)
    parser.add_argument("--project15-orthology-dir", type=Path, required=True)
    parser.add_argument("--project15-annotation-dir", type=Path, required=True)
    parser.add_argument("--connection-review-dir", type=Path, required=True)
    parser.add_argument("--authentic-dir", type=Path, required=True)
    parser.add_argument("--baseline-accession", default="GCF_000005845.2")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest = publish_quality_review(args)
    except (OSError, ValueError, json.JSONDecodeError, LiftoverError, ReviewError) as error:
        print(f"quality review failed: {error}", file=sys.stderr)
        return 1
    print(
        "quality review published: "
        f"{manifest['tests']['passed']} tests passed with no failures"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
