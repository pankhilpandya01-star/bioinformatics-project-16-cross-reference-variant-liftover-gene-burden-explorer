"""Audit the local Project 16 publication candidate without creating a repository."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from variant_liftover_explorer.liftover import LiftoverError
from variant_liftover_explorer.publication_audit import publish_audit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/publication_audit")
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = publish_audit(args.project_root, args.output_dir)
    except (OSError, ValueError, json.JSONDecodeError, LiftoverError) as error:
        print(f"publication audit failed: {error}", file=sys.stderr)
        return 1
    print(
        "publication candidate passed: "
        f"{summary['candidate_files']} files, {summary['candidate_bytes']} bytes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
