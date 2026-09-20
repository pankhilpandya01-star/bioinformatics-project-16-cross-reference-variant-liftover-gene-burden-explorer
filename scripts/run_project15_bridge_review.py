"""Validate the Project 15 tables required by the Project 16 integration layer."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from variant_liftover_explorer.liftover import LiftoverError
from variant_liftover_explorer.project15_bridge import run_project15_bridge_review


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-dir", type=Path, required=True)
    parser.add_argument("--annotation-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        manifest = run_project15_bridge_review(
            args.comparison_dir,
            args.annotation_dir,
            args.output_dir,
        )
    except (OSError, ValueError, LiftoverError) as error:
        print(f"Project 15 bridge review failed: {error}", file=sys.stderr)
        return 1
    print(
        "Project 15 bridge review published: "
        f"{manifest['counts']['comparable_gene_pairs']} comparable pairs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
