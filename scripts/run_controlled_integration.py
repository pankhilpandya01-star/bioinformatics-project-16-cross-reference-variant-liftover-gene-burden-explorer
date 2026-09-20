"""Publish the controlled Project 15 integration demonstration."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from variant_liftover_explorer.controlled_integration import run_controlled_integration
from variant_liftover_explorer.liftover import LiftoverError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        manifest = run_controlled_integration(args.output_dir)
    except (OSError, ValueError, LiftoverError) as error:
        print(f"controlled integration failed: {error}", file=sys.stderr)
        return 1
    print(
        "controlled integration published: "
        f"{manifest['counts']['comparable_gene_pairs']} gene-pair explanations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
