"""Publish the Project 16 controlled liftover demonstration."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from variant_liftover_explorer.controlled_demo import run_controlled_demo
from variant_liftover_explorer.liftover import LiftoverError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        manifest = run_controlled_demo(args.output_dir)
    except (OSError, ValueError, subprocess.SubprocessError, LiftoverError) as error:
        print(f"controlled demonstration failed: {error}", file=sys.stderr)
        return 1
    print(
        f"controlled demonstration published: {manifest['counts']['source_variants']} variants"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
