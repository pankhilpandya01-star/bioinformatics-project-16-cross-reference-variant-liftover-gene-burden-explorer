#!/usr/bin/env python3
"""Run Project 16 Milestone 1 without changing inherited source artifacts."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from variant_liftover_explorer.source_review import main


if __name__ == "__main__":
    raise SystemExit(main())

