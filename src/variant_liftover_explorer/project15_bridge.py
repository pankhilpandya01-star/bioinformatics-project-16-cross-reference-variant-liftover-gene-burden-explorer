"""Read-only validation of the Project 15 tables consumed by Project 16."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import uuid
from pathlib import Path

from .integration import (
    load_project15_gene_pairs,
    load_project15_gene_variants,
    load_project15_review_contexts,
    validate_project15_burdens,
)
from .liftover import LiftoverError


EXPECTED = {
    "comparison_rows": 23_255,
    "comparable_gene_pairs": 6_356,
    "annotation_effects": 357_959,
    "variant_gene_assignments": 327_542,
    "unlinked_effects": 30_382,
    "gene_burden_rows": 28_859,
    "review_contexts": 400,
    "review_candidates": 80,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_project15_bridge_review(
    comparison_dir: Path,
    annotation_dir: Path,
    output_dir: Path,
) -> dict[str, object]:
    comparison_dir = comparison_dir.resolve()
    annotation_dir = annotation_dir.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise LiftoverError("output directory must not already exist")
    paths = {
        "baseline_ortholog_comparison": comparison_dir / "tables" / "baseline_ortholog_comparison.csv",
        "gene_consequence_burden": comparison_dir / "tables" / "gene_consequence_burden.csv",
        "review_candidate_orthology": comparison_dir / "tables" / "review_candidate_orthology.csv",
        "annotation_effects": annotation_dir / "tables" / "annotation_effects.csv",
        "gene_catalog": annotation_dir / "tables" / "gene_catalog.csv",
    }
    for name, path in paths.items():
        if not path.is_file():
            raise LiftoverError(f"missing Project 15 bridge input {name}: {path}")
    hashes_before = {name: _sha256(path) for name, path in paths.items()}
    pairs, comparison_count = load_project15_gene_pairs(paths["baseline_ortholog_comparison"])
    variants, effect_count, unlinked = load_project15_gene_variants(
        paths["annotation_effects"], paths["gene_catalog"]
    )
    burden_count = validate_project15_burdens(paths["gene_consequence_burden"], variants)
    reviews = load_project15_review_contexts(paths["review_candidate_orthology"])
    observed = {
        "comparison_rows": comparison_count,
        "comparable_gene_pairs": len(pairs),
        "annotation_effects": effect_count,
        "variant_gene_assignments": len(variants),
        "unlinked_effects": unlinked,
        "gene_burden_rows": burden_count,
        "review_contexts": len(reviews),
        "review_candidates": len({row["baseline_variant_key"] for row in reviews}),
    }
    if observed != EXPECTED:
        raise LiftoverError(f"Project 15 bridge counts disagree: {observed}")

    stage = output_dir.parent / f".{output_dir.name}.staging-{uuid.uuid4().hex}"
    stage.mkdir(parents=True)
    (stage / "tables").mkdir()
    try:
        with (stage / "tables" / "project15_bridge_counts.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["metric", "expected", "observed", "status"],
                lineterminator="\n",
            )
            writer.writeheader()
            for metric, expected in EXPECTED.items():
                writer.writerow(
                    {
                        "metric": metric,
                        "expected": expected,
                        "observed": observed[metric],
                        "status": "pass",
                    }
                )
        hashes_after = {name: _sha256(path) for name, path in paths.items()}
        if hashes_after != hashes_before:
            raise LiftoverError("a Project 15 bridge input changed during review")
        manifest = {
            "workflow": "variant-liftover-explorer",
            "workflow_version": "0.1.0",
            "milestone": "project15_bridge_review",
            "status": "pass",
            "counts": observed,
            "input_hashes": hashes_before,
            "source_inputs_unchanged": True,
            "atomic_publication": True,
            "interpretation_boundary": (
                "This review validates inherited schemas and accounting; it does not project authentic variants."
            ),
        }
        (stage / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(stage, output_dir)
        return manifest
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
