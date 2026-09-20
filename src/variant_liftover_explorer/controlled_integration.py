"""Controlled Project 15 integration with one example per burden category."""

from __future__ import annotations

import csv
import json
import os
import shutil
import uuid
from collections import Counter
from pathlib import Path

from .integration import (
    BURDEN_CATEGORIES,
    GenePair,
    GeneVariant,
    ProjectionEvidence,
    connect_review_candidates,
    integrate_gene_burdens,
    write_rows,
)
from .liftover import LiftoverError


BASELINE = "BASELINE.1"
ALTERNATIVE = "ALTERNATIVE.1"


def controlled_records() -> tuple[
    list[GenePair], list[GeneVariant], list[ProjectionEvidence], list[dict[str, str]]
]:
    pairs = [
        GenePair(BASELINE, "G_same", "same", ALTERNATIVE, "A_same", "same"),
        GenePair(BASELINE, "G_partial", "partial", ALTERNATIVE, "A_partial", "partial"),
        GenePair(BASELINE, "G_reference", "reference", ALTERNATIVE, "A_reference", "reference"),
        GenePair(BASELINE, "G_filter", "filter", ALTERNATIVE, "A_filter", "filter"),
        GenePair(BASELINE, "G_different", "different", ALTERNATIVE, "A_different", "different"),
        GenePair(BASELINE, "G_unresolved", "unresolved", ALTERNATIVE, "A_unresolved", "unresolved"),
    ]
    variants = [
        GeneVariant(BASELINE, "B:10:A:G", "G_same", "same", "missense_variant", "MODERATE"),
        GeneVariant(BASELINE, "B:20:C:T", "G_same", "same", "synonymous_variant", "LOW"),
        GeneVariant(BASELINE, "B:25:A:T", "G_same", "same", "missense_variant", "MODERATE"),
        GeneVariant(ALTERNATIVE, "A:11:A:G", "A_same", "same", "missense_variant", "MODERATE"),
        GeneVariant(ALTERNATIVE, "A:21:C:T", "A_same", "same", "stop_retained_variant", "LOW"),
        GeneVariant(ALTERNATIVE, "A:26:A:T", "A_same", "same", "stop_gained", "HIGH"),
        GeneVariant(BASELINE, "B:30:G:A", "G_partial", "partial", "missense_variant", "MODERATE"),
        GeneVariant(BASELINE, "B:40:T:C", "G_partial", "partial", "synonymous_variant", "LOW"),
        GeneVariant(ALTERNATIVE, "A:31:G:A", "A_partial", "partial", "missense_variant", "MODERATE"),
        GeneVariant(ALTERNATIVE, "A:45:C:G", "A_partial", "partial", "frameshift_variant", "HIGH"),
        GeneVariant(BASELINE, "B:50:A:C", "G_reference", "reference", "missense_variant", "MODERATE"),
        GeneVariant(BASELINE, "B:60:G:T", "G_filter", "filter", "missense_variant", "MODERATE"),
        GeneVariant(BASELINE, "B:70:C:A", "G_different", "different", "synonymous_variant", "LOW"),
        GeneVariant(ALTERNATIVE, "A:75:T:G", "A_different", "different", "missense_variant", "MODERATE"),
        GeneVariant(BASELINE, "B:80:T:A", "G_unresolved", "unresolved", "frameshift_variant", "HIGH"),
        GeneVariant(ALTERNATIVE, "A:85:G:C", "A_unresolved", "unresolved", "frameshift_variant", "HIGH"),
    ]
    projections = [
        ProjectionEvidence(ALTERNATIVE, "B:10:A:G", "accepted_match", "A:11:A:G", "exact key"),
        ProjectionEvidence(ALTERNATIVE, "B:20:C:T", "accepted_match", "A:21:C:T", "exact key"),
        ProjectionEvidence(ALTERNATIVE, "B:25:A:T", "accepted_match", "A:26:A:T", "exact key"),
        ProjectionEvidence(ALTERNATIVE, "B:30:G:A", "accepted_match", "A:31:G:A", "exact key"),
        ProjectionEvidence(ALTERNATIVE, "B:40:T:C", "callable_no_candidate", "A:41:T:C", "absent"),
        ProjectionEvidence(ALTERNATIVE, "B:50:A:C", "target_reference_match", "", "reference allele"),
        ProjectionEvidence(ALTERNATIVE, "B:60:G:T", "rejected_match", "A:61:G:T", "filtered"),
        ProjectionEvidence(ALTERNATIVE, "B:70:C:A", "callable_no_candidate", "A:71:C:A", "absent"),
        ProjectionEvidence(ALTERNATIVE, "B:80:T:A", "ambiguous_mapping", "", "repeat"),
    ]
    review_contexts = [
        {
            "baseline_variant_key": "B:10:A:G",
            "baseline_gene_id": "G_same",
            "baseline_gene_name": "same",
            "baseline_primary_consequence": "missense_variant",
            "baseline_primary_impact": "MODERATE",
            "alternative_accession": ALTERNATIVE,
            "alternative_gene_ids": "A_same",
            "comparison_scope": "controlled_gene_context",
        },
        {
            "baseline_variant_key": "B:50:A:C",
            "baseline_gene_id": "G_reference",
            "baseline_gene_name": "reference",
            "baseline_primary_consequence": "missense_variant",
            "baseline_primary_impact": "MODERATE",
            "alternative_accession": ALTERNATIVE,
            "alternative_gene_ids": "A_reference",
            "comparison_scope": "controlled_gene_context",
        },
    ]
    return pairs, variants, projections, review_contexts


def run_controlled_integration(output_dir: Path) -> dict[str, object]:
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise LiftoverError("output directory must not already exist")
    stage = output_dir.parent / f".{output_dir.name}.staging-{uuid.uuid4().hex}"
    stage.mkdir(parents=True)
    (stage / "tables").mkdir()
    (stage / "summaries").mkdir()
    try:
        pairs, variants, projections, review_contexts = controlled_records()
        equivalence, explanations = integrate_gene_burdens(pairs, variants, projections)
        review = connect_review_candidates(
            review_contexts,
            projections,
            equivalence,
            explanations,
        )
        counts = Counter(row["burden_explanation"] for row in explanations)
        if set(counts) != set(BURDEN_CATEGORIES) or any(counts[item] != 1 for item in counts):
            raise LiftoverError(f"controlled burden-category accounting failed: {dict(counts)}")
        write_rows(stage / "tables" / "ortholog_variant_equivalence.csv", equivalence)
        write_rows(stage / "tables" / "project15_burden_explanation.csv", explanations)
        write_rows(stage / "tables" / "review_candidate_liftover.csv", review)
        with (stage / "summaries" / "burden_explanation_summary.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=["burden_explanation", "count"], lineterminator="\n")
            writer.writeheader()
            for category in BURDEN_CATEGORIES:
                writer.writerow({"burden_explanation": category, "count": counts[category]})
        effect_counts = Counter(row["effect_relation"] for row in equivalence)
        with (stage / "summaries" / "effect_relation_summary.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=["effect_relation", "count"], lineterminator="\n")
            writer.writeheader()
            for relation, count in sorted(effect_counts.items()):
                writer.writerow({"effect_relation": relation, "count": count})
        manifest = {
            "workflow": "variant-liftover-explorer",
            "workflow_version": "0.1.0",
            "milestone": "controlled_project15_integration",
            "status": "pass",
            "counts": {
                "comparable_gene_pairs": len(explanations),
                "variant_gene_projection_rows": len(equivalence),
                "review_candidate_contexts": len(review),
                "burden_categories": dict(counts),
                "effect_relations": dict(effect_counts),
            },
            "category_precedence": [
                "unresolved_liftover",
                "same_exact_variants",
                "partially_shared_variants",
                "reference_allele_difference",
                "callability_or_filter_difference",
                "different_variants_in_same_gene",
            ],
            "atomic_publication": True,
            "interpretation_boundary": (
                "Gene-level integration explains controlled burden patterns without asserting phenotype."
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
