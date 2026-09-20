from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from variant_liftover_explorer import controlled_integration
from variant_liftover_explorer.controlled_integration import (
    ALTERNATIVE,
    BASELINE,
    controlled_records,
    run_controlled_integration,
)
from variant_liftover_explorer.integration import (
    GenePair,
    GeneVariant,
    ProjectionEvidence,
    connect_review_candidates,
    integrate_gene_burdens,
    validate_project15_burdens,
)
from variant_liftover_explorer.liftover import LiftoverError


class IntegrationUnitTests(unittest.TestCase):
    def test_controlled_records_assign_every_burden_category(self) -> None:
        pairs, variants, projections, _reviews = controlled_records()
        equivalence, explanations = integrate_gene_burdens(pairs, variants, projections)
        self.assertEqual(len(explanations), 6)
        self.assertEqual(
            {row["burden_explanation"] for row in explanations},
            {
                "same_exact_variants",
                "partially_shared_variants",
                "reference_allele_difference",
                "callability_or_filter_difference",
                "different_variants_in_same_gene",
                "unresolved_liftover",
            },
        )
        self.assertEqual(
            {row["effect_relation"] for row in equivalence},
            {
                "same_consequence_and_impact",
                "same_impact_different_consequence",
                "different_impact",
                "target_gene_effect_missing",
            },
        )

    def test_unresolved_projection_takes_conservative_precedence(self) -> None:
        pairs = [GenePair(BASELINE, "G", "gene", ALTERNATIVE, "A", "gene")]
        variants = [
            GeneVariant(BASELINE, "B:1:A:G", "G", "gene", "missense_variant", "MODERATE"),
            GeneVariant(BASELINE, "B:2:C:T", "G", "gene", "synonymous_variant", "LOW"),
            GeneVariant(ALTERNATIVE, "A:1:A:G", "A", "gene", "missense_variant", "MODERATE"),
        ]
        projections = [
            ProjectionEvidence(ALTERNATIVE, "B:1:A:G", "accepted_match", "A:1:A:G"),
            ProjectionEvidence(ALTERNATIVE, "B:2:C:T", "nonreciprocal_mapping"),
        ]
        _equivalence, explanations = integrate_gene_burdens(pairs, variants, projections)
        self.assertEqual(explanations[0]["burden_explanation"], "unresolved_liftover")

    def test_duplicate_and_missing_projections_fail(self) -> None:
        pairs, variants, projections, _reviews = controlled_records()
        with self.assertRaisesRegex(LiftoverError, "duplicate directional projection"):
            integrate_gene_burdens(pairs, variants, projections + [projections[0]])
        with self.assertRaisesRegex(LiftoverError, "missing projection"):
            integrate_gene_burdens(pairs, variants, projections[1:])

    def test_review_candidate_join_is_exact_and_complete(self) -> None:
        pairs, variants, projections, reviews = controlled_records()
        equivalence, explanations = integrate_gene_burdens(pairs, variants, projections)
        output = connect_review_candidates(reviews, projections, equivalence, explanations)
        self.assertEqual(len(output), len(reviews))
        self.assertEqual([row["projection_status"] for row in output], [
            "accepted_match",
            "target_reference_match",
        ])
        self.assertEqual(output[0]["effect_relation"], "same_consequence_and_impact")
        self.assertEqual(output[0]["burden_explanation"], "same_exact_variants")
        self.assertEqual(output[1]["burden_explanation"], "reference_allele_difference")
        with self.assertRaisesRegex(LiftoverError, "lacks projection"):
            connect_review_candidates(
                [{"baseline_variant_key": "B:999:A:C", "alternative_accession": ALTERNATIVE}],
                projections,
            )

    def test_noncomparable_pairs_are_not_explained(self) -> None:
        pair = GenePair(
            BASELINE,
            "G",
            "gene",
            ALTERNATIVE,
            "A",
            "gene",
            comparison_status="low_callability",
            comparable=False,
        )
        equivalence, explanations = integrate_gene_burdens([pair], [], [])
        self.assertEqual(equivalence, [])
        self.assertEqual(explanations, [])

    def test_atomic_integration_cleans_stage_after_failure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "integration"
            with mock.patch.object(
                controlled_integration,
                "integrate_gene_burdens",
                side_effect=LiftoverError("injected integration failure"),
            ):
                with self.assertRaisesRegex(LiftoverError, "injected integration failure"):
                    run_controlled_integration(output)
            self.assertFalse(output.exists())
            self.assertEqual(list(Path(directory).glob(".integration.staging-*")), [])

    def test_burden_validation_rejects_count_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            burden = Path(directory) / "burden.csv"
            burden.write_text(
                "accession,gene_id,distinct_variant_count\n"
                "BASELINE.1,G,2\n",
                encoding="utf-8",
            )
            variants = [
                GeneVariant(
                    BASELINE,
                    "B:1:A:G",
                    "G",
                    "gene",
                    "missense_variant",
                    "MODERATE",
                )
            ]
            with self.assertRaisesRegex(LiftoverError, "burden disagrees"):
                validate_project15_burdens(burden, variants)


if __name__ == "__main__":
    unittest.main()
