from __future__ import annotations

import csv
import json
import unittest
from collections import Counter
from pathlib import Path

from variant_liftover_explorer.integration import BURDEN_CATEGORIES


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "controlled_integration"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class ControlledIntegrationResultTests(unittest.TestCase):
    def test_every_burden_category_is_demonstrated_once(self) -> None:
        explanations = rows(RESULTS / "tables" / "project15_burden_explanation.csv")
        counts = Counter(row["burden_explanation"] for row in explanations)
        self.assertEqual(len(explanations), 6)
        self.assertEqual(counts, Counter({category: 1 for category in BURDEN_CATEGORIES}))

    def test_effect_and_review_context_links_are_accounted(self) -> None:
        equivalence = rows(RESULTS / "tables" / "ortholog_variant_equivalence.csv")
        review = rows(RESULTS / "tables" / "review_candidate_liftover.csv")
        manifest = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(len(equivalence), 9)
        self.assertEqual(len(review), 2)
        self.assertEqual(manifest["counts"]["variant_gene_projection_rows"], len(equivalence))
        self.assertEqual(manifest["counts"]["review_candidate_contexts"], len(review))
        self.assertEqual(
            {row["effect_relation"] for row in equivalence},
            {
                "same_consequence_and_impact",
                "same_impact_different_consequence",
                "different_impact",
                "target_gene_effect_missing",
            },
        )
        self.assertEqual(sum(row["exact_accepted_match"] == "True" for row in review), 1)
        self.assertEqual(review[0]["effect_relation"], "same_consequence_and_impact")
        self.assertEqual(review[0]["burden_explanation"], "same_exact_variants")
        self.assertEqual(review[1]["burden_explanation"], "reference_allele_difference")

    def test_outputs_are_portable_and_atomically_published(self) -> None:
        manifest = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))
        combined = "\n".join(
            path.read_text(encoding="utf-8", errors="replace")
            for path in RESULTS.rglob("*")
            if path.is_file()
        )
        self.assertTrue(manifest["atomic_publication"])
        self.assertNotIn("C:" + "\\Users\\", combined)
        self.assertNotIn("/mnt/" + "c/Users/", combined)
        self.assertNotIn("One" + "Drive", combined)


if __name__ == "__main__":
    unittest.main()
