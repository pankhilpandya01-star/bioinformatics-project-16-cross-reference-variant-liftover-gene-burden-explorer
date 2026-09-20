from __future__ import annotations

import csv
import json
import unittest
from pathlib import Path

from variant_liftover_explorer.project15_bridge import EXPECTED


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "project15_integration_review"


class Project15BridgeResultTests(unittest.TestCase):
    def test_inherited_counts_are_exactly_reproduced(self) -> None:
        manifest = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))
        with (RESULTS / "tables" / "project15_bridge_counts.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            table = {row["metric"]: row for row in csv.DictReader(handle)}
        self.assertEqual(manifest["status"], "pass")
        self.assertEqual(manifest["counts"], EXPECTED)
        for metric, expected in EXPECTED.items():
            self.assertEqual(int(table[metric]["expected"]), expected)
            self.assertEqual(int(table[metric]["observed"]), expected)
            self.assertEqual(table[metric]["status"], "pass")

    def test_hashes_match_the_connection_review(self) -> None:
        bridge = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))
        connection = json.loads(
            (ROOT / "results" / "connection_review" / "run_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        expected = connection["project15_hashes"]
        self.assertEqual(
            bridge["input_hashes"]["baseline_ortholog_comparison"],
            expected["baseline_ortholog_comparison"],
        )
        self.assertEqual(
            bridge["input_hashes"]["gene_consequence_burden"],
            expected["gene_consequence_burden"],
        )
        self.assertEqual(
            bridge["input_hashes"]["review_candidate_orthology"],
            expected["review_candidate_orthology"],
        )
        self.assertEqual(
            bridge["input_hashes"]["annotation_effects"],
            expected["annotation_effects"],
        )
        self.assertTrue(bridge["source_inputs_unchanged"])
        self.assertTrue(bridge["atomic_publication"])

    def test_public_bridge_outputs_contain_no_machine_paths(self) -> None:
        combined = "\n".join(
            path.read_text(encoding="utf-8", errors="replace")
            for path in RESULTS.rglob("*")
            if path.is_file()
        )
        self.assertNotIn("C:" + "\\Users\\", combined)
        self.assertNotIn("/mnt/" + "c/Users/", combined)
        self.assertNotIn("One" + "Drive", combined)


if __name__ == "__main__":
    unittest.main()
