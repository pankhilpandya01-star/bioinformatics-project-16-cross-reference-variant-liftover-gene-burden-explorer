from __future__ import annotations

import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "connection_review"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class ConnectionReviewResultTests(unittest.TestCase):
    def test_project15_totals_are_exactly_reproduced(self) -> None:
        observed = {
            row["metric"]: (int(row["expected"]), int(row["observed"]), row["status"])
            for row in rows(RESULTS / "tables" / "project15_reproduction.csv")
        }
        self.assertEqual(observed["accepted_variants"], (357_581, 357_581, "pass"))
        self.assertEqual(observed["annotation_effects"], (357_959, 357_959, "pass"))
        self.assertEqual(observed["comparable_gene_pairs"], (6_356, 6_356, "pass"))
        self.assertEqual(observed["review_candidate_reference_contexts"], (400, 400, "pass"))

    def test_alignment_accounting_and_coverage_labels_agree(self) -> None:
        manifest = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))
        blocks = rows(RESULTS / "tables" / "alignment_blocks.csv")
        coverage = rows(RESULTS / "tables" / "reciprocal_coverage.csv")
        self.assertEqual(manifest["status"], "pass")
        self.assertEqual(manifest["alignment_counts"]["directions"], 10)
        self.assertEqual(len(blocks), manifest["alignment_counts"]["paf_records"])
        self.assertEqual(
            sum(row["qualifying"] == "true" for row in blocks),
            manifest["alignment_counts"]["qualifying_blocks"],
        )
        self.assertEqual(
            sum(row["reciprocal"] == "true" for row in blocks),
            manifest["alignment_counts"]["reciprocal_blocks"],
        )
        self.assertEqual(len(coverage), 5)
        self.assertTrue(all(row["substantive_comparison"] == "false" for row in coverage))

    def test_commands_and_public_results_are_portable(self) -> None:
        commands = (RESULTS / "portable_commands.txt").read_text(encoding="utf-8")
        self.assertEqual(len(commands.splitlines()), 10)
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
