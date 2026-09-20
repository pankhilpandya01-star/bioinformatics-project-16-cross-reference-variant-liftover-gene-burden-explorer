from __future__ import annotations

import csv
import hashlib
import json
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "quality_review"
AUTHENTIC = ROOT / "results" / "authentic_analysis"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class QualityReviewResultTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))

    def test_quality_review_passed_without_skips_or_partial_outputs(self) -> None:
        tests = self.manifest["tests"]
        self.assertEqual(self.manifest["status"], "pass")
        self.assertEqual(tests["tests"], tests["passed"])
        self.assertEqual(tests["failures"], 0)
        self.assertEqual(tests["errors"], 0)
        self.assertEqual(tests["skipped"], 0)
        self.assertTrue(self.manifest["source_inputs_unchanged"])
        self.assertEqual(self.manifest["partial_outputs_after_injected_failures"], 0)
        self.assertFalse(self.manifest["repository_created"])

    def test_matrix_and_authentic_checksums_agree(self) -> None:
        with (RESULTS / "test_matrix.csv").open(newline="", encoding="utf-8") as handle:
            matrix = list(csv.DictReader(handle))
        self.assertEqual(len(matrix), 10)
        self.assertTrue(all(row["status"] == "pass" for row in matrix))
        checksums = self.manifest["artifact_checksums"]
        self.assertEqual(checksums["test_matrix"], sha256(RESULTS / "test_matrix.csv"))
        self.assertEqual(
            checksums["authentic_manifest"], sha256(AUTHENTIC / "run_manifest.json")
        )
        self.assertEqual(
            checksums["authentic_dashboard"],
            sha256(AUTHENTIC / "dashboard" / "variant_liftover_dashboard.png"),
        )

    def test_quality_outputs_are_portable_and_publication_safe(self) -> None:
        combined = "\n".join(
            path.read_text(encoding="utf-8", errors="replace")
            for path in RESULTS.rglob("*")
            if path.is_file()
        )
        self.assertNotIn("C:" + "\\Users\\", combined)
        self.assertNotIn("/mnt/" + "c/Users/", combined)
        self.assertNotIn("One" + "Drive", combined)
        blocked_authorship_terms = (
            "chat" + "gpt",
            "open" + "ai",
            "artificial " + "intelligence",
            "ai" + "-generated",
        )
        self.assertTrue(
            all(
                re.search(rf"\b{re.escape(term)}\b", combined, re.I) is None
                for term in blocked_authorship_terms
            )
        )


if __name__ == "__main__":
    unittest.main()
