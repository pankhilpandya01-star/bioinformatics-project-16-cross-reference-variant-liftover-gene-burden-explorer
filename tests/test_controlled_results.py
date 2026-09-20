from __future__ import annotations

import csv
import gzip
import json
import unittest
from collections import Counter
from pathlib import Path

from variant_liftover_explorer.controlled_demo import EXPECTED_STATUSES


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "controlled_core"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class ControlledResultTests(unittest.TestCase):
    def test_every_controlled_variant_has_one_expected_status(self) -> None:
        projections = rows(RESULTS / "tables" / "controlled_projections.csv")
        self.assertEqual(len(projections), 12)
        self.assertEqual(len({row["source_variant_key"] for row in projections}), 12)
        self.assertEqual(Counter(row["status"] for row in projections), Counter(EXPECTED_STATUSES))

    def test_reverse_strand_and_indels_are_exact_matches(self) -> None:
        projections = {
            row["source_variant_id"]: row
            for row in rows(RESULTS / "tables" / "controlled_projections.csv")
        }
        self.assertEqual(projections["reverse_snp"]["target_variant_key"], "tgt_reverse:17:T:C")
        self.assertEqual(projections["reverse_snp"]["strand"], "-")
        self.assertEqual(projections["accepted_insertion"]["target_variant_key"], "tgt_plus:10:C:CAA")
        self.assertEqual(projections["accepted_deletion"]["target_variant_key"], "tgt_plus:13:AC:A")
        self.assertEqual(projections["target_reference"]["target_variant_key"], "")
        for identifier in ("reverse_snp", "accepted_insertion", "accepted_deletion"):
            self.assertEqual(projections[identifier]["status"], "accepted_match")
            self.assertEqual(projections[identifier]["reciprocal"], "true")

    def test_manifest_and_compressed_vcf_accounting_agree(self) -> None:
        manifest = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))
        with gzip.open(RESULTS / "reports" / "normalized_projected.vcf.gz", "rt") as handle:
            vcf_records = sum(not line.startswith("#") for line in handle)
        self.assertEqual(manifest["status"], "pass")
        self.assertTrue(manifest["atomic_publication"])
        self.assertEqual(manifest["counts"]["source_variants"], 12)
        self.assertEqual(manifest["counts"]["projected_vcf_records"], vcf_records)
        self.assertEqual(vcf_records, 7)
        self.assertTrue((RESULTS / "reports" / "normalized_projected.vcf.gz.csi").is_file())

    def test_public_outputs_are_portable(self) -> None:
        combined = "\n".join(
            path.read_text(encoding="utf-8", errors="replace")
            for path in RESULTS.rglob("*")
            if path.is_file() and not path.name.endswith((".gz", ".csi"))
        )
        self.assertNotIn("C:" + "\\Users\\", combined)
        self.assertNotIn("/mnt/" + "c/Users/", combined)
        self.assertNotIn("One" + "Drive", combined)


if __name__ == "__main__":
    unittest.main()
