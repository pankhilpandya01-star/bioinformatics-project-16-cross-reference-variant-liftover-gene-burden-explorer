from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "authentic_analysis"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def vcf_keys(path: Path) -> list[str]:
    keys: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            keys.append(f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}")
    return keys


class AuthenticAnalysisResultTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))

    def test_projection_and_baseline_accounting(self) -> None:
        status_counts: Counter[tuple[str, str, str]] = Counter()
        keys: set[tuple[str, str, str]] = set()
        count = 0
        with gzip.open(
            RESULTS / "tables" / "variant_projections.csv.gz",
            "rt",
            newline="",
            encoding="utf-8",
        ) as handle:
            for row in csv.DictReader(handle):
                key = (
                    row["source_accession"],
                    row["target_accession"],
                    row["source_variant_key"],
                )
                self.assertNotIn(key, keys)
                keys.add(key)
                status_counts[(key[0], key[1], row["status"])] += 1
                count += 1
        self.assertEqual(count, 481_473)
        self.assertEqual(count, self.manifest["counts"]["directional_projection_attempts"])
        for direction, expected in self.manifest["direction_status_counts"].items():
            source, target = direction.split("_to_")
            for status, value in expected.items():
                self.assertEqual(status_counts[(source, target, status)], value)

        baseline_counts = Counter(
            row["source_variant_key"]
            for row in rows(RESULTS / "tables" / "baseline_variant_comparison.csv")
        )
        self.assertEqual(sum(baseline_counts.values()), 154_865)
        self.assertEqual(len(baseline_counts), 30_973)
        self.assertEqual(set(baseline_counts.values()), {5})

    def test_clusters_and_gene_pair_explanations_reconcile(self) -> None:
        clusters = rows(RESULTS / "tables" / "variant_equivalence_clusters.csv")
        self.assertEqual(len(clusters), 357_581)
        self.assertEqual(
            len({(row["accession"], row["variant_key"]) for row in clusters}),
            357_581,
        )
        cluster_sizes = Counter(row["cluster_id"] for row in clusters)
        self.assertEqual(len(cluster_sizes), 337_372)
        self.assertEqual(sum(size > 1 for size in cluster_sizes.values()), 14_409)

        explanations = rows(RESULTS / "tables" / "project15_burden_explanation.csv")
        self.assertEqual(len(explanations), 6_356)
        self.assertEqual(
            len(
                {
                    (row["baseline_gene_id"], row["alternative_accession"])
                    for row in explanations
                }
            ),
            6_356,
        )
        observed = Counter(row["burden_explanation"] for row in explanations)
        self.assertEqual(observed, Counter(self.manifest["burden_explanation_counts"]))
        same = [row for row in explanations if row["burden_explanation"] == "same_exact_variants"]
        self.assertEqual(len(same), 225)
        self.assertEqual(sum(int(row["exact_shared_variant_count"]) > 0 for row in same), 24)

    def test_effects_are_reported_only_for_accepted_exact_matches(self) -> None:
        equivalence = rows(RESULTS / "tables" / "ortholog_variant_equivalence.csv")
        self.assertEqual(len(equivalence), 44_097)
        target_effect_rows = [
            row for row in equivalence if row["effect_relation"] != "target_gene_effect_missing"
        ]
        self.assertEqual(len(target_effect_rows), 14_280)
        self.assertTrue(
            all(
                row["projection_status"] == "accepted_match"
                and row["target_in_ortholog_gene"] == "True"
                for row in target_effect_rows
            )
        )
        self.assertEqual(
            Counter(row["effect_relation"] for row in target_effect_rows),
            Counter(
                {
                    "same_consequence_and_impact": 14_201,
                    "same_impact_different_consequence": 2,
                    "different_impact": 77,
                }
            ),
        )

    def test_review_candidates_have_five_complete_contexts(self) -> None:
        review = rows(RESULTS / "tables" / "review_candidate_liftover.csv")
        self.assertEqual(len(review), 400)
        per_variant = Counter(row["baseline_variant_key"] for row in review)
        self.assertEqual(len(per_variant), 80)
        self.assertEqual(set(per_variant.values()), {5})
        self.assertEqual(sum(row["projection_status"] == "accepted_match" for row in review), 43)
        matched_variants = {
            row["baseline_variant_key"]
            for row in review
            if row["projection_status"] == "accepted_match"
        }
        self.assertEqual(len(matched_variants), 32)

    def test_all_directional_vcfs_are_unique_indexed_and_counted(self) -> None:
        observed_exact = 0
        for direction, expected in self.manifest["vcf_counts"].items():
            directory = RESULTS / "per_direction" / direction
            for stem, field in (
                ("projected.normalized.vcf.gz", "projected_vcf_records"),
                ("exact_equivalence.vcf.gz", "exact_vcf_records"),
            ):
                path = directory / stem
                self.assertTrue(path.is_file())
                self.assertTrue(Path(f"{path}.csi").is_file())
                keys = vcf_keys(path)
                self.assertEqual(len(keys), expected[field])
                self.assertEqual(len(keys), len(set(keys)))
                if field == "exact_vcf_records":
                    observed_exact += len(keys)
            self.assertEqual(expected["projected_duplicate_attempts"], 0)
            self.assertEqual(expected["exact_duplicate_attempts"], 0)
        self.assertEqual(observed_exact, self.manifest["counts"]["accepted_projection_edges"])

    def test_hashes_integrity_and_publication_hygiene(self) -> None:
        self.assertEqual(self.manifest["status"], "pass")
        self.assertTrue(self.manifest["source_inputs_unchanged"])
        for relative, expected in self.manifest["output_hashes"].items():
            self.assertEqual(sha256(RESULTS / relative), expected)

        text_files = [
            path
            for path in RESULTS.rglob("*")
            if path.is_file()
            and not path.name.endswith((".gz", ".csi", ".png"))
            and path.name != "run_manifest.json"
        ]
        combined = "\n".join(
            path.read_text(encoding="utf-8", errors="replace") for path in text_files
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
