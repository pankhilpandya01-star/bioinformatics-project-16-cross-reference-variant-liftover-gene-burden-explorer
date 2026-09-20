from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from variant_liftover_explorer.liftover import LiftoverError
from variant_liftover_explorer.publication_audit import publish_audit


class PublicationAuditTests(unittest.TestCase):
    def test_safe_candidate_is_inventoried_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "docs").mkdir()
            (root / "docs" / "methods.md").write_text("# Methods\n", encoding="utf-8")
            (root / "README.md").write_text(
                "# Controlled\n\n[Methods](docs/methods.md)\n", encoding="utf-8"
            )
            (root / "local").mkdir()
            (root / "local" / "private.txt").write_text(
                "ignored local data", encoding="utf-8"
            )
            output = root / "results" / "publication_audit"
            output.parent.mkdir()
            summary = publish_audit(root, output)
            self.assertEqual(summary["status"], "pass")
            self.assertEqual(summary["candidate_files"], 2)
            self.assertTrue(summary["local_storage_excluded"])
            self.assertFalse(summary["repository_created"])
            recorded = json.loads((output / "audit_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(recorded, summary)

    def test_private_path_and_broken_link_fail_without_partial_output(self) -> None:
        cases = {
            "private_path": "native=" + "C:" + "\\Users\\someone\\data\n",
            "broken_link": "[Missing](docs/missing.md)\n",
        }
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            for name, content in cases.items():
                with self.subTest(name=name):
                    root = base / name
                    root.mkdir()
                    (root / "README.md").write_text(content, encoding="utf-8")
                    output = root / "results" / "publication_audit"
                    output.parent.mkdir()
                    with self.assertRaisesRegex(LiftoverError, "publication audit failed"):
                        publish_audit(root, output)
                    self.assertFalse(output.exists())
                    self.assertEqual(list(output.parent.glob(".publication_audit.staging-*")), [])

    def test_output_must_be_inside_project_and_new(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "project"
            root.mkdir()
            (root / "README.md").write_text("# Project\n", encoding="utf-8")
            outside = Path(directory) / "outside"
            with self.assertRaisesRegex(LiftoverError, "inside the project root"):
                publish_audit(root, outside)
            existing = root / "audit"
            existing.mkdir()
            with self.assertRaisesRegex(LiftoverError, "must not already exist"):
                publish_audit(root, existing)


if __name__ == "__main__":
    unittest.main()
