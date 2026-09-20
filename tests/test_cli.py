from __future__ import annotations

import json
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from unittest import mock

from variant_liftover_explorer.cli import build_parser, main
from variant_liftover_explorer.source_review import ReviewError
from variant_liftover_explorer.workflow import WorkflowConfig, WorkflowError, run_workflow


class CliTests(unittest.TestCase):
    def test_parser_exposes_the_approved_defaults(self) -> None:
        args = build_parser().parse_args(
            [
                "--reference-manifest",
                "references.csv",
                "--project15-comparison-dir",
                "comparison",
                "--project15-orthology-dir",
                "orthology",
                "--project15-annotation-dir",
                "annotation",
                "--output-dir",
                "output",
            ]
        )
        self.assertEqual(args.baseline_accession, "GCF_000005845.2")
        self.assertEqual(args.sample_id, "SRR13921545")
        self.assertEqual(args.alignment_preset, "asm20")
        self.assertEqual(args.minimum_alignment_length, 10_000)
        self.assertEqual(args.minimum_alignment_identity, 0.85)
        self.assertEqual(args.minimum_alignment_mapq, 20)
        self.assertEqual(args.ambiguity_score_fraction, 0.95)
        self.assertEqual(args.minimum_comparable_coverage, 0.70)
        self.assertEqual(args.threads, 1)

    def test_version_and_workflow_errors_have_stable_exit_codes(self) -> None:
        stdout = StringIO()
        with redirect_stdout(stdout), self.assertRaises(SystemExit) as raised:
            main(["--version"])
        self.assertEqual(raised.exception.code, 0)
        self.assertIn("variant-liftover-explorer 0.1.0", stdout.getvalue())

        required = [
            "--reference-manifest", "references.csv",
            "--project15-comparison-dir", "comparison",
            "--project15-orthology-dir", "orthology",
            "--project15-annotation-dir", "annotation",
            "--output-dir", "output",
        ]
        stderr = StringIO()
        with mock.patch(
            "variant_liftover_explorer.cli.run_workflow",
            side_effect=WorkflowError("controlled failure"),
        ), redirect_stderr(stderr):
            self.assertEqual(main(required), 2)
        self.assertIn("controlled failure", stderr.getvalue())

    def test_configuration_rejects_invalid_values_and_existing_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = dict(
                reference_manifest=root / "references.csv",
                project15_comparison_dir=root / "comparison",
                project15_orthology_dir=root / "orthology",
                project15_annotation_dir=root / "annotation",
                output_dir=root / "output",
            )
            invalid = (
                WorkflowConfig(**base, alignment_preset="reads"),
                WorkflowConfig(**base, sample_id="unsafe sample"),
                WorkflowConfig(**base, threads=0),
                WorkflowConfig(**base, minimum_comparable_coverage=0),
                WorkflowConfig(**base, minimum_alignment_mapq=256),
            )
            for config in invalid:
                with self.subTest(config=config), self.assertRaises(WorkflowError):
                    config.validate()
            base["output_dir"].mkdir()
            with self.assertRaisesRegex(WorkflowError, "must not already exist"):
                WorkflowConfig(**base).validate()

    def test_public_workflow_publishes_one_complete_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference_manifest = root / "references.csv"
            reference_manifest.write_text("controlled\n", encoding="utf-8")
            output = root / "output"
            config = WorkflowConfig(
                reference_manifest=reference_manifest,
                project15_comparison_dir=root / "comparison",
                project15_orthology_dir=root / "orthology",
                project15_annotation_dir=root / "annotation",
                output_dir=output,
            )

            def fake_connection(args) -> None:
                (args.alignment_dir / "alignments").mkdir(parents=True)
                (args.alignment_dir / "logs").mkdir()
                (args.alignment_dir / "alignments" / "one.paf").write_text(
                    "controlled\n", encoding="utf-8"
                )
                args.output_dir.mkdir()
                (args.output_dir / "run_manifest.json").write_text(
                    json.dumps(
                        {
                            "alignment_counts": {"directions": 10},
                            "source_inputs_unchanged": True,
                        }
                    ),
                    encoding="utf-8",
                )

            def fake_authentic(args) -> dict[str, object]:
                args.output_dir.mkdir()
                (args.output_dir / "dashboard").mkdir()
                (args.output_dir / "dashboard" / "variant_liftover_dashboard.png").write_bytes(
                    b"controlled-dashboard"
                )
                manifest = {
                    "counts": {"directional_projection_attempts": 481_473},
                    "source_inputs_unchanged": True,
                }
                (args.output_dir / "run_manifest.json").write_text(
                    json.dumps(manifest), encoding="utf-8"
                )
                return manifest

            with (
                mock.patch(
                    "variant_liftover_explorer.workflow.connection_review",
                    side_effect=fake_connection,
                ),
                mock.patch(
                    "variant_liftover_explorer.workflow.authentic_analysis",
                    side_effect=fake_authentic,
                ),
            ):
                self.assertEqual(run_workflow(config), output.resolve())
            self.assertTrue((output / "alignments" / "one.paf").is_file())
            self.assertTrue((output / "alignment_logs").is_dir())
            self.assertTrue((output / "connection_review" / "run_manifest.json").is_file())
            self.assertTrue((output / "analysis" / "run_manifest.json").is_file())
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "pass")
            self.assertTrue(manifest["source_inputs_unchanged"])

    def test_public_workflow_failure_removes_all_staged_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference_manifest = root / "references.csv"
            reference_manifest.write_text("controlled\n", encoding="utf-8")
            output = root / "output"
            config = WorkflowConfig(
                reference_manifest=reference_manifest,
                project15_comparison_dir=root / "comparison",
                project15_orthology_dir=root / "orthology",
                project15_annotation_dir=root / "annotation",
                output_dir=output,
            )

            def fail_after_writing(args) -> None:
                args.output_dir.mkdir()
                args.alignment_dir.mkdir()
                raise ReviewError("controlled connection failure")

            with mock.patch(
                "variant_liftover_explorer.workflow.connection_review",
                side_effect=fail_after_writing,
            ), self.assertRaisesRegex(WorkflowError, "controlled connection failure"):
                run_workflow(config)
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".output.staging-*")), [])


if __name__ == "__main__":
    unittest.main()
