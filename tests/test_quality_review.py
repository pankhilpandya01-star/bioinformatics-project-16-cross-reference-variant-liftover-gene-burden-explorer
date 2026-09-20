from __future__ import annotations

import csv
import hashlib
import subprocess
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from variant_liftover_explorer import authentic_analysis as authentic_module
from variant_liftover_explorer.controlled_demo import build_controlled_fixtures, run_controlled_demo
from variant_liftover_explorer.liftover import (
    LiftoverError,
    ProjectionConfig,
    ProjectionDraft,
    Variant,
    normalize_drafts,
    parse_paf,
    project_once,
    project_variants_loaded,
    read_bed,
    read_fasta,
    read_vcf,
    validate_alignment_sequences,
)
from variant_liftover_explorer.quality_review import parse_junit, publish_quality_review
from variant_liftover_explorer.source_review import ReferenceInput


ROOT = Path(__file__).resolve().parents[1]
CONTROLLED_RESULTS = ROOT / "results" / "controlled_core" / "tables" / "controlled_projections.csv"


def controlled_rows() -> dict[str, dict[str, str]]:
    with CONTROLLED_RESULTS.open(newline="", encoding="utf-8") as handle:
        return {row["source_variant_id"]: row for row in csv.DictReader(handle)}


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class InputValidationQualityTests(unittest.TestCase):
    def test_projection_parameters_reject_unsafe_boundaries(self) -> None:
        invalid = (
            ProjectionConfig(minimum_alignment_length=0),
            ProjectionConfig(minimum_alignment_identity=0),
            ProjectionConfig(minimum_alignment_identity=1.01),
            ProjectionConfig(minimum_alignment_mapq=-1),
            ProjectionConfig(minimum_alignment_mapq=256),
            ProjectionConfig(ambiguity_score_fraction=0),
            ProjectionConfig(ambiguity_score_fraction=1.01),
        )
        for config in invalid:
            with self.subTest(config=config), self.assertRaises(LiftoverError):
                config.validate()

    def test_malformed_fasta_records_are_rejected(self) -> None:
        cases = {
            "sequence_before_header": "ACGT\n",
            "invalid_base": ">chr\nACGU\n",
            "empty_record": ">chr\n",
            "duplicate_id": ">chr\nACGT\n>chr\nACGT\n",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, content in cases.items():
                with self.subTest(name=name):
                    path = root / f"{name}.fasta"
                    path.write_text(content, encoding="ascii")
                    with self.assertRaises(LiftoverError):
                        read_fasta(path)

    def test_malformed_vcf_records_are_rejected(self) -> None:
        header = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        cases = {
            "missing_header": "chr\t1\t.\tA\tC\t.\tPASS\t.\n",
            "multiallelic": header + "chr\t1\t.\tA\tC,G\t.\tPASS\t.\n",
            "symbolic": header + "chr\t1\t.\tA\t<DEL>\t.\tPASS\t.\n",
            "bad_position": header + "chr\tnot-an-integer\t.\tA\tC\t.\tPASS\t.\n",
            "zero_position": header + "chr\t0\t.\tA\tC\t.\tPASS\t.\n",
            "equal_alleles": header + "chr\t1\t.\tA\tA\t.\tPASS\t.\n",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, content in cases.items():
                with self.subTest(name=name):
                    path = root / f"{name}.vcf"
                    path.write_text(content, encoding="utf-8")
                    with self.assertRaises(LiftoverError):
                        read_vcf(path)

    def test_empty_vcf_with_a_valid_header_is_supported(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty.vcf"
            path.write_text(
                "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                encoding="utf-8",
            )
            self.assertEqual(read_vcf(path), [])

    def test_malformed_bed_records_are_rejected(self) -> None:
        cases = {
            "missing_column": "chr\t0\n",
            "not_integer": "chr\tA\t5\n",
            "negative": "chr\t-1\t5\n",
            "empty": "chr\t5\t5\n",
            "empty_contig": "\t0\t5\n",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, content in cases.items():
                with self.subTest(name=name):
                    path = root / f"{name}.bed"
                    path.write_text(content, encoding="utf-8")
                    with self.assertRaises(LiftoverError):
                        read_bed(path)

    def test_malformed_paf_numeric_and_tag_fields_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = build_controlled_fixtures(Path(directory) / "fixture")
            original = fixture["forward_paf"].read_text(encoding="utf-8").splitlines()[0]
            cases: dict[str, str] = {}
            fields = original.split("\t")
            bad_integer = fields.copy()
            bad_integer[1] = "length"
            cases["non_integer"] = "\t".join(bad_integer)
            bad_mapq = fields.copy()
            bad_mapq[11] = "300"
            cases["mapq"] = "\t".join(bad_mapq)
            bad_counts = fields.copy()
            bad_counts[9] = "21"
            cases["counts"] = "\t".join(bad_counts)
            cases["tp"] = original.replace("tp:A:P", "tp:A:X", 1)
            for name, row in cases.items():
                with self.subTest(name=name):
                    path = Path(directory) / f"{name}.paf"
                    path.write_text(row + "\n", encoding="utf-8")
                    with self.assertRaises(LiftoverError):
                        parse_paf(path)

    def test_alignment_contig_and_length_disagreements_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = build_controlled_fixtures(Path(directory) / "fixture")
            alignment = parse_paf(fixture["forward_paf"])[0]
            source = read_fasta(fixture["source_fasta"])
            target = read_fasta(fixture["target_fasta"])
            with self.assertRaisesRegex(LiftoverError, "unknown contig"):
                validate_alignment_sequences(
                    replace(alignment, target_name="missing_plasmid"), source, target
                )
            with self.assertRaisesRegex(LiftoverError, "length disagrees"):
                validate_alignment_sequences(
                    replace(alignment, target_length=alignment.target_length + 1), source, target
                )


class ProjectionBoundaryQualityTests(unittest.TestCase):
    def test_controlled_edge_outcomes_remain_explicit(self) -> None:
        records = controlled_rows()
        expected = {
            "accepted_snp": "accepted_match",
            "accepted_insertion": "accepted_match",
            "accepted_deletion": "accepted_match",
            "reverse_snp": "accepted_match",
            "target_reference": "target_reference_match",
            "rejected_snp": "rejected_match",
            "callable_absent": "callable_no_candidate",
            "not_callable": "target_not_callable",
            "repeat_ambiguity": "ambiguous_mapping",
            "one_way_only": "nonreciprocal_mapping",
            "missing_contig_alignment": "nonalignable",
            "gap_context": "complex_alignment_context",
        }
        self.assertEqual({name: records[name]["status"] for name in expected}, expected)
        self.assertEqual(records["reverse_snp"]["strand"], "-")
        self.assertEqual(records["reverse_snp"]["target_contig"], "tgt_reverse")
        self.assertEqual(records["target_reference"]["reciprocal"], "true")

    def test_rearrangement_gap_is_not_forced_through(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = build_controlled_fixtures(Path(directory) / "fixture")
            result = project_once(
                Variant("gap_base", "src_gap", 6, "G", "A"),
                read_fasta(fixture["source_fasta"]),
                read_fasta(fixture["target_fasta"]),
                parse_paf(fixture["forward_paf"]),
                ProjectionConfig(minimum_alignment_length=10),
            )
            self.assertEqual(result.status, "complex_alignment_context")

    def test_circular_origin_spanning_allele_is_rejected_for_presplitting(self) -> None:
        with self.assertRaisesRegex(LiftoverError, "circular-origin"):
            project_once(
                Variant("wrap", "circular", 4, "AA", "A"),
                {"circular": "AAAA"},
                {"target": "AAAA"},
                [],
                ProjectionConfig(minimum_alignment_length=1),
            )

    def test_duplicate_evidence_and_unsorted_callability_are_rejected(self) -> None:
        variant = Variant("v", "chr", 1, "A", "C")
        common = {
            "source_sequences": {"chr": "A"},
            "target_sequences": {"chr": "A"},
            "forward": [],
            "reverse": [],
            "config": ProjectionConfig(minimum_alignment_length=1),
            "source_fasta": Path("source.fasta"),
            "target_fasta": Path("target.fasta"),
            "work_dir": Path("work"),
        }
        with self.assertRaisesRegex(LiftoverError, "duplicate normalized keys"):
            project_variants_loaded(
                [variant, variant],
                accepted_target_keys=set(),
                rejected_target_keys=set(),
                callable_intervals={},
                **common,
            )
        with self.assertRaisesRegex(LiftoverError, "evidence overlap"):
            project_variants_loaded(
                [],
                accepted_target_keys={"chr:1:A:C"},
                rejected_target_keys={"chr:1:A:C"},
                callable_intervals={},
                **common,
            )
        with self.assertRaisesRegex(LiftoverError, "invalid or unsorted"):
            project_variants_loaded(
                [],
                accepted_target_keys=set(),
                rejected_target_keys=set(),
                callable_intervals={"chr": [(5, 10), (0, 4)]},
                **common,
            )

    def test_repeat_indel_is_left_normalized_by_bcftools(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "target.fasta"
            fasta.write_text(">chr\nCAAAAAGT\n", encoding="ascii")
            subprocess.run(["samtools", "faidx", str(fasta)], check=True, capture_output=True)
            source = Variant("repeat_insertion", "source", 5, "A", "AA")
            draft = ProjectionDraft(
                source=source,
                status="projected",
                reason="controlled repeat",
                target_contig="chr",
                target_position=5,
                target_reference="A",
                target_alternate="AA",
            )
            normalized = normalize_drafts([draft], fasta, root / "work")
            self.assertEqual(normalized["P1"].key, "chr:1:C:CA")


class FailureAndCleanupQualityTests(unittest.TestCase):
    def test_normalization_subprocess_failure_is_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "target.fasta"
            fasta.write_text(">chr\nACGT\n", encoding="ascii")
            draft = ProjectionDraft(
                source=Variant("v", "source", 1, "A", "C"),
                status="projected",
                reason="controlled",
                target_contig="chr",
                target_position=1,
                target_reference="A",
                target_alternate="C",
            )
            failed = SimpleNamespace(returncode=2, stderr="controlled failure", stdout="")
            with mock.patch(
                "variant_liftover_explorer.liftover.subprocess.run", return_value=failed
            ), self.assertRaisesRegex(LiftoverError, "controlled failure"):
                normalize_drafts([draft], fasta, root / "work")

    def test_native_failure_leaves_no_controlled_output_or_stage(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "controlled"
            with self.assertRaises(LiftoverError):
                run_controlled_demo(output, samtools="definitely-missing-samtools")
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".controlled.staging-*")), [])

    def test_authentic_failure_removes_output_and_work_stages(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            references: list[ReferenceInput] = []
            accepted_counts = [30_973, 92_665, 125_327, 49_201, 15_264, 44_151]
            inventories: dict[str, dict[str, int]] = {}
            for index, accepted_count in enumerate(accepted_counts):
                accession = f"REF{index}"
                fasta = root / f"{accession}.fasta"
                fasta.write_text(f">{accession}\nA\n", encoding="ascii")
                references.append(
                    ReferenceInput(
                        accession=accession,
                        label=accession,
                        role="baseline" if index == 0 else "alternative",
                        reference_fasta=fasta,
                        accepted_vcf=fasta,
                        evaluated_vcf=fasta,
                        rejected_vcf=fasta,
                        callable_bed=fasta,
                        expected_hashes={"reference_fasta": file_hash(fasta)},
                    )
                )
                inventories[accession] = {"accepted_variants": accepted_count}
            args = SimpleNamespace(
                output_dir=root / "output",
                work_dir=root / "work",
                minimum_alignment_length=10_000,
                minimum_alignment_identity=0.85,
                minimum_alignment_mapq=20,
                ambiguity_score_fraction=0.95,
                minimum_comparable_coverage=0.70,
                reference_manifest=root / "manifest.csv",
                baseline_accession="REF0",
                sample_id="controlled",
                project15_comparison_dir=root / "comparison",
                project15_orthology_dir=root / "orthology",
                project15_annotation_dir=root / "annotation",
                alignment_dir=root / "alignments",
                connection_review_dir=root / "connection",
                alignment_preset="asm20",
                threads=1,
            )
            with (
                mock.patch.object(authentic_module, "load_reference_manifest", return_value=references),
                mock.patch.object(authentic_module, "tool_versions", return_value={}),
                mock.patch.object(
                    authentic_module,
                    "verify_reference_inputs",
                    return_value=([], inventories, {}),
                ),
                mock.patch.object(
                    authentic_module,
                    "reproduce_project15",
                    return_value=([], {}, {}),
                ),
                mock.patch.object(authentic_module, "_source_hashes", return_value={}),
                mock.patch.object(
                    authentic_module,
                    "_run",
                    side_effect=LiftoverError("controlled native failure"),
                ),
                self.assertRaisesRegex(LiftoverError, "controlled native failure"),
            ):
                authentic_module.authentic_analysis(args)
            self.assertFalse(args.output_dir.exists())
            self.assertFalse(args.work_dir.exists())
            self.assertEqual(list(root.glob(".output.staging-*")), [])
            self.assertEqual(list(root.glob(".work.staging-*")), [])

    def test_quality_review_failure_removes_its_stage(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = SimpleNamespace(output_dir=root / "quality", junit_xml=root / "missing.xml")
            with self.assertRaisesRegex(LiftoverError, "missing JUnit"):
                publish_quality_review(args)
            self.assertFalse(args.output_dir.exists())
            self.assertEqual(list(root.glob(".quality.staging-*")), [])

    def test_junit_parser_rejects_failures_and_counts_passes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            passing = root / "passing.xml"
            passing.write_text(
                '<testsuites><testsuite tests="3" failures="0" errors="0" skipped="0"/></testsuites>',
                encoding="utf-8",
            )
            self.assertEqual(parse_junit(passing)["passed"], 3)
            failing = root / "failing.xml"
            failing.write_text(
                '<testsuite tests="2" failures="1" errors="0" skipped="0"/>',
                encoding="utf-8",
            )
            with self.assertRaisesRegex(LiftoverError, "did not pass completely"):
                parse_junit(failing)


if __name__ == "__main__":
    unittest.main()
