from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from variant_liftover_explorer import controlled_demo
from variant_liftover_explorer.controlled_demo import build_controlled_fixtures, run_controlled_demo
from variant_liftover_explorer.liftover import (
    LiftoverError,
    ProjectionConfig,
    Variant,
    parse_paf,
    project_once,
    read_fasta,
    read_vcf,
    reverse_complement,
    validate_alignment_sequences,
)


class LiftoverUnitTests(unittest.TestCase):
    def test_reverse_complement_preserves_iupac_subset(self) -> None:
        self.assertEqual(reverse_complement("ACGTN"), "NACGT")

    def test_variant_rejects_nonliteral_and_equal_alleles(self) -> None:
        with self.assertRaises(LiftoverError):
            Variant("symbolic", "chr", 1, "A", "<DEL>")
        with self.assertRaises(LiftoverError):
            Variant("equal", "chr", 1, "A", "A")

    def test_parser_rejects_cigar_cs_consumption_disagreement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = build_controlled_fixtures(Path(directory) / "fixture")
            path = fixture["forward_paf"]
            text = path.read_text(encoding="utf-8").replace("5=1X14=", "4=1X14=", 1)
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(LiftoverError, "CIGAR and cs"):
                parse_paf(path)

    def test_cs_sequence_validation_detects_wrong_reference(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = build_controlled_fixtures(Path(directory) / "fixture")
            alignment = parse_paf(fixture["forward_paf"])[0]
            source = read_fasta(fixture["source_fasta"])
            target = read_fasta(fixture["target_fasta"])
            target["tgt_plus"] = "T" + target["tgt_plus"][1:]
            with self.assertRaisesRegex(LiftoverError, "cs target match"):
                validate_alignment_sequences(alignment, source, target)

    def test_project_once_rejects_source_ref_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = build_controlled_fixtures(Path(directory) / "fixture")
            source = read_fasta(fixture["source_fasta"])
            target = read_fasta(fixture["target_fasta"])
            alignments = parse_paf(fixture["forward_paf"])
            with self.assertRaisesRegex(LiftoverError, "disagrees with source FASTA"):
                project_once(
                    Variant("wrong", "src_plus", 2, "A", "G"),
                    source,
                    target,
                    alignments,
                    ProjectionConfig(minimum_alignment_length=10),
                )

    def test_gzipped_vcf_reader_rejects_duplicate_keys(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(
                    "##fileformat=VCFv4.3\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr\t2\tA\tC\tG\t.\tPASS\t.\n"
                    "chr\t2\tB\tC\tG\t.\tPASS\t.\n"
                )
            with self.assertRaisesRegex(LiftoverError, "duplicate variant key"):
                read_vcf(path)

    def test_atomic_demo_cleans_stage_after_failure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "controlled"
            with mock.patch.object(
                controlled_demo,
                "project_variants",
                side_effect=LiftoverError("injected failure"),
            ):
                with self.assertRaisesRegex(LiftoverError, "injected failure"):
                    run_controlled_demo(output)
            self.assertFalse(output.exists())
            self.assertEqual(list(Path(directory).glob(".controlled.staging-*")), [])


if __name__ == "__main__":
    unittest.main()
