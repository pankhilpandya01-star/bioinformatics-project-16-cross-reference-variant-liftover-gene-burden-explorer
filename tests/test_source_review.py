from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from variant_liftover_explorer.source_review import (
    ReviewError,
    _merged_bases,
    fasta_lengths,
    read_vcf_keys,
)


class SourceReviewUnitTests(unittest.TestCase):
    def test_merged_bases_deduplicates_overlaps_per_contig(self) -> None:
        self.assertEqual(
            _merged_bases([("chr", 0, 10), ("chr", 5, 20), ("plasmid", 0, 4)]),
            24,
        )

    def test_fasta_lengths_rejects_duplicate_identifiers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference.fasta"
            path.write_text(">chr\nACGT\n>chr\nACGT\n", encoding="ascii")
            with self.assertRaises(ReviewError):
                fasta_lengths(path)

    def test_vcf_partition_reader_rejects_duplicate_keys(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.vcf"
            path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr\t2\t.\tA\tG\t50\tPASS\t.\tGT\t1\n"
                "chr\t2\t.\tA\tG\t50\tPASS\t.\tGT\t1\n",
                encoding="utf-8",
            )
            with self.assertRaises(ReviewError):
                read_vcf_keys(path, "SAMPLE")


if __name__ == "__main__":
    unittest.main()

