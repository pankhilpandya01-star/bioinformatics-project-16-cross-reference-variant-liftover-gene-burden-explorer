"""Deterministic controlled data and atomic publication for the liftover engine."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import uuid
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping

from .liftover import (
    FINAL_STATUSES,
    LiftoverError,
    ProjectionConfig,
    ProjectionResult,
    Variant,
    project_variants,
    reverse_complement,
    write_projection_csv,
)


EXPECTED_STATUSES = {
    "accepted_match": 4,
    "rejected_match": 1,
    "target_reference_match": 1,
    "callable_no_candidate": 1,
    "target_not_callable": 1,
    "ambiguous_mapping": 1,
    "nonreciprocal_mapping": 1,
    "nonalignable": 1,
    "complex_alignment_context": 1,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_fasta(path: Path, records: Mapping[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")


def _paf_row(
    query: str,
    query_length: int,
    strand: str,
    target: str,
    target_length: int,
    matches: int,
    alignment_bases: int,
    mapq: int,
    cigar: str,
    cs: str,
    alignment_type: str = "P",
    score: int | None = None,
) -> str:
    chain_score = alignment_bases * 2 if score is None else score
    return "\t".join(
        [
            query,
            str(query_length),
            "0",
            str(query_length),
            strand,
            target,
            str(target_length),
            "0",
            str(target_length),
            str(matches),
            str(alignment_bases),
            str(mapq),
            f"tp:A:{alignment_type}",
            f"s1:i:{chain_score}",
            "s2:i:0",
            f"cg:Z:{cigar}",
            f"cs:Z:{cs}",
        ]
    )


def build_controlled_fixtures(directory: Path) -> dict[str, object]:
    """Write small explicit genomes, PAFs, VCF-like variants, and callable BED data."""
    directory.mkdir(parents=True)
    target_reverse = "TTGACCGTAACTGGCATCGA"
    source_records = {
        "src_plus": "ACGTATGTACGTACGTACGT",
        "src_reverse": reverse_complement(target_reverse),
        "src_missing": "CCCCCTTTTT",
        "src_repeat": "GATTACAGATTACAGATTA",
        "src_gap": "AAAAAGCCCCC",
        "src_oneway": "AGCTTAGGCTAC",
    }
    target_records = {
        "tgt_plus": "ACGTACGTACGTACGTACGT",
        "tgt_reverse": target_reverse,
        "tgt_repeat_a": source_records["src_repeat"],
        "tgt_repeat_b": source_records["src_repeat"],
        "tgt_gap": "AAAAACCCCC",
        "tgt_oneway": source_records["src_oneway"],
    }
    source_fasta = directory / "source.fasta"
    target_fasta = directory / "target.fasta"
    forward_paf = directory / "source_to_target.paf"
    reverse_paf = directory / "target_to_source.paf"
    callable_bed = directory / "target_callable.bed"
    _write_fasta(source_fasta, source_records)
    _write_fasta(target_fasta, target_records)

    forward_rows = [
        _paf_row(
            "src_plus",
            20,
            "+",
            "tgt_plus",
            20,
            19,
            20,
            60,
            "5=1X14=",
            "=ACGTA*ct=GTACGTACGTACGT",
        ),
        _paf_row(
            "src_reverse",
            20,
            "-",
            "tgt_reverse",
            20,
            20,
            20,
            60,
            "20=",
            f"={target_reverse}",
        ),
        _paf_row(
            "src_repeat",
            19,
            "+",
            "tgt_repeat_a",
            19,
            19,
            19,
            0,
            "19=",
            f"={source_records['src_repeat']}",
            score=38,
        ),
        _paf_row(
            "src_repeat",
            19,
            "+",
            "tgt_repeat_b",
            19,
            19,
            19,
            0,
            "19=",
            f"={source_records['src_repeat']}",
            alignment_type="S",
            score=38,
        ),
        _paf_row(
            "src_gap",
            11,
            "+",
            "tgt_gap",
            10,
            10,
            11,
            60,
            "5=1I5=",
            "=AAAAA+g=CCCCC",
        ),
        _paf_row(
            "src_oneway",
            12,
            "+",
            "tgt_oneway",
            12,
            12,
            12,
            60,
            "12=",
            f"={source_records['src_oneway']}",
        ),
    ]
    reverse_rows = [
        _paf_row(
            "tgt_plus",
            20,
            "+",
            "src_plus",
            20,
            19,
            20,
            60,
            "5=1X14=",
            "=ACGTA*tc=GTACGTACGTACGT",
        ),
        _paf_row(
            "tgt_reverse",
            20,
            "-",
            "src_reverse",
            20,
            20,
            20,
            60,
            "20=",
            f"={source_records['src_reverse']}",
        ),
        _paf_row(
            "tgt_gap",
            10,
            "+",
            "src_gap",
            11,
            10,
            11,
            60,
            "5=1D5=",
            "=AAAAA-g=CCCCC",
        ),
    ]
    forward_paf.write_text("\n".join(forward_rows) + "\n", encoding="utf-8")
    reverse_paf.write_text("\n".join(reverse_rows) + "\n", encoding="utf-8")
    callable_bed.write_text(
        "tgt_plus\t0\t7\n"
        "tgt_plus\t8\t20\n"
        "tgt_reverse\t0\t20\n"
        "tgt_oneway\t0\t12\n",
        encoding="utf-8",
    )

    variants = [
        Variant("accepted_snp", "src_plus", 2, "C", "G"),
        Variant("rejected_snp", "src_plus", 3, "G", "A"),
        Variant("target_reference", "src_plus", 6, "T", "C"),
        Variant("callable_absent", "src_plus", 7, "G", "A"),
        Variant("not_callable", "src_plus", 8, "T", "A"),
        Variant("accepted_insertion", "src_plus", 10, "C", "CAA"),
        Variant("accepted_deletion", "src_plus", 13, "AC", "A"),
        Variant("reverse_snp", "src_reverse", 4, "A", "G"),
        Variant("missing_contig_alignment", "src_missing", 2, "C", "A"),
        Variant("repeat_ambiguity", "src_repeat", 4, "T", "C"),
        Variant("gap_context", "src_gap", 5, "A", "AG"),
        Variant("one_way_only", "src_oneway", 2, "G", "A"),
    ]
    accepted = {
        "tgt_plus:2:C:G",
        "tgt_plus:10:C:CAA",
        "tgt_plus:13:AC:A",
        "tgt_reverse:17:T:C",
    }
    rejected = {"tgt_plus:3:G:A"}
    callable_intervals = {
        "tgt_plus": [(0, 7), (8, 20)],
        "tgt_reverse": [(0, 20)],
        "tgt_oneway": [(0, 12)],
    }
    return {
        "source_fasta": source_fasta,
        "target_fasta": target_fasta,
        "forward_paf": forward_paf,
        "reverse_paf": reverse_paf,
        "callable_bed": callable_bed,
        "variants": variants,
        "accepted": accepted,
        "rejected": rejected,
        "callable_intervals": callable_intervals,
    }


def _run(command: list[str]) -> str:
    try:
        result = subprocess.run(command, capture_output=True, text=True)
    except OSError as error:
        raise LiftoverError(f"command could not start ({command[0]}): {error}") from error
    if result.returncode != 0:
        raise LiftoverError(f"command failed ({command[0]}): {result.stderr.strip()}")
    return (result.stdout or result.stderr).strip()


def _write_status_summary(path: Path, counts: Counter[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["status", "count"], lineterminator="\n")
        writer.writeheader()
        for status in sorted(FINAL_STATUSES):
            writer.writerow({"status": status, "count": counts.get(status, 0)})


def _write_projection_vcf(
    path: Path,
    results: Iterable[ProjectionResult],
    target_sequences: Mapping[str, str],
) -> int:
    included = [
        result
        for result in results
        if result.target_position is not None
        and result.status
        in {"accepted_match", "rejected_match", "callable_no_candidate", "target_not_callable"}
    ]
    contig_order = {contig: index for index, contig in enumerate(target_sequences)}
    included.sort(key=lambda item: (contig_order[item.target_contig], int(item.target_position)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("##fileformat=VCFv4.3\n")
        handle.write('##INFO=<ID=LIFTOVER_STATUS,Number=1,Type=String,Description="Controlled projection status">\n')
        for contig, sequence in target_sequences.items():
            handle.write(f"##contig=<ID={contig},length={len(sequence)}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for result in included:
            handle.write(
                f"{result.target_contig}\t{result.target_position}\t{result.source.identifier}\t"
                f"{result.target_reference}\t{result.target_alternate}\t.\tPASS\t"
                f"LIFTOVER_STATUS={result.status}\n"
            )
    return len(included)


def run_controlled_demo(
    output_dir: Path,
    *,
    bcftools: str = "bcftools",
    samtools: str = "samtools",
) -> dict[str, object]:
    """Run the controlled liftover suite and publish the complete directory atomically."""
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise LiftoverError("output directory must not already exist")
    stage = output_dir.parent / f".{output_dir.name}.staging-{uuid.uuid4().hex}"
    fixtures = stage / "fixtures"
    tables = stage / "tables"
    reports = stage / "reports"
    work = stage / "work"
    stage.mkdir(parents=True)
    tables.mkdir()
    reports.mkdir()
    work.mkdir()
    try:
        fixture = build_controlled_fixtures(fixtures)
        source_fasta = fixture["source_fasta"]
        target_fasta = fixture["target_fasta"]
        _run([samtools, "faidx", str(source_fasta)])
        _run([samtools, "faidx", str(target_fasta)])
        config = ProjectionConfig(
            minimum_alignment_length=10,
            minimum_alignment_identity=0.85,
            minimum_alignment_mapq=20,
            ambiguity_score_fraction=0.95,
        )
        results = project_variants(
            fixture["variants"],
            source_fasta,
            target_fasta,
            fixture["forward_paf"],
            fixture["reverse_paf"],
            fixture["accepted"],
            fixture["rejected"],
            fixture["callable_intervals"],
            config,
            work,
        )
        counts = Counter(result.status for result in results)
        observed = {status: counts.get(status, 0) for status in EXPECTED_STATUSES}
        if observed != EXPECTED_STATUSES or sum(counts.values()) != 12:
            raise LiftoverError(f"controlled status accounting failed: {dict(counts)}")
        write_projection_csv(tables / "controlled_projections.csv", results)
        _write_status_summary(tables / "status_summary.csv", counts)

        from .liftover import read_fasta

        target_sequences = read_fasta(target_fasta)
        vcf_path = reports / "normalized_projected.vcf"
        vcf_count = _write_projection_vcf(vcf_path, results, target_sequences)
        compressed_vcf = reports / "normalized_projected.vcf.gz"
        _run(
            [
                bcftools,
                "view",
                "--no-version",
                "-Oz",
                "-o",
                str(compressed_vcf),
                str(vcf_path),
            ]
        )
        _run([bcftools, "index", "-f", "--csi", str(compressed_vcf)])
        _run([bcftools, "view", "-h", str(compressed_vcf)])

        versions = {
            "bcftools": _run([bcftools, "--version"]).splitlines()[0],
            "samtools": _run([samtools, "--version"]).splitlines()[0],
        }
        with (reports / "tool_versions.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["tool", "version"], lineterminator="\n")
            writer.writeheader()
            writer.writerows({"tool": key, "version": value} for key, value in versions.items())

        source_hashes = {
            path.name: _sha256(path)
            for path in (
                source_fasta,
                target_fasta,
                fixture["forward_paf"],
                fixture["reverse_paf"],
                fixture["callable_bed"],
            )
        }
        manifest = {
            "workflow": "variant-liftover-explorer",
            "workflow_version": "0.1.0",
            "milestone": "controlled_liftover_engine",
            "status": "pass",
            "parameters": {
                "minimum_alignment_length": config.minimum_alignment_length,
                "minimum_alignment_identity": config.minimum_alignment_identity,
                "minimum_alignment_mapq": config.minimum_alignment_mapq,
                "ambiguity_score_fraction": config.ambiguity_score_fraction,
            },
            "counts": {
                "source_variants": len(results),
                "projected_vcf_records": vcf_count,
                "reciprocal": sum(result.reciprocal for result in results),
                "normalized": sum(result.normalized for result in results),
                "statuses": observed,
            },
            "tools": versions,
            "fixture_checksums": source_hashes,
            "atomic_publication": True,
            "interpretation_boundary": (
                "Controlled examples validate coordinate and allele mechanics; they are not biological findings."
            ),
        }
        (stage / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (stage / "portable_commands.txt").write_text(
            "samtools faidx <SOURCE_FASTA>\n"
            "samtools faidx <TARGET_FASTA>\n"
            "bcftools norm --no-version -f <TARGET_FASTA> -m -any -Ov -o <NORMALIZED_VCF> <PROVISIONAL_VCF>\n"
            "bcftools view --no-version -Oz -o <COMPRESSED_VCF> <NORMALIZED_VCF>\n"
            "bcftools index -f --csi <COMPRESSED_VCF>\n",
            encoding="utf-8",
        )
        shutil.rmtree(work)
        os.replace(stage, output_dir)
        return manifest
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
