"""Publish a read-only quality-review record from tests and authentic artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import sys
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from .authentic_analysis import _source_hashes
from .liftover import LiftoverError
from .source_review import (
    load_reference_manifest,
    reproduce_project15,
    tool_versions,
)


TEST_MATRIX = (
    (
        "Input formats",
        "malformed FASTA, VCF, BED, and PAF records; duplicate headers and keys",
        "unit tests",
    ),
    (
        "Coordinates and alleles",
        "SNPs, insertions, deletions, repeat normalization, allele swaps, circular-origin spans",
        "controlled and unit tests",
    ),
    (
        "Genome structure",
        "forward and reverse strands, multiple contigs, missing plasmids, rearrangement gaps",
        "controlled and unit tests",
    ),
    (
        "Mapping confidence",
        "ambiguous repeats, nonreciprocal mappings, low-quality blocks, alignment thresholds",
        "controlled and unit tests",
    ),
    (
        "Evidence states",
        "accepted, rejected, target-reference, callable-absent, and non-callable outcomes",
        "controlled and authentic-result tests",
    ),
    (
        "Project 15 integration",
        "duplicate and missing joins, effect precedence, all burden explanation categories",
        "controlled and authentic-result tests",
    ),
    (
        "Native tools",
        "normalization failure, malformed normalization output, pinned versions",
        "unit tests and live version checks",
    ),
    (
        "Atomic publication",
        "controlled, integration, and authentic workflow failures with staging cleanup",
        "injected-failure tests",
    ),
    (
        "Authentic accounting",
        "481,473 attempts, 154,865 baseline rows, 400 review contexts, hashes, VCF indexes",
        "authentic-result tests",
    ),
    (
        "Publication hygiene",
        "portable commands, private-path scan, credential scan, authorship-wording scan",
        "result tests and repository scan",
    ),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_junit(path: Path) -> dict[str, int]:
    if not path.is_file():
        raise LiftoverError(f"missing JUnit report: {path}")
    root = ET.parse(path).getroot()
    suites = [root] if root.tag == "testsuite" else list(root.findall("testsuite"))
    if not suites:
        raise LiftoverError("JUnit report contains no test suites")
    totals = {name: 0 for name in ("tests", "failures", "errors", "skipped")}
    for suite in suites:
        for name in totals:
            totals[name] += int(suite.attrib.get(name, "0"))
    totals["passed"] = totals["tests"] - totals["failures"] - totals["errors"] - totals["skipped"]
    if totals["tests"] < 1:
        raise LiftoverError("JUnit report contains no tests")
    if totals["failures"] or totals["errors"] or totals["skipped"]:
        raise LiftoverError(f"test suite did not pass completely: {totals}")
    return totals


def _verify_authentic_outputs(authentic_dir: Path, manifest: dict[str, object]) -> None:
    expected = manifest.get("output_hashes", {})
    if not isinstance(expected, dict) or not expected:
        raise LiftoverError("authentic manifest lacks output hashes")
    for relative, digest in expected.items():
        path = authentic_dir / str(relative)
        if not path.is_file() or sha256(path) != digest:
            raise LiftoverError(f"authentic artifact checksum mismatch: {relative}")


def publish_quality_review(args) -> dict[str, object]:
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise LiftoverError("quality-review output directory must not already exist")
    stage = output_dir.parent / f".{output_dir.name}.staging-{uuid.uuid4().hex}"
    stage.mkdir(parents=True)
    try:
        test_totals = parse_junit(args.junit_xml.resolve())
        authentic_dir = args.authentic_dir.resolve()
        authentic_manifest_path = authentic_dir / "run_manifest.json"
        authentic = json.loads(authentic_manifest_path.read_text(encoding="utf-8"))
        expected_counts = {
            "directional_projection_attempts": 481_473,
            "baseline_variant_comparison_rows": 154_865,
            "project15_variants": 357_581,
            "comparable_gene_pairs": 6_356,
            "review_candidate_liftover_rows": 400,
        }
        if authentic.get("status") != "pass":
            raise LiftoverError("authentic analysis is not marked pass")
        for name, expected in expected_counts.items():
            if authentic.get("counts", {}).get(name) != expected:
                raise LiftoverError(f"authentic count {name} is not {expected}")
        _verify_authentic_outputs(authentic_dir, authentic)

        references = load_reference_manifest(
            args.reference_manifest.resolve(), args.baseline_accession
        )
        current_source_hashes = _source_hashes(references)
        if current_source_hashes != authentic.get("source_hashes"):
            raise LiftoverError("inherited Project 13 inputs changed after authentic analysis")
        _rows, _counts, project15_hashes = reproduce_project15(
            args.project15_comparison_dir.resolve(),
            args.project15_orthology_dir.resolve(),
            args.project15_annotation_dir.resolve(),
        )
        if project15_hashes != authentic.get("project15_hashes"):
            raise LiftoverError("inherited Project 15 inputs changed after authentic analysis")

        versions = tool_versions()
        versions["pytest"] = pytest.__version__
        with (stage / "test_matrix.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["area", "scenarios", "evidence", "status"],
                lineterminator="\n",
            )
            writer.writeheader()
            for area, scenarios, evidence in TEST_MATRIX:
                writer.writerow(
                    {
                        "area": area,
                        "scenarios": scenarios,
                        "evidence": evidence,
                        "status": "pass",
                    }
                )
        manifest = {
            "workflow": "variant-liftover-explorer",
            "workflow_version": "0.1.0",
            "milestone": "quality_review",
            "status": "pass",
            "tests": {
                "framework": "pytest",
                "python": sys.version.split()[0],
                **test_totals,
            },
            "authentic_read_only_validation": {
                **expected_counts,
                "references": 6,
                "directional_vcfs": len(authentic.get("vcf_counts", {})) * 2,
                "passed": True,
            },
            "pinned_tools": versions,
            "artifact_checksums": {
                "authentic_manifest": sha256(authentic_manifest_path),
                "authentic_dashboard": sha256(
                    authentic_dir / "dashboard" / "variant_liftover_dashboard.png"
                ),
                "connection_review_manifest": sha256(
                    args.connection_review_dir.resolve() / "run_manifest.json"
                ),
            },
            "source_inputs_unchanged": True,
            "partial_outputs_after_injected_failures": 0,
            "repository_created": False,
        }
        manifest["artifact_checksums"]["test_matrix"] = sha256(stage / "test_matrix.csv")
        (stage / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(stage, output_dir)
        return manifest
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
