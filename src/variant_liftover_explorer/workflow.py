"""End-to-end public workflow for Project 16."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import uuid
from dataclasses import dataclass
from pathlib import Path

from . import __version__
from .authentic_analysis import authentic_analysis
from .liftover import LiftoverError, ProjectionConfig
from .source_review import ReviewError, connection_review, sha256


SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9_.-]+$")


class WorkflowError(RuntimeError):
    """Raised when the public workflow cannot publish a complete result."""


@dataclass(frozen=True)
class WorkflowConfig:
    reference_manifest: Path
    project15_comparison_dir: Path
    project15_orthology_dir: Path
    project15_annotation_dir: Path
    output_dir: Path
    baseline_accession: str = "GCF_000005845.2"
    sample_id: str = "SRR13921545"
    alignment_preset: str = "asm20"
    minimum_alignment_length: int = 10_000
    minimum_alignment_identity: float = 0.85
    minimum_alignment_mapq: int = 20
    ambiguity_score_fraction: float = 0.95
    minimum_comparable_coverage: float = 0.70
    threads: int = 1

    def validate(self) -> None:
        if self.output_dir.resolve().exists():
            raise WorkflowError("output directory must not already exist")
        if not SAFE_IDENTIFIER.fullmatch(self.baseline_accession):
            raise WorkflowError("baseline accession contains unsafe characters")
        if not SAFE_IDENTIFIER.fullmatch(self.sample_id):
            raise WorkflowError("sample ID contains unsafe characters")
        if self.alignment_preset not in {"asm5", "asm10", "asm20"}:
            raise WorkflowError("alignment preset must be asm5, asm10, or asm20")
        if self.threads < 1:
            raise WorkflowError("threads must be at least one")
        if not 0 < self.minimum_comparable_coverage <= 1:
            raise WorkflowError("minimum comparable coverage must be in (0, 1]")
        try:
            ProjectionConfig(
                minimum_alignment_length=self.minimum_alignment_length,
                minimum_alignment_identity=self.minimum_alignment_identity,
                minimum_alignment_mapq=self.minimum_alignment_mapq,
                ambiguity_score_fraction=self.ambiguity_score_fraction,
            ).validate()
        except LiftoverError as error:
            raise WorkflowError(str(error)) from error


def run_workflow(config: WorkflowConfig) -> Path:
    """Run alignment and authentic analysis, then publish one atomic directory."""
    config.validate()
    output_dir = config.output_dir.resolve()
    stage = output_dir.parent / f".{output_dir.name}.staging-{uuid.uuid4().hex}"
    stage.mkdir(parents=True)
    connection_dir = stage / "connection_review"
    native_dir = stage / "native_alignment_work"
    analysis_dir = stage / "analysis"
    work_dir = stage / "normalization_work"
    try:
        shared = {
            "reference_manifest": config.reference_manifest.resolve(),
            "project15_comparison_dir": config.project15_comparison_dir.resolve(),
            "project15_orthology_dir": config.project15_orthology_dir.resolve(),
            "project15_annotation_dir": config.project15_annotation_dir.resolve(),
            "baseline_accession": config.baseline_accession,
            "sample_id": config.sample_id,
            "alignment_preset": config.alignment_preset,
            "minimum_alignment_length": config.minimum_alignment_length,
            "minimum_alignment_identity": config.minimum_alignment_identity,
            "minimum_alignment_mapq": config.minimum_alignment_mapq,
            "minimum_comparable_coverage": config.minimum_comparable_coverage,
            "threads": config.threads,
        }
        connection_review(
            argparse.Namespace(
                **shared,
                alignment_dir=native_dir,
                output_dir=connection_dir,
            )
        )
        authentic_manifest = authentic_analysis(
            argparse.Namespace(
                **shared,
                ambiguity_score_fraction=config.ambiguity_score_fraction,
                connection_review_dir=connection_dir,
                alignment_dir=native_dir / "alignments",
                work_dir=work_dir,
                output_dir=analysis_dir,
            )
        )
        connection_manifest_path = connection_dir / "run_manifest.json"
        connection_manifest = json.loads(connection_manifest_path.read_text(encoding="utf-8"))

        # Native PAFs and logs are kept for a full user run, but the portfolio
        # repository may exclude them and retain their recorded checksums.
        os.replace(native_dir / "alignments", stage / "alignments")
        os.replace(native_dir / "logs", stage / "alignment_logs")
        native_dir.rmdir()
        manifest = {
            "workflow": "variant-liftover-explorer",
            "workflow_version": __version__,
            "status": "pass",
            "baseline_accession": config.baseline_accession,
            "sample_id": config.sample_id,
            "parameters": {
                "alignment_preset": config.alignment_preset,
                "minimum_alignment_length": config.minimum_alignment_length,
                "minimum_alignment_identity": config.minimum_alignment_identity,
                "minimum_alignment_mapq": config.minimum_alignment_mapq,
                "ambiguity_score_fraction": config.ambiguity_score_fraction,
                "minimum_comparable_coverage": config.minimum_comparable_coverage,
                "threads": config.threads,
            },
            "counts": authentic_manifest["counts"],
            "alignment_counts": connection_manifest["alignment_counts"],
            "source_inputs_unchanged": (
                authentic_manifest["source_inputs_unchanged"]
                and connection_manifest["source_inputs_unchanged"]
            ),
            "artifacts": {
                "connection_review": "connection_review",
                "alignments": "alignments",
                "alignment_logs": "alignment_logs",
                "analysis": "analysis",
            },
            "checksums": {
                "reference_manifest": sha256(config.reference_manifest.resolve()),
                "connection_manifest": sha256(connection_manifest_path),
                "analysis_manifest": sha256(analysis_dir / "run_manifest.json"),
                "dashboard": sha256(
                    analysis_dir / "dashboard" / "variant_liftover_dashboard.png"
                ),
            },
            "interpretation_boundary": (
                "Exact normalized matches require reciprocal coordinate and allele validation; "
                "unresolved projections are not biological absence."
            ),
        }
        (stage / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(stage, output_dir)
        return output_dir
    except (OSError, ValueError, json.JSONDecodeError, LiftoverError, ReviewError) as error:
        shutil.rmtree(stage, ignore_errors=True)
        raise WorkflowError(str(error)) from error
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
