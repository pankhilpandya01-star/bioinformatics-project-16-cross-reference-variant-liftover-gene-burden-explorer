"""Milestone 4 authentic cross-reference projection and Project 15 explanation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import uuid
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .integration import (
    BURDEN_CATEGORIES,
    ProjectionEvidence,
    connect_review_candidates,
    integrate_gene_burdens,
    load_project15_gene_pairs,
    load_project15_gene_variants,
    load_project15_review_contexts,
    validate_project15_burdens,
    write_rows,
)
from .liftover import (
    FINAL_STATUSES,
    LiftoverError,
    ProjectionConfig,
    ProjectionResult,
    Variant,
    parse_paf,
    project_variants_loaded,
    read_bed,
    read_fasta,
    read_vcf,
    validate_alignment_sequences,
)
from .source_review import (
    EXPECTED_PROJECT15_COUNTS,
    ReferenceInput,
    load_reference_manifest,
    reproduce_project15,
    sha256,
    tool_versions,
    verify_reference_inputs,
    write_csv,
)


PROJECTION_FIELDS = [
    "source_accession",
    "source_label",
    "target_accession",
    "target_label",
    "source_variant_id",
    "source_variant_key",
    "source_contig",
    "source_position",
    "source_reference",
    "source_alternate",
    "source_variant_type",
    "status",
    "reason",
    "alignment_id",
    "strand",
    "target_contig",
    "target_position",
    "target_reference",
    "target_alternate",
    "target_variant_key",
    "reciprocal",
    "normalized",
]


@dataclass(frozen=True)
class Direction:
    source: ReferenceInput
    target: ReferenceInput

    @property
    def name(self) -> str:
        return f"{self.source.accession}_to_{self.target.accession}"


class UnionFind:
    def __init__(self, nodes: Iterable[tuple[str, str]]) -> None:
        self.parent = {node: node for node in nodes}
        self.rank = {node: 0 for node in self.parent}

    def find(self, node: tuple[str, str]) -> tuple[str, str]:
        parent = self.parent[node]
        if parent != node:
            self.parent[node] = self.find(parent)
        return self.parent[node]

    def union(self, left: tuple[str, str], right: tuple[str, str]) -> None:
        root_left, root_right = self.find(left), self.find(right)
        if root_left == root_right:
            return
        if self.rank[root_left] < self.rank[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if self.rank[root_left] == self.rank[root_right]:
            self.rank[root_left] += 1


def _run(command: list[str]) -> subprocess.CompletedProcess[str]:
    try:
        result = subprocess.run(command, capture_output=True, text=True)
    except OSError as error:
        raise LiftoverError(f"command could not start ({command[0]}): {error}") from error
    if result.returncode != 0:
        raise LiftoverError(f"command failed ({command[0]}): {result.stderr.strip()}")
    return result


def _projection_row(
    source: ReferenceInput,
    target: ReferenceInput,
    result: ProjectionResult,
) -> dict[str, object]:
    return {
        "source_accession": source.accession,
        "source_label": source.label,
        "target_accession": target.accession,
        "target_label": target.label,
        "source_variant_id": result.source.identifier,
        "source_variant_key": result.source.key,
        "source_contig": result.source.contig,
        "source_position": result.source.position,
        "source_reference": result.source.reference,
        "source_alternate": result.source.alternate,
        "source_variant_type": result.source.variant_type,
        "status": result.status,
        "reason": result.reason,
        "alignment_id": result.alignment_id,
        "strand": result.strand,
        "target_contig": result.target_contig,
        "target_position": "" if result.target_position is None else result.target_position,
        "target_reference": result.target_reference,
        "target_alternate": result.target_alternate,
        "target_variant_key": result.target_key,
        "reciprocal": str(result.reciprocal).lower(),
        "normalized": str(result.normalized).lower(),
    }


def _source_hashes(references: Sequence[ReferenceInput]) -> dict[str, str]:
    values: dict[str, str] = {}
    for reference in references:
        paths = {
            "reference_fasta": reference.reference_fasta,
            "accepted_vcf": reference.accepted_vcf,
            "evaluated_vcf": reference.evaluated_vcf,
            "rejected_vcf": reference.rejected_vcf,
            "callable_bed": reference.callable_bed,
            "accepted_vcf_index": Path(f"{reference.accepted_vcf}.csi"),
            "evaluated_vcf_index": Path(f"{reference.evaluated_vcf}.csi"),
            "rejected_vcf_index": Path(f"{reference.rejected_vcf}.csi"),
        }
        for name, path in paths.items():
            values[f"{reference.accession}:{name}"] = sha256(path)
    return values


def _write_vcf(
    path: Path,
    results: Sequence[ProjectionResult],
    sequences: Mapping[str, str],
    statuses: set[str],
) -> tuple[int, int]:
    grouped: dict[str, list[ProjectionResult]] = defaultdict(list)
    for result in results:
        if result.status in statuses and result.target_key:
            grouped[result.target_key].append(result)
    contig_order = {contig: index for index, contig in enumerate(sequences)}
    representatives = [items[0] for items in grouped.values()]
    representatives.sort(
        key=lambda item: (
            contig_order[item.target_contig],
            int(item.target_position),
            item.target_reference,
            item.target_alternate,
        )
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("##fileformat=VCFv4.3\n")
        handle.write('##INFO=<ID=LIFTOVER_STATUS,Number=1,Type=String,Description="Projection evidence status">\n')
        handle.write('##INFO=<ID=SOURCE_COUNT,Number=1,Type=Integer,Description="Source alleles yielding this normalized target key">\n')
        for contig, sequence in sequences.items():
            handle.write(f"##contig=<ID={contig},length={len(sequence)}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for index, result in enumerate(representatives, start=1):
            handle.write(
                f"{result.target_contig}\t{result.target_position}\tLIFT{index}\t"
                f"{result.target_reference}\t{result.target_alternate}\t.\tPASS\t"
                f"LIFTOVER_STATUS={result.status};SOURCE_COUNT={len(grouped[result.target_key])}\n"
            )
    duplicate_attempts = sum(len(items) - 1 for items in grouped.values())
    return len(representatives), duplicate_attempts


def _compress_and_validate_vcf(
    raw_vcf: Path,
    target_fasta: Path,
    compressed_vcf: Path,
    log_path: Path,
) -> int:
    result = _run(
        [
            "bcftools",
            "norm",
            "--no-version",
            "-c",
            "e",
            "-d",
            "exact",
            "-f",
            str(target_fasta),
            "-Oz",
            "-o",
            str(compressed_vcf),
            str(raw_vcf),
        ]
    )
    log_path.write_text(result.stderr, encoding="utf-8")
    _run(["bcftools", "index", "-f", "--csi", str(compressed_vcf)])
    count = int(_run(["bcftools", "index", "-n", str(compressed_vcf)]).stdout.strip() or 0)
    raw_vcf.unlink()
    return count


def _write_direction_vcfs(
    directory: Path,
    target: ReferenceInput,
    target_fasta: Path,
    target_sequences: Mapping[str, str],
    results: Sequence[ProjectionResult],
) -> dict[str, int]:
    directory.mkdir(parents=True)
    projected_raw = directory / "projected.normalized.vcf"
    exact_raw = directory / "exact_equivalence.vcf"
    projected_unique, projected_duplicates = _write_vcf(
        projected_raw,
        results,
        target_sequences,
        {"accepted_match", "rejected_match", "callable_no_candidate", "target_not_callable"},
    )
    exact_unique, exact_duplicates = _write_vcf(
        exact_raw,
        results,
        target_sequences,
        {"accepted_match"},
    )
    projected_count = _compress_and_validate_vcf(
        projected_raw,
        target_fasta,
        directory / "projected.normalized.vcf.gz",
        directory / "projected.bcftools.log.txt",
    )
    exact_count = _compress_and_validate_vcf(
        exact_raw,
        target_fasta,
        directory / "exact_equivalence.vcf.gz",
        directory / "exact_equivalence.bcftools.log.txt",
    )
    if projected_count != projected_unique or exact_count != exact_unique:
        raise LiftoverError(f"VCF normalization changed unique record totals for {target.accession}")
    return {
        "projected_vcf_records": projected_count,
        "projected_duplicate_attempts": projected_duplicates,
        "exact_vcf_records": exact_count,
        "exact_duplicate_attempts": exact_duplicates,
    }


def _write_clusters(
    path: Path,
    accepted_variants: Mapping[str, Sequence[Variant]],
    accepted_edges: Sequence[tuple[tuple[str, str], tuple[str, str]]],
    baseline_accession: str,
) -> dict[str, int]:
    nodes = [(accession, variant.key) for accession, variants in accepted_variants.items() for variant in variants]
    union = UnionFind(nodes)
    for left, right in accepted_edges:
        if left not in union.parent or right not in union.parent:
            raise LiftoverError("accepted projection edge refers to a non-accepted variant")
        union.union(left, right)
    components: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    for node in nodes:
        components[union.find(node)].append(node)
    ordered_components = sorted((sorted(items) for items in components.values()), key=lambda items: items[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "cluster_id",
            "accession",
            "variant_key",
            "cluster_size",
            "reference_count",
            "includes_baseline",
            "exact_edge_supported",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for index, members in enumerate(ordered_components, start=1):
            accessions = {item[0] for item in members}
            for accession, variant_key in members:
                writer.writerow(
                    {
                        "cluster_id": f"EQ{index:07d}",
                        "accession": accession,
                        "variant_key": variant_key,
                        "cluster_size": len(members),
                        "reference_count": len(accessions),
                        "includes_baseline": str(baseline_accession in accessions).lower(),
                        "exact_edge_supported": str(len(members) > 1).lower(),
                    }
                )
    return {
        "clusters": len(ordered_components),
        "cluster_members": len(nodes),
        "multi_reference_clusters": sum(len({item[0] for item in members}) > 1 for members in ordered_components),
        "variants_in_multi_reference_clusters": sum(
            len(members) for members in ordered_components if len({item[0] for item in members}) > 1
        ),
    }


def _summary_rows(
    direction_counts: Mapping[tuple[str, str], Counter[str]],
    labels: Mapping[str, str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (source, target), counts in direction_counts.items():
        total = sum(counts.values())
        for status in sorted(FINAL_STATUSES):
            count = counts[status]
            rows.append(
                {
                    "source_accession": source,
                    "source_label": labels[source],
                    "target_accession": target,
                    "target_label": labels[target],
                    "status": status,
                    "count": count,
                    "percent": count * 100 / total if total else 0,
                }
            )
    return rows


def _write_dashboard(
    path: Path,
    coverage_path: Path,
    baseline_accession: str,
    alternatives: Sequence[ReferenceInput],
    direction_counts: Mapping[tuple[str, str], Counter[str]],
    explanations: Sequence[Mapping[str, object]],
    reviews: Sequence[Mapping[str, object]],
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    with coverage_path.open(newline="", encoding="utf-8") as handle:
        coverage = {row["alternative_accession"]: row for row in csv.DictReader(handle)}
    labels = [item.label.replace(" / ", "\n") for item in alternatives]
    accessions = [item.accession for item in alternatives]
    figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)

    x = np.arange(len(alternatives))
    baseline_cov = [float(coverage[item]["baseline_reciprocal_percent"]) for item in accessions]
    alternative_cov = [float(coverage[item]["alternative_reciprocal_percent"]) for item in accessions]
    width = 0.38
    axes[0, 0].bar(x - width / 2, baseline_cov, width, label="MG1655")
    axes[0, 0].bar(x + width / 2, alternative_cov, width, label="Alternative")
    axes[0, 0].axhline(70, color="#555555", linestyle="--", linewidth=1)
    axes[0, 0].set_title("Reciprocal genome coverage")
    axes[0, 0].set_ylabel("Genome covered (%)")
    axes[0, 0].set_xticks(x, labels, fontsize=8)
    axes[0, 0].legend(frameon=False, fontsize=8)

    projection_groups = [
        ("accepted_match", "Accepted exact", "#2b8cbe"),
        ("rejected_match", "Rejected exact", "#a6bddb"),
        ("target_reference_match", "Target reference", "#7fcdbb"),
        ("callable_no_candidate", "Callable absent", "#fdae6b"),
        ("target_not_callable", "Not callable", "#e6550d"),
        ("unresolved", "Unresolved", "#969696"),
    ]
    bottoms = np.zeros(len(alternatives))
    for status, label, color in projection_groups:
        values = []
        for accession in accessions:
            counts = direction_counts[(baseline_accession, accession)]
            if status == "unresolved":
                value = sum(
                    counts[item]
                    for item in (
                        "ambiguous_mapping",
                        "nonreciprocal_mapping",
                        "nonalignable",
                        "complex_alignment_context",
                    )
                )
            else:
                value = counts[status]
            values.append(value)
        axes[0, 1].bar(x, values, bottom=bottoms, label=label, color=color)
        bottoms += np.array(values)
    axes[0, 1].set_title("MG1655 variant projection outcomes")
    axes[0, 1].set_ylabel("Baseline variants")
    axes[0, 1].set_xticks(x, labels, fontsize=8)
    axes[0, 1].legend(frameon=False, fontsize=7, ncol=2)

    explanation_counts = Counter(
        (str(row["alternative_accession"]), str(row["burden_explanation"]))
        for row in explanations
    )
    bottoms = np.zeros(len(alternatives))
    colors = ["#238b45", "#74c476", "#9ecae1", "#fdae6b", "#756bb1", "#969696"]
    for category, color in zip(BURDEN_CATEGORIES, colors):
        values = [explanation_counts[(accession, category)] for accession in accessions]
        axes[1, 0].bar(x, values, bottom=bottoms, label=category.replace("_", " "), color=color)
        bottoms += np.array(values)
    axes[1, 0].set_title("Comparable-gene burden explanations")
    axes[1, 0].set_ylabel("Gene pairs")
    axes[1, 0].set_xticks(x, labels, fontsize=8)
    axes[1, 0].legend(frameon=False, fontsize=6.5, ncol=2)

    review_counts = Counter(
        (str(row["alternative_accession"]), str(row["projection_status"])) for row in reviews
    )
    review_statuses = [
        "accepted_match",
        "rejected_match",
        "target_reference_match",
        "callable_no_candidate",
        "target_not_callable",
        "ambiguous_mapping",
        "nonreciprocal_mapping",
        "nonalignable",
        "complex_alignment_context",
    ]
    bottoms = np.zeros(len(alternatives))
    palette = plt.get_cmap("tab10")
    for index, status in enumerate(review_statuses):
        values = [review_counts[(accession, status)] for accession in accessions]
        if not any(values):
            continue
        axes[1, 1].bar(x, values, bottom=bottoms, label=status.replace("_", " "), color=palette(index))
        bottoms += np.array(values)
    axes[1, 1].set_title("Review-candidate liftover outcomes")
    axes[1, 1].set_ylabel("Candidate contexts")
    axes[1, 1].set_xticks(x, labels, fontsize=8)
    axes[1, 1].legend(frameon=False, fontsize=6.5, ncol=2)

    figure.suptitle("Cross-reference variant equivalence and Project 15 explanations", fontsize=15)
    figure.savefig(path, dpi=200)
    plt.close(figure)


def _hash_outputs(root: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.name == "run_manifest.json":
            continue
        output[path.relative_to(root).as_posix()] = sha256(path)
    return output


def authentic_analysis(args) -> dict[str, object]:
    output_dir = args.output_dir.resolve()
    work_dir = args.work_dir.resolve()
    if output_dir.exists() or work_dir.exists():
        raise LiftoverError("output and work directories must not already exist")
    config = ProjectionConfig(
        minimum_alignment_length=args.minimum_alignment_length,
        minimum_alignment_identity=args.minimum_alignment_identity,
        minimum_alignment_mapq=args.minimum_alignment_mapq,
        ambiguity_score_fraction=args.ambiguity_score_fraction,
    )
    config.validate()
    references = load_reference_manifest(args.reference_manifest.resolve(), args.baseline_accession)
    versions = tool_versions()
    _validation, inventories, source_hashes_before = verify_reference_inputs(references, args.sample_id)
    _reproduction, project15_counts, project15_hashes = reproduce_project15(
        args.project15_comparison_dir.resolve(),
        args.project15_orthology_dir.resolve(),
        args.project15_annotation_dir.resolve(),
    )
    source_hashes_before = _source_hashes(references)
    baseline = next(item for item in references if item.role == "baseline")
    alternatives = [item for item in references if item.role == "alternative"]
    labels = {item.accession: item.label for item in references}
    expected_attempts = sum(
        int(inventories[baseline.accession]["accepted_variants"])
        + int(inventories[item.accession]["accepted_variants"])
        for item in alternatives
    )
    expected_baseline_rows = int(inventories[baseline.accession]["accepted_variants"]) * len(alternatives)
    if expected_attempts != 481_473 or expected_baseline_rows != 154_865:
        raise LiftoverError("inherited variant counts do not yield the required projection totals")

    run_id = uuid.uuid4().hex
    output_stage = output_dir.parent / f".{output_dir.name}.staging-{run_id}"
    work_stage = work_dir.parent / f".{work_dir.name}.staging-{run_id}"
    output_stage.mkdir(parents=True)
    work_stage.mkdir(parents=True)
    for name in ("tables", "summaries", "per_direction", "dashboard", "reports"):
        (output_stage / name).mkdir()
    reference_cache = work_stage / "reference_cache"
    normalization_work = work_stage / "normalization"
    reference_cache.mkdir()
    normalization_work.mkdir()
    try:
        copied_fastas: dict[str, Path] = {}
        sequences: dict[str, dict[str, str]] = {}
        accepted_variants: dict[str, list[Variant]] = {}
        accepted_keys: dict[str, set[str]] = {}
        rejected_keys: dict[str, set[str]] = {}
        callable_intervals: dict[str, dict[str, list[tuple[int, int]]]] = {}
        for reference in references:
            copied = reference_cache / f"{reference.accession}.fasta"
            shutil.copyfile(reference.reference_fasta, copied)
            if sha256(copied) != reference.expected_hashes["reference_fasta"]:
                raise LiftoverError(f"reference cache checksum mismatch for {reference.accession}")
            _run(["samtools", "faidx", str(copied)])
            copied_fastas[reference.accession] = copied
            sequences[reference.accession] = read_fasta(copied)
            accepted_variants[reference.accession] = read_vcf(reference.accepted_vcf)
            accepted_keys[reference.accession] = {
                item.key for item in accepted_variants[reference.accession]
            }
            rejected_keys[reference.accession] = {
                item.key for item in read_vcf(reference.rejected_vcf)
            }
            callable_intervals[reference.accession] = read_bed(reference.callable_bed)

        projection_path = output_stage / "tables" / "variant_projections.csv.gz"
        baseline_path = output_stage / "tables" / "baseline_variant_comparison.csv"
        direction_counts: dict[tuple[str, str], Counter[str]] = {}
        accepted_edges: list[tuple[tuple[str, str], tuple[str, str]]] = []
        baseline_evidence: list[ProjectionEvidence] = []
        vcf_counts: dict[str, dict[str, int]] = {}
        with gzip.open(projection_path, "wt", newline="", encoding="utf-8") as projection_handle, baseline_path.open(
            "w", newline="", encoding="utf-8"
        ) as baseline_handle:
            projection_writer = csv.DictWriter(
                projection_handle, fieldnames=PROJECTION_FIELDS, lineterminator="\n"
            )
            baseline_writer = csv.DictWriter(
                baseline_handle, fieldnames=PROJECTION_FIELDS, lineterminator="\n"
            )
            projection_writer.writeheader()
            baseline_writer.writeheader()

            for alternative in alternatives:
                paf_baseline_to_alt = (
                    args.alignment_dir.resolve()
                    / f"{baseline.accession}_to_{alternative.accession}.paf"
                )
                paf_alt_to_baseline = (
                    args.alignment_dir.resolve()
                    / f"{alternative.accession}_to_{baseline.accession}.paf"
                )
                forward = parse_paf(paf_baseline_to_alt)
                reverse = parse_paf(paf_alt_to_baseline)
                for alignment in forward:
                    validate_alignment_sequences(
                        alignment,
                        sequences[baseline.accession],
                        sequences[alternative.accession],
                    )
                for alignment in reverse:
                    validate_alignment_sequences(
                        alignment,
                        sequences[alternative.accession],
                        sequences[baseline.accession],
                    )

                pair_work = normalization_work / alternative.accession
                baseline_results = project_variants_loaded(
                    accepted_variants[baseline.accession],
                    sequences[baseline.accession],
                    sequences[alternative.accession],
                    forward,
                    reverse,
                    accepted_keys[alternative.accession],
                    rejected_keys[alternative.accession],
                    callable_intervals[alternative.accession],
                    config,
                    copied_fastas[baseline.accession],
                    copied_fastas[alternative.accession],
                    pair_work / "baseline_to_alternative",
                )
                direction_counts[(baseline.accession, alternative.accession)] = Counter(
                    item.status for item in baseline_results
                )
                for result in baseline_results:
                    row = _projection_row(baseline, alternative, result)
                    projection_writer.writerow(row)
                    baseline_writer.writerow(row)
                    baseline_evidence.append(
                        ProjectionEvidence(
                            alternative_accession=alternative.accession,
                            source_variant_key=result.source.key,
                            status=result.status,
                            target_variant_key=result.target_key,
                            reason=result.reason,
                        )
                    )
                    if result.status == "accepted_match":
                        accepted_edges.append(
                            (
                                (baseline.accession, result.source.key),
                                (alternative.accession, result.target_key),
                            )
                        )
                baseline_direction = f"{baseline.accession}_to_{alternative.accession}"
                vcf_counts[baseline_direction] = _write_direction_vcfs(
                    output_stage / "per_direction" / baseline_direction,
                    alternative,
                    copied_fastas[alternative.accession],
                    sequences[alternative.accession],
                    baseline_results,
                )

                alternative_results = project_variants_loaded(
                    accepted_variants[alternative.accession],
                    sequences[alternative.accession],
                    sequences[baseline.accession],
                    reverse,
                    forward,
                    accepted_keys[baseline.accession],
                    rejected_keys[baseline.accession],
                    callable_intervals[baseline.accession],
                    config,
                    copied_fastas[alternative.accession],
                    copied_fastas[baseline.accession],
                    pair_work / "alternative_to_baseline",
                )
                direction_counts[(alternative.accession, baseline.accession)] = Counter(
                    item.status for item in alternative_results
                )
                for result in alternative_results:
                    projection_writer.writerow(_projection_row(alternative, baseline, result))
                    if result.status == "accepted_match":
                        accepted_edges.append(
                            (
                                (alternative.accession, result.source.key),
                                (baseline.accession, result.target_key),
                            )
                        )
                alternative_direction = f"{alternative.accession}_to_{baseline.accession}"
                vcf_counts[alternative_direction] = _write_direction_vcfs(
                    output_stage / "per_direction" / alternative_direction,
                    baseline,
                    copied_fastas[baseline.accession],
                    sequences[baseline.accession],
                    alternative_results,
                )
                shutil.rmtree(pair_work, ignore_errors=True)
                print(
                    f"completed {baseline.accession} <-> {alternative.accession}: "
                    f"{len(baseline_results) + len(alternative_results)} attempts",
                    flush=True,
                )

        observed_attempts = sum(sum(counts.values()) for counts in direction_counts.values())
        if observed_attempts != expected_attempts or len(baseline_evidence) != expected_baseline_rows:
            raise LiftoverError("directional projection totals do not match the acceptance requirements")

        cluster_counts = _write_clusters(
            output_stage / "tables" / "variant_equivalence_clusters.csv",
            accepted_variants,
            accepted_edges,
            baseline.accession,
        )

        comparison_dir = args.project15_comparison_dir.resolve()
        annotation_dir = args.project15_annotation_dir.resolve()
        gene_pairs, comparison_rows = load_project15_gene_pairs(
            comparison_dir / "tables" / "baseline_ortholog_comparison.csv"
        )
        gene_variants, effect_rows, unlinked_effects = load_project15_gene_variants(
            annotation_dir / "tables" / "annotation_effects.csv",
            annotation_dir / "tables" / "gene_catalog.csv",
        )
        burden_rows = validate_project15_burdens(
            comparison_dir / "tables" / "gene_consequence_burden.csv",
            gene_variants,
        )
        review_contexts = load_project15_review_contexts(
            comparison_dir / "tables" / "review_candidate_orthology.csv"
        )
        equivalence_rows, explanation_rows = integrate_gene_burdens(
            gene_pairs,
            gene_variants,
            baseline_evidence,
        )
        review_rows = connect_review_candidates(
            review_contexts,
            baseline_evidence,
            equivalence_rows,
            explanation_rows,
        )
        if len(explanation_rows) != 6_356 or len(review_rows) != 400:
            raise LiftoverError("Project 15 integration output totals are incorrect")
        write_rows(output_stage / "tables" / "ortholog_variant_equivalence.csv", equivalence_rows)
        write_rows(output_stage / "tables" / "project15_burden_explanation.csv", explanation_rows)
        write_rows(output_stage / "tables" / "review_candidate_liftover.csv", review_rows)

        projection_summary = _summary_rows(direction_counts, labels)
        write_csv(output_stage / "summaries" / "projection_status_summary.csv", projection_summary)
        explanation_summary = Counter(
            (str(row["alternative_accession"]), str(row["burden_explanation"]))
            for row in explanation_rows
        )
        write_csv(
            output_stage / "summaries" / "burden_explanation_summary.csv",
            [
                {
                    "alternative_accession": accession,
                    "alternative_label": labels[accession],
                    "burden_explanation": category,
                    "gene_pair_count": explanation_summary[(accession, category)],
                }
                for accession in [item.accession for item in alternatives]
                for category in BURDEN_CATEGORIES
            ],
        )
        effect_summary = Counter(str(row["effect_relation"]) for row in equivalence_rows)
        write_csv(
            output_stage / "summaries" / "effect_relation_summary.csv",
            [{"effect_relation": key, "variant_gene_rows": value} for key, value in sorted(effect_summary.items())],
        )
        review_summary = Counter(
            (str(row["alternative_accession"]), str(row["projection_status"])) for row in review_rows
        )
        write_csv(
            output_stage / "summaries" / "review_liftover_summary.csv",
            [
                {
                    "alternative_accession": accession,
                    "alternative_label": labels[accession],
                    "projection_status": status,
                    "review_candidate_count": review_summary[(accession, status)],
                }
                for accession in [item.accession for item in alternatives]
                for status in sorted(FINAL_STATUSES)
            ],
        )
        write_csv(
            output_stage / "summaries" / "reference_vcf_summary.csv",
            [
                {
                    "direction": f"{source.accession}_to_{target.accession}",
                    "source_accession": source.accession,
                    "source_label": source.label,
                    "target_accession": target.accession,
                    "target_label": target.label,
                    **vcf_counts[f"{source.accession}_to_{target.accession}"],
                }
                for item in alternatives
                for source, target in ((baseline, item), (item, baseline))
            ],
        )
        write_csv(
            output_stage / "summaries" / "equivalence_cluster_summary.csv",
            [{"metric": key, "count": value} for key, value in cluster_counts.items()],
        )

        _write_dashboard(
            output_stage / "dashboard" / "variant_liftover_dashboard.png",
            args.connection_review_dir.resolve() / "tables" / "reciprocal_coverage.csv",
            baseline.accession,
            alternatives,
            direction_counts,
            explanation_rows,
            review_rows,
        )

        source_hashes_after = _source_hashes(references)
        if source_hashes_after != source_hashes_before:
            raise LiftoverError("one or more inherited Project 13 source files changed")
        project15_after = {
            "annotation_effects": sha256(annotation_dir / "tables" / "annotation_effects.csv"),
            "baseline_ortholog_comparison": sha256(
                comparison_dir / "tables" / "baseline_ortholog_comparison.csv"
            ),
            "review_candidate_orthology": sha256(
                comparison_dir / "tables" / "review_candidate_orthology.csv"
            ),
            "gene_consequence_burden": sha256(
                comparison_dir / "tables" / "gene_consequence_burden.csv"
            ),
        }
        for name, value in project15_after.items():
            if value != project15_hashes[name]:
                raise LiftoverError(f"Project 15 input checksum changed: {name}")

        (output_stage / "portable_commands.txt").write_text(
            "minimap2 -cx asm20 --cs=long --eqx --secondary=yes -N 5 -t 1 <TARGET_FASTA> <SOURCE_FASTA> > <PAF>\n"
            "bcftools norm --no-version -f <TARGET_FASTA> -m -any -Ov -o <NORMALIZED_VCF> <PROVISIONAL_VCF>\n"
            "bcftools norm --no-version -c e -d exact -f <TARGET_FASTA> -Oz -o <PROJECTED_VCF_GZ> <PROJECTED_VCF>\n"
            "bcftools index -f --csi <PROJECTED_VCF_GZ>\n",
            encoding="utf-8",
        )
        manifest = {
            "workflow": "variant-liftover-explorer",
            "workflow_version": "0.1.0",
            "milestone": "authentic_analysis",
            "status": "pass",
            "baseline_accession": baseline.accession,
            "sample_id": args.sample_id,
            "parameters": {
                "alignment_preset": args.alignment_preset,
                "minimum_alignment_length": config.minimum_alignment_length,
                "minimum_alignment_identity": config.minimum_alignment_identity,
                "minimum_alignment_mapq": config.minimum_alignment_mapq,
                "ambiguity_score_fraction": config.ambiguity_score_fraction,
                "minimum_comparable_coverage": args.minimum_comparable_coverage,
                "threads": args.threads,
            },
            "counts": {
                "directional_projection_attempts": observed_attempts,
                "baseline_variant_comparison_rows": len(baseline_evidence),
                "accepted_projection_edges": len(accepted_edges),
                "project15_variants": EXPECTED_PROJECT15_COUNTS["accepted_variants"],
                "project15_annotation_effects": effect_rows,
                "project15_variant_gene_assignments": len(gene_variants),
                "project15_unlinked_effects": unlinked_effects,
                "project15_gene_burden_rows": burden_rows,
                "comparable_gene_pairs": len(explanation_rows),
                "ortholog_variant_equivalence_rows": len(equivalence_rows),
                "review_candidate_liftover_rows": len(review_rows),
                **cluster_counts,
            },
            "direction_status_counts": {
                f"{source}_to_{target}": dict(counts)
                for (source, target), counts in direction_counts.items()
            },
            "burden_explanation_counts": dict(Counter(row["burden_explanation"] for row in explanation_rows)),
            "review_projection_counts": dict(Counter(row["projection_status"] for row in review_rows)),
            "vcf_counts": vcf_counts,
            "tools": versions,
            "source_hashes": source_hashes_before,
            "project15_hashes": project15_hashes,
            "source_inputs_unchanged": True,
            "interpretation_boundary": (
                "Exact allele equivalence is reference- and alignment-dependent; unresolved results are not biological absence."
            ),
        }
        manifest["output_hashes"] = _hash_outputs(output_stage)
        (output_stage / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        shutil.rmtree(work_stage)
        os.replace(output_stage, output_dir)
        return manifest
    except Exception:
        shutil.rmtree(output_stage, ignore_errors=True)
        shutil.rmtree(work_stage, ignore_errors=True)
        raise
