"""Milestone 1 source validation and reciprocal whole-genome alignment review."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


EXPECTED_PROJECT15_COUNTS = {
    "references": 6,
    "accepted_variants": 357_581,
    "annotation_effects": 357_959,
    "baseline_alternative_gene_comparisons": 23_255,
    "comparable_gene_pairs": 6_356,
    "project14_review_candidates": 80,
    "review_candidate_reference_contexts": 400,
}

REQUIRED_MANIFEST_COLUMNS = (
    "accession",
    "label",
    "role",
    "reference_fasta",
    "accepted_vcf",
    "evaluated_vcf",
    "rejected_vcf",
    "callable_bed",
    "reference_sha256",
    "accepted_vcf_sha256",
    "evaluated_vcf_sha256",
    "rejected_vcf_sha256",
    "callable_bed_sha256",
)

CS_TOKEN = re.compile(
    r"=[A-Za-z]+|:[0-9]+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+|~[A-Za-z]{2}[0-9]+[A-Za-z]{2}"
)
CIGAR_TOKEN = re.compile(r"[0-9]+[MIDNSHP=X]")


class ReviewError(RuntimeError):
    """Raised when a Milestone 1 invariant is not satisfied."""


@dataclass(frozen=True)
class ReferenceInput:
    accession: str
    label: str
    role: str
    reference_fasta: Path
    accepted_vcf: Path
    evaluated_vcf: Path
    rejected_vcf: Path
    callable_bed: Path
    expected_hashes: dict[str, str]


@dataclass
class AlignmentBlock:
    block_id: str
    direction: str
    source_accession: str
    target_accession: str
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    alignment_bases: int
    mapq: int
    alignment_type: str
    chain_score: int | None
    secondary_score: int | None
    identity: float
    qualifying: bool
    reciprocal: bool = False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve(base: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def load_reference_manifest(path: Path, baseline_accession: str) -> list[ReferenceInput]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != list(REQUIRED_MANIFEST_COLUMNS):
            raise ReviewError("reference manifest columns do not match the required schema")
        rows = list(reader)
    if len(rows) != 6:
        raise ReviewError(f"expected six references, found {len(rows)}")
    accessions = [row["accession"] for row in rows]
    if len(set(accessions)) != len(accessions):
        raise ReviewError("reference manifest contains duplicate accessions")
    if accessions.count(baseline_accession) != 1:
        raise ReviewError("baseline accession must occur exactly once")
    if sum(row["role"] == "baseline" for row in rows) != 1:
        raise ReviewError("reference manifest must contain exactly one baseline role")
    base = path.resolve().parent
    references: list[ReferenceInput] = []
    for row in rows:
        if row["role"] not in {"baseline", "alternative"}:
            raise ReviewError(f"invalid role for {row['accession']}: {row['role']}")
        if (row["accession"] == baseline_accession) != (row["role"] == "baseline"):
            raise ReviewError("baseline accession and role disagree")
        references.append(
            ReferenceInput(
                accession=row["accession"],
                label=row["label"],
                role=row["role"],
                reference_fasta=_resolve(base, row["reference_fasta"]),
                accepted_vcf=_resolve(base, row["accepted_vcf"]),
                evaluated_vcf=_resolve(base, row["evaluated_vcf"]),
                rejected_vcf=_resolve(base, row["rejected_vcf"]),
                callable_bed=_resolve(base, row["callable_bed"]),
                expected_hashes={
                    "reference_fasta": row["reference_sha256"],
                    "accepted_vcf": row["accepted_vcf_sha256"],
                    "evaluated_vcf": row["evaluated_vcf_sha256"],
                    "rejected_vcf": row["rejected_vcf_sha256"],
                    "callable_bed": row["callable_bed_sha256"],
                },
            )
        )
    return references


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    with path.open(encoding="ascii") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if not name or name in lengths:
                    raise ReviewError(f"invalid or duplicate FASTA identifier at {path}:{line_number}")
                lengths[name] = 0
            elif name is None:
                raise ReviewError(f"FASTA sequence precedes header at {path}:{line_number}")
            else:
                sequence = line.upper()
                if re.search(r"[^ACGTN]", sequence):
                    raise ReviewError(f"invalid FASTA base at {path}:{line_number}")
                lengths[name] += len(sequence)
    if not lengths or any(length == 0 for length in lengths.values()):
        raise ReviewError(f"empty FASTA record in {path}")
    return lengths


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def read_vcf_keys(path: Path, sample_id: str) -> set[tuple[str, int, str, str]]:
    keys: set[tuple[str, int, str, str]] = set()
    header_seen = False
    with _open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            if raw.startswith("##"):
                continue
            if raw.startswith("#CHROM"):
                columns = raw.rstrip("\n").split("\t")
                samples = columns[9:]
                if samples != [sample_id]:
                    raise ReviewError(f"{path} contains samples {samples}, expected [{sample_id}]")
                header_seen = True
                continue
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise ReviewError(f"malformed VCF row at {path}:{line_number}")
            chrom, pos_text, _identifier, ref, alt = fields[:5]
            if "," in alt or alt.startswith("<") or "[" in alt or "]" in alt:
                raise ReviewError(f"non-biallelic or symbolic allele at {path}:{line_number}")
            key = (chrom, int(pos_text), ref.upper(), alt.upper())
            if key in keys:
                raise ReviewError(f"duplicate VCF key at {path}:{line_number}: {key}")
            keys.add(key)
    if not header_seen:
        raise ReviewError(f"missing #CHROM header in {path}")
    return keys


def csv_rows(path: Path) -> Iterator[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def verify_reference_inputs(
    references: list[ReferenceInput], sample_id: str
) -> tuple[list[dict[str, object]], dict[str, dict[str, object]], dict[str, str]]:
    validation_rows: list[dict[str, object]] = []
    inventories: dict[str, dict[str, object]] = {}
    source_hashes: dict[str, str] = {}
    for reference in references:
        files = {
            "reference_fasta": reference.reference_fasta,
            "accepted_vcf": reference.accepted_vcf,
            "evaluated_vcf": reference.evaluated_vcf,
            "rejected_vcf": reference.rejected_vcf,
            "callable_bed": reference.callable_bed,
        }
        for kind, file_path in files.items():
            if not file_path.is_file():
                raise ReviewError(f"missing {kind} for {reference.accession}: {file_path}")
            observed = sha256(file_path)
            expected = reference.expected_hashes[kind]
            if observed != expected:
                raise ReviewError(f"checksum mismatch for {reference.accession} {kind}")
            source_hashes[f"{reference.accession}:{kind}"] = observed
            validation_rows.append(
                {
                    "accession": reference.accession,
                    "role": reference.role,
                    "input_kind": kind,
                    "sha256": observed,
                    "status": "pass",
                }
            )
        index_paths = {
            "accepted_vcf_index": Path(f"{reference.accepted_vcf}.csi"),
            "evaluated_vcf_index": Path(f"{reference.evaluated_vcf}.csi"),
            "rejected_vcf_index": Path(f"{reference.rejected_vcf}.csi"),
        }
        for kind, index_path in index_paths.items():
            if not index_path.is_file():
                raise ReviewError(f"missing {kind} for {reference.accession}")
            observed = sha256(index_path)
            source_hashes[f"{reference.accession}:{kind}"] = observed
            validation_rows.append(
                {
                    "accession": reference.accession,
                    "role": reference.role,
                    "input_kind": kind,
                    "sha256": observed,
                    "status": "pass",
                }
            )
        accepted = read_vcf_keys(reference.accepted_vcf, sample_id)
        evaluated = read_vcf_keys(reference.evaluated_vcf, sample_id)
        rejected = read_vcf_keys(reference.rejected_vcf, sample_id)
        if accepted & rejected:
            raise ReviewError(f"accepted and rejected VCFs overlap for {reference.accession}")
        if accepted | rejected != evaluated:
            raise ReviewError(f"evaluated VCF is not the accepted/rejected union for {reference.accession}")
        lengths = fasta_lengths(reference.reference_fasta)
        inventories[reference.accession] = {
            "label": reference.label,
            "role": reference.role,
            "contigs": len(lengths),
            "bases": sum(lengths.values()),
            "lengths": lengths,
            "accepted_variants": len(accepted),
            "evaluated_variants": len(evaluated),
            "rejected_variants": len(rejected),
        }
    return validation_rows, inventories, source_hashes


def reproduce_project15(
    comparison_dir: Path, orthology_dir: Path, annotation_dir: Path
) -> tuple[list[dict[str, object]], dict[str, int], dict[str, str]]:
    comparison_manifest_path = comparison_dir / "run_manifest.json"
    orthology_manifest_path = orthology_dir / "run_manifest.json"
    annotation_manifest_path = annotation_dir / "run_manifest.json"
    for path in (comparison_manifest_path, orthology_manifest_path, annotation_manifest_path):
        if not path.is_file():
            raise ReviewError(f"missing Project 15 manifest: {path}")
    comparison_manifest = json.loads(comparison_manifest_path.read_text(encoding="utf-8"))
    recorded = comparison_manifest.get("counts", {})
    for name, expected in EXPECTED_PROJECT15_COUNTS.items():
        if recorded.get(name) != expected:
            raise ReviewError(f"Project 15 manifest count {name} is not {expected}")

    effects_path = annotation_dir / "tables" / "annotation_effects.csv"
    comparisons_path = comparison_dir / "tables" / "baseline_ortholog_comparison.csv"
    reviews_path = comparison_dir / "tables" / "review_candidate_orthology.csv"
    burden_path = comparison_dir / "tables" / "gene_consequence_burden.csv"
    orthology_path = orthology_dir / "tables" / "orthology_relationships.csv"
    for path in (effects_path, comparisons_path, reviews_path, burden_path, orthology_path):
        if not path.is_file():
            raise ReviewError(f"missing Project 15 table: {path}")

    effect_count = sum(1 for _ in csv_rows(effects_path))
    comparisons = list(csv_rows(comparisons_path))
    reviews = list(csv_rows(reviews_path))
    comparable_count = sum(row["comparable"].lower() == "true" for row in comparisons)
    unique_review_candidates = {row["baseline_variant_key"] for row in reviews}
    observed = {
        "accepted_variants": int(recorded["accepted_variants"]),
        "annotation_effects": effect_count,
        "baseline_alternative_gene_comparisons": len(comparisons),
        "comparable_gene_pairs": comparable_count,
        "project14_review_candidates": len(unique_review_candidates),
        "review_candidate_reference_contexts": len(reviews),
    }
    for name, expected in EXPECTED_PROJECT15_COUNTS.items():
        if name == "references":
            continue
        if observed[name] != expected:
            raise ReviewError(f"Project 15 table count {name} is {observed[name]}, expected {expected}")

    rows = [
        {
            "metric": name,
            "expected": expected,
            "observed": EXPECTED_PROJECT15_COUNTS["references"] if name == "references" else observed[name],
            "status": "pass",
        }
        for name, expected in EXPECTED_PROJECT15_COUNTS.items()
    ]
    hashes = {
        "comparison_manifest": sha256(comparison_manifest_path),
        "orthology_manifest": sha256(orthology_manifest_path),
        "annotation_manifest": sha256(annotation_manifest_path),
        "annotation_effects": sha256(effects_path),
        "baseline_ortholog_comparison": sha256(comparisons_path),
        "review_candidate_orthology": sha256(reviews_path),
        "gene_consequence_burden": sha256(burden_path),
        "orthology_relationships": sha256(orthology_path),
    }
    return rows, observed, hashes


def tool_versions() -> dict[str, str]:
    def first_line(command: list[str]) -> str:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return (result.stdout or result.stderr).splitlines()[0].strip()

    import matplotlib

    versions = {
        "python": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "minimap2": first_line(["minimap2", "--version"]),
        "bcftools": first_line(["bcftools", "--version"]),
        "samtools": first_line(["samtools", "--version"]),
        "matplotlib": matplotlib.__version__,
    }
    expected_prefixes = {
        "python": "3.12.",
        "minimap2": "2.31",
        "bcftools": "bcftools 1.24",
        "samtools": "samtools 1.24",
        "matplotlib": "3.11.1",
    }
    for tool, prefix in expected_prefixes.items():
        if not versions[tool].startswith(prefix):
            raise ReviewError(f"unexpected {tool} version: {versions[tool]}")
    return versions


def _parse_tag(fields: list[str], tag: str) -> str | None:
    prefix = f"{tag}:"
    for field in fields:
        if field.startswith(prefix):
            return field.split(":", 2)[2]
    return None


def _validate_alignment_string(value: str, token_pattern: re.Pattern[str], name: str) -> None:
    tokens = token_pattern.findall(value)
    if not tokens or "".join(tokens) != value:
        raise ReviewError(f"invalid {name} string in PAF")


def parse_paf(
    path: Path,
    direction: str,
    source_accession: str,
    target_accession: str,
    minimum_length: int,
    minimum_identity: float,
    minimum_mapq: int,
) -> list[AlignmentBlock]:
    blocks: list[AlignmentBlock] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ReviewError(f"malformed PAF row at {path}:{line_number}")
            tags = fields[12:]
            cs = _parse_tag(tags, "cs")
            cigar = _parse_tag(tags, "cg")
            if cs is None or cigar is None:
                raise ReviewError(f"PAF row lacks cs or cg tag at {path}:{line_number}")
            _validate_alignment_string(cs, CS_TOKEN, "cs")
            _validate_alignment_string(cigar, CIGAR_TOKEN, "CIGAR")
            alignment_bases = int(fields[10])
            identity = int(fields[9]) / alignment_bases if alignment_bases else 0.0
            alignment_type = _parse_tag(tags, "tp") or ""
            chain_score_text = _parse_tag(tags, "s1")
            secondary_score_text = _parse_tag(tags, "s2")
            block = AlignmentBlock(
                block_id=f"{direction}:{line_number}",
                direction=direction,
                source_accession=source_accession,
                target_accession=target_accession,
                query_name=fields[0],
                query_length=int(fields[1]),
                query_start=int(fields[2]),
                query_end=int(fields[3]),
                strand=fields[4],
                target_name=fields[5],
                target_length=int(fields[6]),
                target_start=int(fields[7]),
                target_end=int(fields[8]),
                matches=int(fields[9]),
                alignment_bases=alignment_bases,
                mapq=int(fields[11]),
                alignment_type=alignment_type,
                chain_score=int(chain_score_text) if chain_score_text else None,
                secondary_score=int(secondary_score_text) if secondary_score_text else None,
                identity=identity,
                qualifying=(
                    alignment_type == "P"
                    and alignment_bases >= minimum_length
                    and identity >= minimum_identity
                    and int(fields[11]) >= minimum_mapq
                ),
            )
            blocks.append(block)
    if not blocks:
        raise ReviewError(f"alignment produced no PAF records: {path}")
    if not any(block.qualifying for block in blocks):
        raise ReviewError(f"alignment produced no qualifying primary blocks: {path}")
    return blocks


def _overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    return max(0, min(end_a, end_b) - max(start_a, start_b))


def _reciprocal_matches(forward: list[AlignmentBlock], reverse: list[AlignmentBlock]) -> list[tuple[str, str]]:
    matches: list[tuple[str, str]] = []
    for left in (block for block in forward if block.qualifying):
        for right in (block for block in reverse if block.qualifying):
            if left.query_name != right.target_name or left.target_name != right.query_name:
                continue
            query_overlap = _overlap(left.query_start, left.query_end, right.target_start, right.target_end)
            target_overlap = _overlap(left.target_start, left.target_end, right.query_start, right.query_end)
            query_min = min(left.query_end - left.query_start, right.target_end - right.target_start)
            target_min = min(left.target_end - left.target_start, right.query_end - right.query_start)
            if query_min and target_min and query_overlap / query_min >= 0.5 and target_overlap / target_min >= 0.5:
                left.reciprocal = True
                right.reciprocal = True
                matches.append((left.block_id, right.block_id))
    return matches


def _merged_bases(intervals: Iterable[tuple[str, int, int]]) -> int:
    grouped: dict[str, list[tuple[int, int]]] = {}
    for contig, start, end in intervals:
        grouped.setdefault(contig, []).append((start, end))
    total = 0
    for contig_intervals in grouped.values():
        current_start = current_end = None
        for start, end in sorted(contig_intervals):
            if current_start is None:
                current_start, current_end = start, end
            elif start <= current_end:
                current_end = max(current_end, end)
            else:
                total += current_end - current_start
                current_start, current_end = start, end
        if current_start is not None:
            total += current_end - current_start
    return total


def run_alignment(
    source: ReferenceInput,
    target: ReferenceInput,
    stage: Path,
    preset: str,
    threads: int,
    minimum_length: int,
    minimum_identity: float,
    minimum_mapq: int,
) -> tuple[Path, list[AlignmentBlock], str]:
    direction = f"{source.accession}_to_{target.accession}"
    paf_path = stage / "alignments" / f"{direction}.paf"
    log_path = stage / "logs" / f"{direction}.stderr.txt"
    command = [
        "minimap2",
        "-c",
        "-x",
        preset,
        "--cs=long",
        "--eqx",
        "--secondary=yes",
        "-N",
        "5",
        "-t",
        str(threads),
        str(target.reference_fasta),
        str(source.reference_fasta),
    ]
    with paf_path.open("w", encoding="utf-8", newline="") as stdout_handle:
        result = subprocess.run(command, stdout=stdout_handle, stderr=subprocess.PIPE, text=True)
    sanitized_stderr = result.stderr.replace(str(source.reference_fasta), "<SOURCE_FASTA>")
    sanitized_stderr = sanitized_stderr.replace(str(target.reference_fasta), "<TARGET_FASTA>")
    log_path.write_text(sanitized_stderr, encoding="utf-8")
    if result.returncode != 0:
        raise ReviewError(f"minimap2 failed for {direction}; see {log_path}")
    blocks = parse_paf(
        paf_path,
        direction,
        source.accession,
        target.accession,
        minimum_length,
        minimum_identity,
        minimum_mapq,
    )
    portable = (
        f"minimap2 -cx {preset} --cs=long --eqx --secondary=yes -N 5 -t {threads} "
        f"<TARGET_FASTA:{target.accession}> <SOURCE_FASTA:{source.accession}> > <PAF:{direction}>"
    )
    return paf_path, blocks, portable


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    if not rows and fieldnames is None:
        raise ReviewError(f"cannot infer columns for empty CSV: {path}")
    columns = fieldnames or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _block_row(block: AlignmentBlock) -> dict[str, object]:
    return {
        "block_id": block.block_id,
        "direction": block.direction,
        "source_accession": block.source_accession,
        "target_accession": block.target_accession,
        "query_name": block.query_name,
        "query_length": block.query_length,
        "query_start": block.query_start,
        "query_end": block.query_end,
        "strand": block.strand,
        "target_name": block.target_name,
        "target_length": block.target_length,
        "target_start": block.target_start,
        "target_end": block.target_end,
        "matches": block.matches,
        "alignment_bases": block.alignment_bases,
        "identity": f"{block.identity:.8f}",
        "mapq": block.mapq,
        "alignment_type": block.alignment_type,
        "chain_score": "" if block.chain_score is None else block.chain_score,
        "secondary_score": "" if block.secondary_score is None else block.secondary_score,
        "qualifying": str(block.qualifying).lower(),
        "reciprocal": str(block.reciprocal).lower(),
    }


def connection_review(args: argparse.Namespace) -> None:
    output_dir = args.output_dir.resolve()
    alignment_dir = args.alignment_dir.resolve()
    if output_dir.exists() or alignment_dir.exists():
        raise ReviewError("output and alignment directories must not already exist")
    if args.threads < 1:
        raise ReviewError("threads must be at least one")
    if args.minimum_alignment_length < 1:
        raise ReviewError("minimum alignment length must be positive")
    if not 0 < args.minimum_alignment_identity <= 1:
        raise ReviewError("minimum alignment identity must be in (0, 1]")
    if not 0 <= args.minimum_alignment_mapq <= 255:
        raise ReviewError("minimum alignment MAPQ must be between 0 and 255")
    if not 0 < args.minimum_comparable_coverage <= 1:
        raise ReviewError("minimum comparable coverage must be in (0, 1]")

    references = load_reference_manifest(args.reference_manifest.resolve(), args.baseline_accession)
    versions = tool_versions()
    validation_rows, inventories, source_hashes_before = verify_reference_inputs(references, args.sample_id)
    reproduction_rows, project15_counts, project15_hashes = reproduce_project15(
        args.project15_comparison_dir.resolve(),
        args.project15_orthology_dir.resolve(),
        args.project15_annotation_dir.resolve(),
    )
    accepted_total = sum(int(item["accepted_variants"]) for item in inventories.values())
    if accepted_total != EXPECTED_PROJECT15_COUNTS["accepted_variants"]:
        raise ReviewError(f"accepted VCF total is {accepted_total}, expected 357581")

    run_id = uuid.uuid4().hex
    output_stage = output_dir.parent / f".{output_dir.name}.staging-{run_id}"
    alignment_stage = alignment_dir.parent / f".{alignment_dir.name}.staging-{run_id}"
    output_stage.mkdir(parents=True)
    (output_stage / "tables").mkdir()
    (output_stage / "reports").mkdir()
    (alignment_stage / "alignments").mkdir(parents=True)
    (alignment_stage / "logs").mkdir()
    baseline = next(reference for reference in references if reference.role == "baseline")
    alternatives = [reference for reference in references if reference.role == "alternative"]
    all_blocks: list[AlignmentBlock] = []
    alignment_sets: dict[str, list[AlignmentBlock]] = {}
    alignment_paths: dict[str, Path] = {}
    portable_commands: list[str] = []
    coverage_rows: list[dict[str, object]] = []
    try:
        for alternative in alternatives:
            for source, target in ((alternative, baseline), (baseline, alternative)):
                paf_path, blocks, portable = run_alignment(
                    source,
                    target,
                    alignment_stage,
                    args.alignment_preset,
                    args.threads,
                    args.minimum_alignment_length,
                    args.minimum_alignment_identity,
                    args.minimum_alignment_mapq,
                )
                direction = f"{source.accession}_to_{target.accession}"
                alignment_sets[direction] = blocks
                alignment_paths[direction] = paf_path
                all_blocks.extend(blocks)
                portable_commands.append(portable)

            forward_name = f"{alternative.accession}_to_{baseline.accession}"
            reverse_name = f"{baseline.accession}_to_{alternative.accession}"
            forward = alignment_sets[forward_name]
            reverse = alignment_sets[reverse_name]
            reciprocal_pairs = _reciprocal_matches(forward, reverse)
            baseline_intervals: list[tuple[str, int, int]] = []
            alternative_intervals: list[tuple[str, int, int]] = []
            for block in forward:
                if block.reciprocal:
                    alternative_intervals.append((block.query_name, block.query_start, block.query_end))
                    baseline_intervals.append((block.target_name, block.target_start, block.target_end))
            for block in reverse:
                if block.reciprocal:
                    baseline_intervals.append((block.query_name, block.query_start, block.query_end))
                    alternative_intervals.append((block.target_name, block.target_start, block.target_end))
            baseline_covered = _merged_bases(baseline_intervals)
            alternative_covered = _merged_bases(alternative_intervals)
            baseline_bases = int(inventories[baseline.accession]["bases"])
            alternative_bases = int(inventories[alternative.accession]["bases"])
            baseline_fraction = baseline_covered / baseline_bases
            alternative_fraction = alternative_covered / alternative_bases
            coverage_rows.append(
                {
                    "baseline_accession": baseline.accession,
                    "alternative_accession": alternative.accession,
                    "alternative_label": alternative.label,
                    "forward_paf_records": len(forward),
                    "reverse_paf_records": len(reverse),
                    "forward_qualifying_blocks": sum(block.qualifying for block in forward),
                    "reverse_qualifying_blocks": sum(block.qualifying for block in reverse),
                    "reciprocal_block_pairs": len(reciprocal_pairs),
                    "baseline_bases": baseline_bases,
                    "baseline_reciprocal_bases": baseline_covered,
                    "baseline_reciprocal_percent": f"{baseline_fraction * 100:.8f}",
                    "alternative_bases": alternative_bases,
                    "alternative_reciprocal_bases": alternative_covered,
                    "alternative_reciprocal_percent": f"{alternative_fraction * 100:.8f}",
                    "substantive_comparison": str(
                        baseline_fraction >= args.minimum_comparable_coverage
                        and alternative_fraction >= args.minimum_comparable_coverage
                    ).lower(),
                }
            )

        source_hashes_after: dict[str, str] = {}
        for reference in references:
            files = {
                "reference_fasta": reference.reference_fasta,
                "accepted_vcf": reference.accepted_vcf,
                "evaluated_vcf": reference.evaluated_vcf,
                "rejected_vcf": reference.rejected_vcf,
                "callable_bed": reference.callable_bed,
                "accepted_vcf_index": Path(f"{reference.accepted_vcf}.csi"),
                "evaluated_vcf_index": Path(f"{reference.evaluated_vcf}.csi"),
                "rejected_vcf_index": Path(f"{reference.rejected_vcf}.csi"),
            }
            for kind, path in files.items():
                source_hashes_after[f"{reference.accession}:{kind}"] = sha256(path)
        if source_hashes_after != source_hashes_before:
            raise ReviewError("one or more inherited source files changed during alignment")

        write_csv(output_stage / "tables" / "source_validation.csv", validation_rows)
        write_csv(output_stage / "tables" / "project15_reproduction.csv", reproduction_rows)
        write_csv(output_stage / "tables" / "alignment_blocks.csv", [_block_row(block) for block in all_blocks])
        write_csv(output_stage / "tables" / "reciprocal_coverage.csv", coverage_rows)
        write_csv(
            output_stage / "reports" / "tool_versions.csv",
            [{"tool": tool, "version": version} for tool, version in versions.items()],
        )
        alignment_hash_rows = [
            {
                "direction": direction,
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
                "paf_records": len(alignment_sets[direction]),
            }
            for direction, path in alignment_paths.items()
        ]
        write_csv(output_stage / "reports" / "alignment_checksums.csv", alignment_hash_rows)
        (output_stage / "portable_commands.txt").write_text(
            "\n".join(portable_commands) + "\n", encoding="utf-8"
        )
        manifest = {
            "workflow": "variant-liftover-explorer",
            "workflow_version": "0.1.0",
            "milestone": "connection_review",
            "status": "pass",
            "baseline_accession": args.baseline_accession,
            "sample_id": args.sample_id,
            "parameters": {
                "alignment_preset": args.alignment_preset,
                "minimum_alignment_length": args.minimum_alignment_length,
                "minimum_alignment_identity": args.minimum_alignment_identity,
                "minimum_alignment_mapq": args.minimum_alignment_mapq,
                "minimum_comparable_coverage": args.minimum_comparable_coverage,
                "threads": args.threads,
                "reciprocal_overlap_fraction": 0.5,
            },
            "project15_reproduction": {
                **project15_counts,
                "references": EXPECTED_PROJECT15_COUNTS["references"],
            },
            "reference_inventory": inventories,
            "alignment_counts": {
                "directions": len(alignment_paths),
                "paf_records": len(all_blocks),
                "qualifying_blocks": sum(block.qualifying for block in all_blocks),
                "reciprocal_blocks": sum(block.reciprocal for block in all_blocks),
                "substantive_alternatives": sum(
                    row["substantive_comparison"] == "true" for row in coverage_rows
                ),
            },
            "tools": versions,
            "source_hashes": source_hashes_before,
            "project15_hashes": project15_hashes,
            "alignment_hashes": {row["direction"]: row["sha256"] for row in alignment_hash_rows},
            "source_inputs_unchanged": True,
            "interpretation_boundary": (
                "Reciprocal whole-genome coverage establishes coordinate feasibility; "
                "it does not yet establish exact variant equivalence."
            ),
        }
        (output_stage / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(alignment_stage, alignment_dir)
        try:
            os.replace(output_stage, output_dir)
        except Exception:
            shutil.rmtree(alignment_dir, ignore_errors=True)
            raise
    except Exception:
        shutil.rmtree(output_stage, ignore_errors=True)
        shutil.rmtree(alignment_stage, ignore_errors=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--project15-comparison-dir", type=Path, required=True)
    parser.add_argument("--project15-orthology-dir", type=Path, required=True)
    parser.add_argument("--project15-annotation-dir", type=Path, required=True)
    parser.add_argument("--baseline-accession", default="GCF_000005845.2")
    parser.add_argument("--sample-id", default="SRR13921545")
    parser.add_argument("--alignment-preset", default="asm20")
    parser.add_argument("--minimum-alignment-length", type=int, default=10_000)
    parser.add_argument("--minimum-alignment-identity", type=float, default=0.85)
    parser.add_argument("--minimum-alignment-mapq", type=int, default=20)
    parser.add_argument("--minimum-comparable-coverage", type=float, default=0.70)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        connection_review(args)
    except (OSError, subprocess.SubprocessError, ReviewError, ValueError, json.JSONDecodeError) as error:
        print(f"connection review failed: {error}", file=sys.stderr)
        return 1
    print(f"connection review published to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

