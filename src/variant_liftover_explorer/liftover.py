"""Conservative, allele-aware variant projection through Minimap2 PAF alignments."""

from __future__ import annotations

import csv
import gzip
import re
import subprocess
from bisect import bisect_right
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence, TextIO


CS_PATTERN = re.compile(
    r"=[A-Za-z]+|:[0-9]+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+|~[A-Za-z]{2}[0-9]+[A-Za-z]{2}"
)
CIGAR_PATTERN = re.compile(r"[0-9]+[MIDNSHP=X]")
DNA_PATTERN = re.compile(r"^[ACGTN]+$")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

FINAL_STATUSES = {
    "accepted_match",
    "rejected_match",
    "target_reference_match",
    "callable_no_candidate",
    "target_not_callable",
    "ambiguous_mapping",
    "nonreciprocal_mapping",
    "nonalignable",
    "complex_alignment_context",
}


class LiftoverError(RuntimeError):
    """Raised when an input or external normalization result is invalid."""


@dataclass(frozen=True)
class Variant:
    """A normalized biallelic small variant using one-based VCF coordinates."""

    identifier: str
    contig: str
    position: int
    reference: str
    alternate: str

    def __post_init__(self) -> None:
        if self.position < 1:
            raise LiftoverError("variant position must be positive")
        if not self.reference or not self.alternate:
            raise LiftoverError("variant alleles must be nonempty")
        if not DNA_PATTERN.fullmatch(self.reference) or not DNA_PATTERN.fullmatch(self.alternate):
            raise LiftoverError("only literal DNA alleles are supported")
        if self.reference == self.alternate:
            raise LiftoverError("reference and alternate alleles must differ")

    @property
    def key(self) -> str:
        return f"{self.contig}:{self.position}:{self.reference}:{self.alternate}"

    @property
    def start0(self) -> int:
        return self.position - 1

    @property
    def end0(self) -> int:
        return self.start0 + len(self.reference)

    @property
    def variant_type(self) -> str:
        if len(self.reference) == len(self.alternate) == 1:
            return "snp"
        if len(self.reference) < len(self.alternate):
            return "insertion"
        if len(self.reference) > len(self.alternate):
            return "deletion"
        return "substitution"


@dataclass(frozen=True)
class MapSegment:
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    strand: str

    def contains(self, query_position: int) -> bool:
        return self.query_start <= query_position < self.query_end

    def map_position(self, query_position: int) -> int:
        if not self.contains(query_position):
            raise LiftoverError("query position falls outside mapping segment")
        if self.strand == "+":
            return self.target_start + query_position - self.query_start
        return self.target_start + self.query_end - 1 - query_position


@dataclass(frozen=True)
class PafAlignment:
    alignment_id: str
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
    chain_score: int
    secondary_score: int
    cigar: str
    cs: str
    segments: tuple[MapSegment, ...]
    segment_starts: tuple[int, ...]
    gap_boundaries: frozenset[int]

    @property
    def identity(self) -> float:
        return self.matches / self.alignment_bases if self.alignment_bases else 0.0

    def contains_span(self, contig: str, start: int, end: int) -> bool:
        return contig == self.query_name and self.query_start <= start and end <= self.query_end

    def map_span(self, start: int, end: int) -> tuple[int, ...] | None:
        mapped: list[int] = []
        for query_position in range(start, end):
            index = bisect_right(self.segment_starts, query_position) - 1
            if index < 0 or not self.segments[index].contains(query_position):
                return None
            segment = self.segments[index]
            mapped.append(segment.map_position(query_position))
        return tuple(mapped)


@dataclass(frozen=True)
class ProjectionDraft:
    source: Variant
    status: str
    reason: str
    alignment_id: str = ""
    strand: str = ""
    target_contig: str = ""
    target_position: int | None = None
    target_reference: str = ""
    target_alternate: str = ""
    source_span_start0: int | None = None
    source_span_end0: int | None = None
    target_span_start0: int | None = None
    target_span_end0: int | None = None
    normalized: bool = False
    reciprocal: bool = False

    @property
    def target_key(self) -> str:
        if (
            self.target_position is None
            or not self.target_contig
            or self.target_reference == self.target_alternate
        ):
            return ""
        return (
            f"{self.target_contig}:{self.target_position}:"
            f"{self.target_reference}:{self.target_alternate}"
        )


@dataclass(frozen=True)
class ProjectionResult:
    source: Variant
    status: str
    reason: str
    alignment_id: str
    strand: str
    target_contig: str
    target_position: int | None
    target_reference: str
    target_alternate: str
    reciprocal: bool
    normalized: bool

    @property
    def target_key(self) -> str:
        if self.target_position is None or self.target_reference == self.target_alternate:
            return ""
        return (
            f"{self.target_contig}:{self.target_position}:"
            f"{self.target_reference}:{self.target_alternate}"
        )


@dataclass(frozen=True)
class ProjectionConfig:
    minimum_alignment_length: int = 10_000
    minimum_alignment_identity: float = 0.85
    minimum_alignment_mapq: int = 20
    ambiguity_score_fraction: float = 0.95

    def validate(self) -> None:
        if self.minimum_alignment_length < 1:
            raise LiftoverError("minimum alignment length must be positive")
        if not 0 < self.minimum_alignment_identity <= 1:
            raise LiftoverError("minimum alignment identity must be in (0, 1]")
        if not 0 <= self.minimum_alignment_mapq <= 255:
            raise LiftoverError("minimum alignment MAPQ must be between 0 and 255")
        if not 0 < self.ambiguity_score_fraction <= 1:
            raise LiftoverError("ambiguity score fraction must be in (0, 1]")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1].upper()


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="ascii") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current or current in sequences:
                    raise LiftoverError(f"invalid or duplicate FASTA ID at {path}:{line_number}")
                sequences[current] = []
            elif current is None:
                raise LiftoverError(f"sequence before FASTA header at {path}:{line_number}")
            else:
                value = line.upper()
                if not DNA_PATTERN.fullmatch(value):
                    raise LiftoverError(f"invalid FASTA sequence at {path}:{line_number}")
                sequences[current].append(value)
    joined = {name: "".join(parts) for name, parts in sequences.items()}
    if not joined or any(not sequence for sequence in joined.values()):
        raise LiftoverError(f"empty FASTA record in {path}")
    return joined


def parse_cs(cs: str) -> list[tuple[str, str]]:
    tokens = CS_PATTERN.findall(cs)
    if not tokens or "".join(tokens) != cs:
        raise LiftoverError("invalid long cs tag")
    operations: list[tuple[str, str]] = []
    for token in tokens:
        code, payload = token[0], token[1:]
        if code == "~":
            raise LiftoverError("splice operations are invalid for genome liftover")
        operations.append((code, payload))
    return operations


def _tag(fields: Sequence[str], name: str, default: str = "") -> str:
    prefix = f"{name}:"
    for field in fields:
        if field.startswith(prefix):
            return field.split(":", 2)[2]
    return default


def _segments_from_cs(
    query_start: int,
    query_end: int,
    target_start: int,
    strand: str,
    cs: str,
) -> tuple[tuple[MapSegment, ...], frozenset[int], int, int]:
    query_cursor = query_start if strand == "+" else query_end
    target_cursor = target_start
    segments: list[MapSegment] = []
    gap_boundaries: set[int] = set()
    query_consumed = 0
    target_consumed = 0
    for code, payload in parse_cs(cs):
        if code == "=":
            length = len(payload)
            query_length = target_length = length
        elif code == ":":
            length = int(payload)
            query_length = target_length = length
        elif code == "*":
            query_length = target_length = 1
        elif code == "+":
            query_length, target_length = len(payload), 0
        elif code == "-":
            query_length, target_length = 0, len(payload)
        else:
            raise LiftoverError(f"unsupported cs operation: {code}")

        if query_length and target_length:
            if strand == "+":
                segment_query_start = query_cursor
                segment_query_end = query_cursor + query_length
                query_cursor += query_length
            else:
                segment_query_start = query_cursor - query_length
                segment_query_end = query_cursor
                query_cursor -= query_length
            segments.append(
                MapSegment(
                    query_start=segment_query_start,
                    query_end=segment_query_end,
                    target_start=target_cursor,
                    target_end=target_cursor + target_length,
                    strand=strand,
                )
            )
            target_cursor += target_length
        elif query_length:
            if strand == "+":
                gap_start = query_cursor
                gap_end = query_cursor + query_length
                query_cursor = gap_end
            else:
                gap_start = query_cursor - query_length
                gap_end = query_cursor
                query_cursor = gap_start
            gap_boundaries.update({gap_start, gap_end})
        elif target_length:
            gap_boundaries.add(query_cursor)
            target_cursor += target_length
        query_consumed += query_length
        target_consumed += target_length
    expected_query_cursor = query_end if strand == "+" else query_start
    if query_cursor != expected_query_cursor:
        raise LiftoverError("cs query consumption disagrees with PAF coordinates")
    ordered_segments = tuple(sorted(segments, key=lambda item: item.query_start))
    return ordered_segments, frozenset(gap_boundaries), query_consumed, target_consumed


def _cigar_consumption(cigar: str) -> tuple[int, int]:
    """Return query and target consumption for an alignment CIGAR."""
    query = target = 0
    for token in CIGAR_PATTERN.findall(cigar):
        length, operation = int(token[:-1]), token[-1]
        if operation in "MIS=X":
            query += length
        if operation in "MDN=X":
            target += length
    return query, target


def parse_paf(path: Path) -> list[PafAlignment]:
    alignments: list[PafAlignment] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise LiftoverError(f"malformed PAF row at {path}:{line_number}")
            cigar = _tag(fields[12:], "cg")
            cs = _tag(fields[12:], "cs")
            if not cigar or not cs:
                raise LiftoverError(f"PAF row lacks cg or cs at {path}:{line_number}")
            if "".join(CIGAR_PATTERN.findall(cigar)) != cigar:
                raise LiftoverError(f"invalid CIGAR at {path}:{line_number}")
            if fields[4] not in {"+", "-"}:
                raise LiftoverError(f"invalid PAF strand at {path}:{line_number}")
            try:
                query_start, query_end = int(fields[2]), int(fields[3])
                target_start, target_end = int(fields[7]), int(fields[8])
                query_length, target_length = int(fields[1]), int(fields[6])
                matches, alignment_bases, mapq = (
                    int(fields[9]),
                    int(fields[10]),
                    int(fields[11]),
                )
                chain_score = int(_tag(fields[12:], "s1", "0"))
                secondary_score = int(_tag(fields[12:], "s2", "0"))
            except ValueError as error:
                raise LiftoverError(
                    f"non-integer PAF field at {path}:{line_number}"
                ) from error
            if query_length < 1 or target_length < 1:
                raise LiftoverError(f"invalid PAF sequence length at {path}:{line_number}")
            if not 0 <= query_start < query_end <= query_length:
                raise LiftoverError(f"invalid query coordinates at {path}:{line_number}")
            if not 0 <= target_start < target_end <= target_length:
                raise LiftoverError(f"invalid target coordinates at {path}:{line_number}")
            if alignment_bases < 1 or not 0 <= matches <= alignment_bases:
                raise LiftoverError(f"invalid PAF alignment counts at {path}:{line_number}")
            if not 0 <= mapq <= 255:
                raise LiftoverError(f"invalid PAF MAPQ at {path}:{line_number}")
            alignment_type = _tag(fields[12:], "tp")
            if alignment_type not in {"P", "S"}:
                raise LiftoverError(f"invalid or missing PAF tp tag at {path}:{line_number}")
            if chain_score < 0 or secondary_score < 0:
                raise LiftoverError(f"negative PAF alignment score at {path}:{line_number}")
            segments, gap_boundaries, query_consumed, target_consumed = _segments_from_cs(
                query_start,
                query_end,
                target_start,
                fields[4],
                cs,
            )
            if query_consumed != query_end - query_start:
                raise LiftoverError(f"cs query length disagrees at {path}:{line_number}")
            if target_consumed != target_end - target_start:
                raise LiftoverError(f"cs target length disagrees at {path}:{line_number}")
            cigar_query, cigar_target = _cigar_consumption(cigar)
            if cigar_query != query_consumed or cigar_target != target_consumed:
                raise LiftoverError(f"CIGAR and cs consumption disagree at {path}:{line_number}")
            alignments.append(
                PafAlignment(
                    alignment_id=f"{path.name}:{line_number}",
                    query_name=fields[0],
                    query_length=query_length,
                    query_start=query_start,
                    query_end=query_end,
                    strand=fields[4],
                    target_name=fields[5],
                    target_length=target_length,
                    target_start=target_start,
                    target_end=target_end,
                    matches=matches,
                    alignment_bases=alignment_bases,
                    mapq=mapq,
                    alignment_type=alignment_type,
                    chain_score=chain_score,
                    secondary_score=secondary_score,
                    cigar=cigar,
                    cs=cs,
                    segments=segments,
                    segment_starts=tuple(segment.query_start for segment in segments),
                    gap_boundaries=gap_boundaries,
                )
            )
    if not alignments:
        raise LiftoverError(f"PAF contains no records: {path}")
    return alignments


def validate_alignment_sequences(
    alignment: PafAlignment,
    query_sequences: Mapping[str, str],
    target_sequences: Mapping[str, str],
) -> None:
    if alignment.query_name not in query_sequences or alignment.target_name not in target_sequences:
        raise LiftoverError(f"alignment {alignment.alignment_id} references an unknown contig")
    query = query_sequences[alignment.query_name]
    target = target_sequences[alignment.target_name]
    if len(query) != alignment.query_length or len(target) != alignment.target_length:
        raise LiftoverError(f"alignment {alignment.alignment_id} length disagrees with FASTA")
    oriented_query = query[alignment.query_start : alignment.query_end]
    if alignment.strand == "-":
        oriented_query = reverse_complement(oriented_query)
    query_cursor = target_cursor = 0
    target_slice = target[alignment.target_start : alignment.target_end]
    for code, payload in parse_cs(alignment.cs):
        if code == "=":
            length = len(payload)
            if target_slice[target_cursor : target_cursor + length].upper() != payload.upper():
                raise LiftoverError(f"cs target match disagrees in {alignment.alignment_id}")
            if oriented_query[query_cursor : query_cursor + length].upper() != payload.upper():
                raise LiftoverError(f"cs query match disagrees in {alignment.alignment_id}")
            query_cursor += length
            target_cursor += length
        elif code == ":":
            length = int(payload)
            if (
                target_slice[target_cursor : target_cursor + length].upper()
                != oriented_query[query_cursor : query_cursor + length].upper()
            ):
                raise LiftoverError(f"short cs match disagrees in {alignment.alignment_id}")
            query_cursor += length
            target_cursor += length
        elif code == "*":
            if target_slice[target_cursor].lower() != payload[0].lower():
                raise LiftoverError(f"cs substitution target disagrees in {alignment.alignment_id}")
            if oriented_query[query_cursor].lower() != payload[1].lower():
                raise LiftoverError(f"cs substitution query disagrees in {alignment.alignment_id}")
            query_cursor += 1
            target_cursor += 1
        elif code == "+":
            length = len(payload)
            if oriented_query[query_cursor : query_cursor + length].lower() != payload.lower():
                raise LiftoverError(f"cs insertion disagrees in {alignment.alignment_id}")
            query_cursor += length
        elif code == "-":
            length = len(payload)
            if target_slice[target_cursor : target_cursor + length].lower() != payload.lower():
                raise LiftoverError(f"cs deletion disagrees in {alignment.alignment_id}")
            target_cursor += length
    if query_cursor != len(oriented_query) or target_cursor != len(target_slice):
        raise LiftoverError(f"cs sequence consumption disagrees in {alignment.alignment_id}")


def _select_alignment(
    variant: Variant,
    alignments: Sequence[PafAlignment],
    config: ProjectionConfig,
) -> tuple[PafAlignment | None, str]:
    config.validate()
    containing = [
        alignment
        for alignment in alignments
        if alignment.contains_span(variant.contig, variant.start0, variant.end0)
        and alignment.alignment_bases >= config.minimum_alignment_length
        and alignment.identity >= config.minimum_alignment_identity
    ]
    if not containing:
        return None, "nonalignable"
    primaries = [
        alignment
        for alignment in containing
        if alignment.alignment_type == "P" and alignment.mapq >= config.minimum_alignment_mapq
    ]
    score_sorted = sorted(containing, key=lambda item: (item.chain_score, item.mapq), reverse=True)
    if not primaries:
        if len(score_sorted) > 1 and score_sorted[1].chain_score >= config.ambiguity_score_fraction * score_sorted[0].chain_score:
            return None, "ambiguous_mapping"
        return None, "nonalignable"
    best = max(primaries, key=lambda item: (item.chain_score, item.mapq, item.alignment_bases))
    competitors = [
        alignment
        for alignment in containing
        if alignment.alignment_id != best.alignment_id
        and alignment.chain_score >= config.ambiguity_score_fraction * best.chain_score
    ]
    if competitors:
        return None, "ambiguous_mapping"
    return best, ""


def _reference_matches(variant: Variant, sequences: Mapping[str, str]) -> None:
    sequence = sequences.get(variant.contig)
    if sequence is None:
        raise LiftoverError(f"variant {variant.identifier} uses unknown contig {variant.contig}")
    if variant.end0 > len(sequence):
        raise LiftoverError(
            f"variant {variant.identifier} extends beyond its source contig; "
            "circular-origin alleles must be split and normalized before liftover"
        )
    observed = sequence[variant.start0 : variant.end0].upper()
    if observed != variant.reference.upper():
        raise LiftoverError(
            f"variant {variant.identifier} REF {variant.reference} disagrees with source FASTA {observed}"
        )


def project_once(
    variant: Variant,
    source_sequences: Mapping[str, str],
    target_sequences: Mapping[str, str],
    alignments: Sequence[PafAlignment],
    config: ProjectionConfig,
) -> ProjectionDraft:
    _reference_matches(variant, source_sequences)
    alignment, failure = _select_alignment(variant, alignments, config)
    if alignment is None:
        return ProjectionDraft(source=variant, status=failure, reason=failure.replace("_", " "))
    mapped = alignment.map_span(variant.start0, variant.end0)
    if mapped is None or len(set(mapped)) != len(mapped):
        return ProjectionDraft(
            source=variant,
            status="complex_alignment_context",
            reason="source REF span crosses an alignment gap",
            alignment_id=alignment.alignment_id,
            strand=alignment.strand,
        )
    ordered = sorted(mapped)
    if any(right - left != 1 for left, right in zip(ordered, ordered[1:])):
        return ProjectionDraft(
            source=variant,
            status="complex_alignment_context",
            reason="source REF span is not contiguous in the target",
            alignment_id=alignment.alignment_id,
            strand=alignment.strand,
        )
    if variant.variant_type == "insertion" and variant.end0 in alignment.gap_boundaries:
        return ProjectionDraft(
            source=variant,
            status="complex_alignment_context",
            reason="insertion anchor coincides with a reference-alignment gap",
            alignment_id=alignment.alignment_id,
            strand=alignment.strand,
        )
    target_start = ordered[0]
    target_end = ordered[-1] + 1
    target_sequence = target_sequences.get(alignment.target_name)
    if target_sequence is None:
        raise LiftoverError(f"alignment target contig is absent: {alignment.target_name}")
    target_reference = target_sequence[target_start:target_end].upper()
    if len(variant.reference) == len(variant.alternate):
        target_alternate = (
            variant.alternate.upper()
            if alignment.strand == "+"
            else reverse_complement(variant.alternate)
        )
    elif variant.variant_type == "insertion":
        if len(variant.reference) != 1 or not variant.alternate.startswith(variant.reference):
            return ProjectionDraft(
                source=variant,
                status="complex_alignment_context",
                reason="insertion is not represented with one retained source anchor",
                alignment_id=alignment.alignment_id,
                strand=alignment.strand,
            )
        inserted = variant.alternate[len(variant.reference) :]
        target_alternate = (
            target_reference + inserted.upper()
            if alignment.strand == "+"
            else reverse_complement(inserted) + target_reference
        )
    elif variant.variant_type == "deletion":
        if len(variant.alternate) != 1 or not variant.reference.startswith(variant.alternate):
            return ProjectionDraft(
                source=variant,
                status="complex_alignment_context",
                reason="deletion is not represented with one retained source anchor",
                alignment_id=alignment.alignment_id,
                strand=alignment.strand,
            )
        retained_target_position = mapped[0]
        target_alternate = target_sequence[retained_target_position].upper()
    else:
        return ProjectionDraft(
            source=variant,
            status="complex_alignment_context",
            reason="unsupported complex allele",
            alignment_id=alignment.alignment_id,
            strand=alignment.strand,
        )
    status = "target_reference_match" if target_reference == target_alternate else "projected"
    return ProjectionDraft(
        source=variant,
        status=status,
        reason=(
            "sample allele equals the target reference"
            if status == "target_reference_match"
            else "projected allele requires normalization"
        ),
        alignment_id=alignment.alignment_id,
        strand=alignment.strand,
        target_contig=alignment.target_name,
        target_position=target_start + 1,
        target_reference=target_reference,
        target_alternate=target_alternate,
        source_span_start0=variant.start0,
        source_span_end0=variant.end0,
        target_span_start0=target_start,
        target_span_end0=target_end,
    )


def normalize_drafts(
    drafts: Sequence[ProjectionDraft],
    target_fasta: Path,
    work_dir: Path,
    bcftools: str = "bcftools",
) -> dict[str, Variant]:
    candidates = [draft for draft in drafts if draft.status == "projected"]
    if not candidates:
        return {}
    sequences = read_fasta(target_fasta)
    work_dir.mkdir(parents=True, exist_ok=True)
    input_vcf = work_dir / "provisional.vcf"
    output_vcf = work_dir / "normalized.vcf"
    with input_vcf.open("w", encoding="utf-8", newline="") as handle:
        handle.write("##fileformat=VCFv4.3\n")
        for contig, sequence in sequences.items():
            handle.write(f"##contig=<ID={contig},length={len(sequence)}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for index, draft in enumerate(candidates, start=1):
            handle.write(
                f"{draft.target_contig}\t{draft.target_position}\tP{index}\t"
                f"{draft.target_reference}\t{draft.target_alternate}\t.\tPASS\t.\n"
            )
    command = [
        bcftools,
        "norm",
        "--no-version",
        "-f",
        str(target_fasta),
        "-m",
        "-any",
        "-Ov",
        "-o",
        str(output_vcf),
        str(input_vcf),
    ]
    try:
        result = subprocess.run(command, capture_output=True, text=True)
    except OSError as error:
        raise LiftoverError(f"bcftools normalization could not start: {error}") from error
    if result.returncode != 0:
        raise LiftoverError(f"bcftools normalization failed: {result.stderr.strip()}")
    normalized: dict[str, Variant] = {}
    with output_vcf.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise LiftoverError("bcftools normalization returned a malformed VCF row")
            identifier = fields[2]
            if identifier in normalized:
                raise LiftoverError(f"normalization split projection {identifier}")
            normalized[identifier] = Variant(
                identifier=identifier,
                contig=fields[0],
                position=int(fields[1]),
                reference=fields[3].upper(),
                alternate=fields[4].upper(),
            )
    expected_identifiers = {f"P{index}" for index in range(1, len(candidates) + 1)}
    if set(normalized) != expected_identifiers:
        raise LiftoverError("normalization did not preserve one row per projection")
    return normalized


def _callable_positions(variant: Variant) -> range:
    # VCF indels are evaluated over their represented REF span, including the anchor.
    return range(variant.start0, variant.end0)


def is_callable(
    variant: Variant,
    callable_intervals: Mapping[str, Sequence[tuple[int, int]]],
    target_sequences: Mapping[str, str],
) -> bool:
    sequence = target_sequences.get(variant.contig)
    if sequence is None:
        return False
    positions = list(_callable_positions(variant))
    if not positions or positions[-1] >= len(sequence):
        return False
    intervals = callable_intervals.get(variant.contig, ())
    for position in positions:
        index = bisect_right(intervals, (position, float("inf"))) - 1
        if index < 0 or not intervals[index][0] <= position < intervals[index][1]:
            return False
    return True


def _replace_with_normalized(
    drafts: Sequence[ProjectionDraft], normalized: Mapping[str, Variant]
) -> list[ProjectionDraft]:
    output: list[ProjectionDraft] = []
    projected_index = 0
    for draft in drafts:
        if draft.status != "projected":
            output.append(draft)
            continue
        projected_index += 1
        variant = normalized[f"P{projected_index}"]
        output.append(
            replace(
                draft,
                target_contig=variant.contig,
                target_position=variant.position,
                target_reference=variant.reference,
                target_alternate=variant.alternate,
                normalized=True,
            )
        )
    return output


def _reciprocal_coordinate_matches(
    draft: ProjectionDraft,
    reverse_alignments: Sequence[PafAlignment],
    config: ProjectionConfig,
) -> tuple[bool, str]:
    """Validate that every mapped target REF base returns to the source REF span."""
    probe_alternate = "A" if draft.target_reference[0] != "A" else "C"
    probe = Variant(
        identifier=f"C_{draft.source.identifier}",
        contig=draft.target_contig,
        position=int(draft.target_position),
        reference=draft.target_reference,
        alternate=probe_alternate,
    )
    reverse_alignment, failure = _select_alignment(probe, reverse_alignments, config)
    if reverse_alignment is None:
        return False, f"reverse coordinate mapping returned {failure}"
    mapped_back = reverse_alignment.map_span(probe.start0, probe.end0)
    expected = tuple(range(draft.source.start0, draft.source.end0))
    if mapped_back is None or tuple(sorted(mapped_back)) != expected:
        return False, "reverse coordinate mapping did not recover the complete source REF span"
    return True, ""


def project_variants(
    variants: Sequence[Variant],
    source_fasta: Path,
    target_fasta: Path,
    forward_paf: Path,
    reverse_paf: Path,
    accepted_target_keys: set[str],
    rejected_target_keys: set[str],
    callable_intervals: Mapping[str, Sequence[tuple[int, int]]],
    config: ProjectionConfig,
    work_dir: Path,
) -> list[ProjectionResult]:
    source_sequences = read_fasta(source_fasta)
    target_sequences = read_fasta(target_fasta)
    forward = parse_paf(forward_paf)
    reverse = parse_paf(reverse_paf)
    for alignment in forward:
        validate_alignment_sequences(alignment, source_sequences, target_sequences)
    for alignment in reverse:
        validate_alignment_sequences(alignment, target_sequences, source_sequences)

    return project_variants_loaded(
        variants,
        source_sequences,
        target_sequences,
        forward,
        reverse,
        accepted_target_keys,
        rejected_target_keys,
        callable_intervals,
        config,
        source_fasta,
        target_fasta,
        work_dir,
    )


def project_variants_loaded(
    variants: Sequence[Variant],
    source_sequences: Mapping[str, str],
    target_sequences: Mapping[str, str],
    forward: Sequence[PafAlignment],
    reverse: Sequence[PafAlignment],
    accepted_target_keys: set[str],
    rejected_target_keys: set[str],
    callable_intervals: Mapping[str, Sequence[tuple[int, int]]],
    config: ProjectionConfig,
    source_fasta: Path,
    target_fasta: Path,
    work_dir: Path,
) -> list[ProjectionResult]:
    """Project a batch using already parsed and sequence-validated alignments."""

    source_keys = [variant.key for variant in variants]
    if len(source_keys) != len(set(source_keys)):
        raise LiftoverError("source variants contain duplicate normalized keys")
    overlap = accepted_target_keys & rejected_target_keys
    if overlap:
        raise LiftoverError(
            f"accepted and rejected target evidence overlap: {sorted(overlap)[0]}"
        )
    for contig, intervals in callable_intervals.items():
        previous_end = -1
        for start, end in intervals:
            if start < 0 or end <= start or start < previous_end:
                raise LiftoverError(f"callable intervals are invalid or unsorted for {contig}")
            previous_end = end

    drafts = [project_once(variant, source_sequences, target_sequences, forward, config) for variant in variants]
    normalized = normalize_drafts(drafts, target_fasta, work_dir / "forward_norm")
    drafts = _replace_with_normalized(drafts, normalized)
    back_drafts: list[ProjectionDraft] = []
    back_indices: list[int] = []
    back_failures: dict[int, str] = {}
    coordinate_failures: dict[int, str] = {}
    for index, draft in enumerate(drafts):
        if draft.status == "projected":
            reciprocal_variant = Variant(
                identifier=f"R{index + 1}",
                contig=draft.target_contig,
                position=int(draft.target_position),
                reference=draft.target_reference,
                alternate=draft.target_alternate,
            )
            back_draft = project_once(
                reciprocal_variant,
                target_sequences,
                source_sequences,
                reverse,
                config,
            )
            if back_draft.status == "projected":
                back_indices.append(index)
                back_drafts.append(back_draft)
            else:
                back_failures[index] = f"reverse allele projection returned {back_draft.status}"
        elif draft.status == "target_reference_match":
            reciprocal, reason = _reciprocal_coordinate_matches(draft, reverse, config)
            if not reciprocal:
                coordinate_failures[index] = reason

    normalized_back: dict[int, Variant] = {}
    if back_drafts:
        back_values = normalize_drafts(back_drafts, source_fasta, work_dir / "reverse_norm")
        for order, draft_index in enumerate(back_indices, start=1):
            normalized_back[draft_index] = back_values[f"P{order}"]

    results: list[ProjectionResult] = []
    for index, draft in enumerate(drafts):
        if draft.status not in {"projected", "target_reference_match"}:
            results.append(
                ProjectionResult(
                    source=draft.source,
                    status=draft.status,
                    reason=draft.reason,
                    alignment_id=draft.alignment_id,
                    strand=draft.strand,
                    target_contig=draft.target_contig,
                    target_position=draft.target_position,
                    target_reference=draft.target_reference,
                    target_alternate=draft.target_alternate,
                    reciprocal=False,
                    normalized=draft.normalized,
                )
            )
            continue

        if draft.status == "projected":
            if index in back_failures:
                results.append(
                    ProjectionResult(
                        source=draft.source,
                        status="nonreciprocal_mapping",
                        reason=back_failures[index],
                        alignment_id=draft.alignment_id,
                        strand=draft.strand,
                        target_contig=draft.target_contig,
                        target_position=draft.target_position,
                        target_reference=draft.target_reference,
                        target_alternate=draft.target_alternate,
                        reciprocal=False,
                        normalized=True,
                    )
                )
                continue
            back_normalized = normalized_back[index]
            if back_normalized.key != draft.source.key:
                results.append(
                    ProjectionResult(
                        source=draft.source,
                        status="nonreciprocal_mapping",
                        reason=f"reverse normalized key {back_normalized.key} differs from source key",
                        alignment_id=draft.alignment_id,
                        strand=draft.strand,
                        target_contig=draft.target_contig,
                        target_position=draft.target_position,
                        target_reference=draft.target_reference,
                        target_alternate=draft.target_alternate,
                        reciprocal=False,
                        normalized=True,
                    )
                )
                continue
        else:
            if index in coordinate_failures:
                results.append(
                    ProjectionResult(
                        source=draft.source,
                        status="nonreciprocal_mapping",
                        reason=coordinate_failures[index],
                        alignment_id=draft.alignment_id,
                        strand=draft.strand,
                        target_contig=draft.target_contig,
                        target_position=draft.target_position,
                        target_reference=draft.target_reference,
                        target_alternate=draft.target_alternate,
                        reciprocal=False,
                        normalized=False,
                    )
                )
                continue

        if draft.status == "target_reference_match":
            final_status = "target_reference_match"
            reason = "sample allele equals the target reference after reciprocal coordinate validation"
        else:
            target_variant = Variant(
                identifier=draft.source.identifier,
                contig=draft.target_contig,
                position=int(draft.target_position),
                reference=draft.target_reference,
                alternate=draft.target_alternate,
            )
            if target_variant.key in accepted_target_keys:
                final_status = "accepted_match"
                reason = "normalized projected key occurs in the target accepted VCF"
            elif target_variant.key in rejected_target_keys:
                final_status = "rejected_match"
                reason = "normalized projected key occurs in the target rejected VCF"
            elif is_callable(target_variant, callable_intervals, target_sequences):
                final_status = "callable_no_candidate"
                reason = "target locus is callable but no evaluated target candidate matches"
            else:
                final_status = "target_not_callable"
                reason = "required target reference positions are not all callable"
        results.append(
            ProjectionResult(
                source=draft.source,
                status=final_status,
                reason=reason,
                alignment_id=draft.alignment_id,
                strand=draft.strand,
                target_contig=draft.target_contig,
                target_position=draft.target_position,
                target_reference=draft.target_reference,
                target_alternate=draft.target_alternate,
                reciprocal=True,
                normalized=draft.normalized,
            )
        )
    if len(results) != len(variants):
        raise LiftoverError("projection accounting failed")
    if any(result.status not in FINAL_STATUSES for result in results):
        raise LiftoverError("projection produced an unknown final status")
    return results


def read_vcf(path: Path) -> list[Variant]:
    variants: list[Variant] = []
    seen: set[str] = set()
    header_seen = False
    handle: TextIO
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if raw.startswith("#CHROM"):
                if header_seen:
                    raise LiftoverError(f"duplicate #CHROM header in {path}")
                header_seen = True
                continue
            if raw.startswith("#") or not raw.strip():
                continue
            if not header_seen:
                raise LiftoverError(f"VCF data precedes #CHROM header in {path}:{line_number}")
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise LiftoverError(f"malformed VCF row at {path}:{line_number}")
            if "," in fields[4] or fields[4].startswith("<") or "[" in fields[4] or "]" in fields[4]:
                raise LiftoverError(f"non-biallelic or symbolic allele at {path}:{line_number}")
            try:
                position = int(fields[1])
            except ValueError as error:
                raise LiftoverError(f"invalid VCF position at {path}:{line_number}") from error
            identifier = fields[2] if fields[2] != "." else f"V{len(variants) + 1}"
            variant = Variant(identifier, fields[0], position, fields[3].upper(), fields[4].upper())
            if variant.key in seen:
                raise LiftoverError(f"duplicate variant key in {path}: {variant.key}")
            seen.add(variant.key)
            variants.append(variant)
    if not header_seen:
        raise LiftoverError(f"missing #CHROM header in {path}")
    return variants


def read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise LiftoverError(f"malformed BED row at {path}:{line_number}")
            if not fields[0]:
                raise LiftoverError(f"empty BED contig at {path}:{line_number}")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise LiftoverError(f"invalid BED coordinate at {path}:{line_number}") from error
            if start < 0 or end <= start:
                raise LiftoverError(f"invalid BED interval at {path}:{line_number}")
            intervals.setdefault(fields[0], []).append((start, end))
    for contig, values in intervals.items():
        ordered = sorted(values)
        merged: list[tuple[int, int]] = []
        for start, end in ordered:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        intervals[contig] = merged
    return intervals


def write_projection_csv(path: Path, results: Iterable[ProjectionResult]) -> None:
    fieldnames = [
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
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow(
                {
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
            )
