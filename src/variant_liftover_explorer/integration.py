"""Connect exact variant projections to Project 15 gene-level comparisons."""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .liftover import FINAL_STATUSES, LiftoverError


BURDEN_CATEGORIES = (
    "same_exact_variants",
    "partially_shared_variants",
    "reference_allele_difference",
    "callability_or_filter_difference",
    "different_variants_in_same_gene",
    "unresolved_liftover",
)

UNRESOLVED_STATUSES = {
    "ambiguous_mapping",
    "nonreciprocal_mapping",
    "nonalignable",
    "complex_alignment_context",
}

IMPACT_ORDER = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


@dataclass(frozen=True)
class GeneVariant:
    accession: str
    variant_key: str
    gene_id: str
    gene_name: str
    consequence: str
    impact: str
    warnings: str = ""

    def __post_init__(self) -> None:
        if not self.accession or not self.variant_key or not self.gene_id:
            raise LiftoverError("gene-variant records require accession, variant key, and gene ID")


@dataclass(frozen=True)
class GenePair:
    baseline_accession: str
    baseline_gene_id: str
    baseline_gene_name: str
    alternative_accession: str
    alternative_gene_id: str
    alternative_gene_name: str
    comparison_status: str = "comparable"
    comparable: bool = True

    def __post_init__(self) -> None:
        required = (
            self.baseline_accession,
            self.baseline_gene_id,
            self.alternative_accession,
            self.alternative_gene_id,
        )
        if any(not value for value in required):
            raise LiftoverError("gene-pair records require both accessions and gene IDs")


@dataclass(frozen=True)
class ProjectionEvidence:
    alternative_accession: str
    source_variant_key: str
    status: str
    target_variant_key: str = ""
    reason: str = ""

    def __post_init__(self) -> None:
        if self.status not in FINAL_STATUSES:
            raise LiftoverError(f"unknown projection status: {self.status}")
        if self.status in {"accepted_match", "rejected_match"} and not self.target_variant_key:
            raise LiftoverError(f"{self.status} requires a target variant key")


def _unique_index(items: Iterable[object], key_function, label: str) -> dict[object, object]:
    index: dict[object, object] = {}
    for item in items:
        key = key_function(item)
        if key in index:
            raise LiftoverError(f"duplicate {label}: {key}")
        index[key] = item
    return index


def _effect_relation(source: GeneVariant, target: GeneVariant | None) -> str:
    if target is None:
        return "target_gene_effect_missing"
    if source.consequence == target.consequence and source.impact == target.impact:
        return "same_consequence_and_impact"
    if source.impact == target.impact:
        return "same_impact_different_consequence"
    return "different_impact"


def _category(
    baseline_keys: set[str],
    alternative_keys: set[str],
    evidences: Sequence[ProjectionEvidence],
    exact_source_keys: set[str],
    exact_target_keys: set[str],
) -> str:
    statuses = {evidence.status for evidence in evidences}
    # Uncertainty wins so a partial match cannot overstate a gene-pair conclusion.
    if statuses & UNRESOLVED_STATUSES:
        return "unresolved_liftover"
    if exact_source_keys == baseline_keys and exact_target_keys == alternative_keys:
        return "same_exact_variants"
    if exact_source_keys:
        return "partially_shared_variants"
    if "target_reference_match" in statuses:
        return "reference_allele_difference"
    if statuses & {"rejected_match", "target_not_callable"}:
        return "callability_or_filter_difference"
    return "different_variants_in_same_gene"


def integrate_gene_burdens(
    gene_pairs: Sequence[GenePair],
    gene_variants: Sequence[GeneVariant],
    projections: Sequence[ProjectionEvidence],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Return per-variant ortholog evidence and one explanation per comparable pair."""
    comparable_pairs = [pair for pair in gene_pairs if pair.comparable]
    _unique_index(
        comparable_pairs,
        lambda pair: (pair.baseline_gene_id, pair.alternative_accession),
        "comparable gene pair",
    )
    variant_index = _unique_index(
        gene_variants,
        lambda item: (item.accession, item.variant_key, item.gene_id),
        "variant-gene assignment",
    )
    projection_index = _unique_index(
        projections,
        lambda item: (item.alternative_accession, item.source_variant_key),
        "directional projection",
    )
    variants_by_gene: dict[tuple[str, str], dict[str, GeneVariant]] = {}
    for item in variant_index.values():
        variants_by_gene.setdefault((item.accession, item.gene_id), {})[item.variant_key] = item

    equivalence_rows: list[dict[str, object]] = []
    explanation_rows: list[dict[str, object]] = []
    for pair in comparable_pairs:
        baseline = variants_by_gene.get((pair.baseline_accession, pair.baseline_gene_id), {})
        alternative = variants_by_gene.get((pair.alternative_accession, pair.alternative_gene_id), {})
        exact_source_keys: set[str] = set()
        exact_target_keys: set[str] = set()
        pair_evidence: list[ProjectionEvidence] = []
        for source_key in sorted(baseline):
            evidence_key = (pair.alternative_accession, source_key)
            if evidence_key not in projection_index:
                raise LiftoverError(f"missing projection for {source_key} to {pair.alternative_accession}")
            evidence = projection_index[evidence_key]
            pair_evidence.append(evidence)
            # A coordinate-compatible target allele is not an exact equivalent until
            # the reciprocal projection succeeds. Keep target-gene consequences
            # hidden for unresolved or rejected mappings to avoid overstating them.
            target_effect = (
                alternative.get(evidence.target_variant_key)
                if evidence.status == "accepted_match"
                else None
            )
            target_in_ortholog = target_effect is not None
            if evidence.status == "accepted_match" and target_in_ortholog:
                exact_source_keys.add(source_key)
                exact_target_keys.add(evidence.target_variant_key)
            source_effect = baseline[source_key]
            equivalence_rows.append(
                {
                    "baseline_accession": pair.baseline_accession,
                    "baseline_gene_id": pair.baseline_gene_id,
                    "baseline_gene_name": pair.baseline_gene_name,
                    "baseline_variant_key": source_key,
                    "baseline_consequence": source_effect.consequence,
                    "baseline_impact": source_effect.impact,
                    "alternative_accession": pair.alternative_accession,
                    "alternative_gene_id": pair.alternative_gene_id,
                    "alternative_gene_name": pair.alternative_gene_name,
                    "projection_status": evidence.status,
                    "projection_reason": evidence.reason,
                    "target_variant_key": evidence.target_variant_key,
                    "target_in_ortholog_gene": target_in_ortholog,
                    "target_consequence": "" if target_effect is None else target_effect.consequence,
                    "target_impact": "" if target_effect is None else target_effect.impact,
                    "effect_relation": _effect_relation(source_effect, target_effect),
                }
            )

        status_counts = Counter(item.status for item in pair_evidence)
        category = _category(
            set(baseline),
            set(alternative),
            pair_evidence,
            exact_source_keys,
            exact_target_keys,
        )
        if category not in BURDEN_CATEGORIES:
            raise LiftoverError(f"unknown burden category: {category}")
        explanation_rows.append(
            {
                "baseline_accession": pair.baseline_accession,
                "baseline_gene_id": pair.baseline_gene_id,
                "baseline_gene_name": pair.baseline_gene_name,
                "alternative_accession": pair.alternative_accession,
                "alternative_gene_id": pair.alternative_gene_id,
                "alternative_gene_name": pair.alternative_gene_name,
                "comparison_status": pair.comparison_status,
                "baseline_variant_count": len(baseline),
                "alternative_variant_count": len(alternative),
                "exact_shared_variant_count": len(exact_source_keys),
                "target_only_variant_count": len(set(alternative) - exact_target_keys),
                "accepted_match_count": status_counts["accepted_match"],
                "rejected_match_count": status_counts["rejected_match"],
                "target_reference_match_count": status_counts["target_reference_match"],
                "callable_no_candidate_count": status_counts["callable_no_candidate"],
                "target_not_callable_count": status_counts["target_not_callable"],
                "unresolved_projection_count": sum(
                    status_counts[status] for status in UNRESOLVED_STATUSES
                ),
                "burden_explanation": category,
            }
        )
    return equivalence_rows, explanation_rows


def connect_review_candidates(
    review_contexts: Sequence[Mapping[str, str]],
    projections: Sequence[ProjectionEvidence],
    equivalence_rows: Sequence[Mapping[str, object]] = (),
    explanation_rows: Sequence[Mapping[str, object]] = (),
) -> list[dict[str, object]]:
    """Attach one projection status to every inherited review-candidate context."""
    projection_index = _unique_index(
        projections,
        lambda item: (item.alternative_accession, item.source_variant_key),
        "directional projection",
    )
    equivalence_index = _unique_index(
        equivalence_rows,
        lambda row: (
            row["alternative_accession"],
            row["baseline_gene_id"],
            row["baseline_variant_key"],
        ),
        "review-candidate equivalence context",
    )
    explanation_index = _unique_index(
        explanation_rows,
        lambda row: (row["alternative_accession"], row["baseline_gene_id"]),
        "review-candidate burden explanation",
    )
    output: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    for context in review_contexts:
        source_key = context.get("baseline_variant_key", "")
        accession = context.get("alternative_accession", "")
        key = (source_key, accession)
        if not source_key or not accession:
            raise LiftoverError("review context lacks baseline_variant_key or alternative_accession")
        if key in seen:
            raise LiftoverError(f"duplicate review context: {key}")
        seen.add(key)
        evidence = projection_index.get((accession, source_key))
        if evidence is None:
            raise LiftoverError(f"review context lacks projection: {key}")
        gene_id = context.get("baseline_gene_id", "")
        equivalence = equivalence_index.get((accession, gene_id, source_key))
        explanation = explanation_index.get((accession, gene_id))
        output.append(
            {
                **context,
                "projection_status": evidence.status,
                "projection_reason": evidence.reason,
                "target_variant_key": evidence.target_variant_key,
                "exact_accepted_match": evidence.status == "accepted_match",
                "target_in_ortholog_gene": (
                    "" if equivalence is None else equivalence["target_in_ortholog_gene"]
                ),
                "target_consequence": (
                    "" if equivalence is None else equivalence["target_consequence"]
                ),
                "target_impact": "" if equivalence is None else equivalence["target_impact"],
                "effect_relation": "" if equivalence is None else equivalence["effect_relation"],
                "burden_explanation": (
                    "" if explanation is None else explanation["burden_explanation"]
                ),
            }
        )
    return output


def write_rows(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise LiftoverError(f"cannot write empty table without a declared schema: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise LiftoverError(f"missing Project 15 table: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _required_columns(rows: Sequence[Mapping[str, str]], required: set[str], label: str) -> None:
    if not rows:
        raise LiftoverError(f"Project 15 {label} table is empty")
    missing = required - set(rows[0])
    if missing:
        raise LiftoverError(f"Project 15 {label} table lacks columns: {sorted(missing)}")


def load_project15_gene_pairs(path: Path) -> tuple[list[GenePair], int]:
    rows = _read_csv(path)
    _required_columns(
        rows,
        {
            "baseline_accession",
            "baseline_gene_id",
            "baseline_gene_name",
            "alternative_accession",
            "alternative_gene_ids",
            "alternative_gene_names",
            "relationship",
            "comparison_status",
            "comparable",
        },
        "baseline ortholog comparison",
    )
    pairs: list[GenePair] = []
    for row in rows:
        comparable_text = row["comparable"].lower()
        if comparable_text not in {"true", "false"}:
            raise LiftoverError("Project 15 comparable column is not boolean")
        if comparable_text == "false":
            continue
        alternative_ids = [item for item in row["alternative_gene_ids"].split(";") if item]
        alternative_names = [item for item in row["alternative_gene_names"].split(";") if item]
        if row["relationship"] != "one_to_one" or len(alternative_ids) != 1:
            raise LiftoverError("a comparable Project 15 pair is not one-to-one")
        pairs.append(
            GenePair(
                baseline_accession=row["baseline_accession"],
                baseline_gene_id=row["baseline_gene_id"],
                baseline_gene_name=row["baseline_gene_name"],
                alternative_accession=row["alternative_accession"],
                alternative_gene_id=alternative_ids[0],
                alternative_gene_name=alternative_names[0] if alternative_names else "",
                comparison_status=row["comparison_status"],
                comparable=True,
            )
        )
    _unique_index(
        pairs,
        lambda pair: (pair.baseline_gene_id, pair.alternative_accession),
        "Project 15 comparable pair",
    )
    return pairs, len(rows)


def load_project15_gene_variants(
    effects_path: Path,
    gene_catalog_path: Path,
) -> tuple[list[GeneVariant], int, int]:
    effects = _read_csv(effects_path)
    genes = _read_csv(gene_catalog_path)
    _required_columns(
        effects,
        {
            "accession",
            "variant_key",
            "effect_order",
            "consequence",
            "impact",
            "gene_name",
            "gene_id",
            "feature_id",
            "warnings",
        },
        "annotation effects",
    )
    _required_columns(genes, {"accession", "gene_id"}, "gene catalog")
    valid_genes = {(row["accession"], row["gene_id"]) for row in genes}
    selected: dict[tuple[str, str, str], dict[str, str]] = {}
    unlinked = 0
    for effect in effects:
        if (effect["accession"], effect["gene_id"]) not in valid_genes:
            unlinked += 1
            continue
        if effect["impact"] not in IMPACT_ORDER:
            raise LiftoverError(f"unknown Project 15 impact: {effect['impact']}")
        key = (effect["accession"], effect["variant_key"], effect["gene_id"])
        rank = (
            IMPACT_ORDER[effect["impact"]],
            int(effect["effect_order"]),
            effect["feature_id"],
        )
        existing = selected.get(key)
        if existing is None:
            selected[key] = effect
            continue
        existing_rank = (
            IMPACT_ORDER[existing["impact"]],
            int(existing["effect_order"]),
            existing["feature_id"],
        )
        if rank < existing_rank:
            selected[key] = effect
    variants = [
        GeneVariant(
            accession=row["accession"],
            variant_key=row["variant_key"],
            gene_id=row["gene_id"],
            gene_name=row["gene_name"],
            consequence=row["consequence"],
            impact=row["impact"],
            warnings=row["warnings"],
        )
        for row in selected.values()
    ]
    return variants, len(effects), unlinked


def validate_project15_burdens(
    burden_path: Path,
    variants: Sequence[GeneVariant],
) -> int:
    burdens = _read_csv(burden_path)
    _required_columns(
        burdens,
        {"accession", "gene_id", "distinct_variant_count"},
        "gene consequence burden",
    )
    counts = Counter((item.accession, item.gene_id) for item in variants)
    seen: set[tuple[str, str]] = set()
    for row in burdens:
        key = (row["accession"], row["gene_id"])
        if key in seen:
            raise LiftoverError(f"duplicate Project 15 burden row: {key}")
        seen.add(key)
        if counts[key] != int(row["distinct_variant_count"]):
            raise LiftoverError(f"Project 15 burden disagrees with effects for {key}")
    if set(counts) - seen:
        raise LiftoverError("variant-gene assignments contain genes absent from the burden table")
    return len(burdens)


def load_project15_review_contexts(path: Path) -> list[dict[str, str]]:
    rows = _read_csv(path)
    _required_columns(
        rows,
        {"baseline_variant_key", "baseline_gene_id", "alternative_accession"},
        "review-candidate orthology",
    )
    keys = [(row["baseline_variant_key"], row["alternative_accession"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise LiftoverError("Project 15 review contexts contain duplicate variant/reference keys")
    return rows
