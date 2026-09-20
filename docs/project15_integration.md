# Milestone 3 — Project 15 integration review

## Outcome

The integration layer is complete. It connects exact projection evidence to
Project 15's one-to-one orthologs, selected variant–gene consequences, burden
tables, and review-candidate contexts.

The controlled demonstration contains six comparable gene pairs, with one pair
assigned to each allowed explanation:

| Explanation | Meaning in the controlled example |
|---|---|
| `same_exact_variants` | Every variant assigned to either gene is paired by an accepted exact projection |
| `partially_shared_variants` | At least one exact variant is shared, but one or both genes contain additional variants |
| `reference_allele_difference` | No exact variant is shared and the sample allele equals the alternative reference |
| `callability_or_filter_difference` | No exact variant is shared and the target evidence was rejected or not callable |
| `different_variants_in_same_gene` | The orthologs contain variant evidence, but no exact projected allele is shared |
| `unresolved_liftover` | At least one baseline variant has ambiguous, nonreciprocal, nonalignable, or complex mapping |

The row-level examples are in
[`project15_burden_explanation.csv`](../results/controlled_integration/tables/project15_burden_explanation.csv).

## Deterministic precedence

Some gene pairs can contain more than one kind of evidence. The workflow uses
this conservative precedence:

1. `unresolved_liftover`
2. `same_exact_variants`
3. `partially_shared_variants`
4. `reference_allele_difference`
5. `callability_or_filter_difference`
6. `different_variants_in_same_gene`

Unresolved evidence therefore cannot be hidden by one successfully projected
variant. A pair is called `same_exact_variants` only when the accepted projected
source keys and target keys exhaust both gene-level variant sets.

## Consequence connection

Project 15 may contain overlapping effects for one variant and gene. The bridge
reproduces Project 15's deterministic selection rule:

1. `HIGH`, `MODERATE`, `LOW`, then `MODIFIER` impact;
2. original effect order;
3. feature identifier.

Each projected variant can then receive one of four effect relationships:

- `same_consequence_and_impact`
- `same_impact_different_consequence`
- `different_impact`
- `target_gene_effect_missing`

An exact DNA allele does not imply an identical predicted consequence. The
controlled table includes examples where the exact allele has the same effect,
the same impact but a different consequence term, and a different impact.

## Review-candidate connection

Each inherited review context retains its baseline variant and gene fields and
gains:

- the exact projection status and target key;
- whether the target variant affects the one-to-one ortholog;
- the target consequence and impact when present;
- the effect relationship; and
- the gene-pair burden explanation.

The controlled output is
[`review_candidate_liftover.csv`](../results/controlled_integration/tables/review_candidate_liftover.csv).
The authentic output will contain five alternative-reference rows for each of
the 80 Project 14 review candidates.

## Actual Project 15 bridge audit

The read-only bridge audit reproduced:

| Inherited item | Count |
|---|---:|
| Baseline-to-alternative gene rows | 23,255 |
| Comparable one-to-one gene pairs | 6,356 |
| Annotation effects | 357,959 |
| Distinct variant–gene assignments | 327,542 |
| Effects without an exact catalog gene key | 30,382 |
| Gene-burden rows | 28,859 |
| Review-candidate contexts | 400 |
| Distinct review candidates | 80 |

Every variant–gene count agrees with the Project 15 burden table. The key input
hashes also match the earlier connection review, and all source files remained
unchanged.

## Verification and boundary

The complete suite now contains 30 passing tests. Milestone 3 adds category
precedence, duplicate and missing projection checks, burden-count agreement,
review-context completeness, effect-relation joins, noncomparable-pair handling,
atomic cleanup, exact inherited counts, checksum continuity, and output
portability.

The controlled categories demonstrate logic, not biological findings. The
actual Project 13 variants have not yet been projected or used to revise any
Project 15 conclusion. That authentic analysis is Milestone 4.
