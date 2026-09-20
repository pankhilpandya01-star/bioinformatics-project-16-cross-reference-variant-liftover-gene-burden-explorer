# Milestone 2 — Liftover engine review

## Outcome

The controlled liftover demonstration passed. Twelve source variants were each
assigned exactly one final status, and no unresolved example was converted into
an exact-match claim.

| Status | Count | Controlled example |
|---|---:|---|
| `accepted_match` | 4 | Forward SNP, insertion, deletion, and reverse-strand SNP |
| `rejected_match` | 1 | Exact projected allele present only in the target rejected set |
| `target_reference_match` | 1 | Source alternate allele equals the target reference base |
| `callable_no_candidate` | 1 | Callable target locus without a matching evaluated candidate |
| `target_not_callable` | 1 | Projected target REF span is outside the callable intervals |
| `ambiguous_mapping` | 1 | Two equal-scoring repeat mappings |
| `nonreciprocal_mapping` | 1 | Forward mapping exists but no reverse mapping recovers the allele |
| `nonalignable` | 1 | No qualifying alignment covers the source REF span |
| `complex_alignment_context` | 1 | An insertion anchor coincides with a reference-alignment gap |

The complete row-level evidence is in
[`controlled_projections.csv`](../results/controlled_core/tables/controlled_projections.csv).

## Conservative decision path

For each source variant, the engine:

1. Confirms that the VCF REF allele exactly matches the source FASTA.
2. Requires one qualifying alignment to cover the complete REF span.
3. Rejects near-best competing alignments at or above 95% of the best score.
4. Translates every REF base through the long `cs` operations.
5. Reverse-complements alleles when the mapping strand is negative.
6. Withholds variants crossing gaps or unsupported complex contexts.
7. Normalizes the target allele with BCFtools against the target FASTA.
8. Projects the normalized target allele through the reverse alignment.
9. Requires the reverse-normalized key to equal the original source key.
10. Only then checks the target accepted set, rejected set, and callable BED.

The `target_reference_match` case uses complete reciprocal coordinate validation
because an equal REF and ALT is not a VCF variant and cannot be normalized as
one.

## Controlled allele checks

The accepted indel examples remained normalized as:

```text
tgt_plus:10:C:CAA
tgt_plus:13:AC:A
```

The reverse-strand source SNP `src_reverse:4:A:G` became
`tgt_reverse:17:T:C`, then returned to the original key through the reverse
alignment.

The published controlled VCF contains seven non-reference target alleles:
four accepted matches, one rejected match, one callable target without a
candidate, and one non-callable target. It is BGZF-compressed, CSI-indexed, and
readable with BCFtools 1.24.

## Atomic output behavior

The demonstration writes into a uniquely named staging directory. The final
output directory appears only after projection, normalization, VCF compression,
indexing, accounting, and version checks succeed. An injected-failure test
confirms that neither the final directory nor a staging directory remains after
a failed run.

## Verification

The complete offline suite contains 17 passing tests. Milestone 2 adds checks
for malformed PAF/CIGAR data, `cs` sequence disagreement, source REF mismatch,
duplicate compressed-VCF keys, reverse complementation, all final status
counts, reverse-strand and indel equivalence, compressed VCF accounting, output
portability, and atomic cleanup.

## Boundary

These controlled sequences verify coordinate and allele mechanics. They are not
biological findings. Milestone 3 will connect exact-equivalence results to
Project 15 orthologous genes, predicted consequences, burden comparisons, and
the 80 review candidates.
