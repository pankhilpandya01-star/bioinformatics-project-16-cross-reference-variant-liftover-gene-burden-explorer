# Milestone 4 — Authentic cross-reference analysis

## Purpose

This milestone asks whether MG1655-relative variants describe the same
normalized DNA alleles when they are translated into each of the five
alternative reference coordinate systems. It then uses those exact matches to
explain the gene-burden comparisons from Project 15.

An unmatched allele is not automatically absent from the sample. The workflow
keeps separate outcomes for different target reference alleles, rejected
candidates, non-callable loci, nonreciprocal mappings, and regions that cannot
be translated conservatively.

## Complete accounting

The authentic run processed all ten reference directions and reconciled every
required record:

| Metric | Result |
|---|---:|
| Directional projection attempts | 481,473 |
| MG1655-variant-by-alternative rows | 154,865 |
| Accepted directional exact-match edges | 40,418 |
| Variants represented in equivalence clusters | 357,581 |
| Equivalence clusters | 337,372 |
| Multi-reference clusters | 14,409 |
| Project 15 comparable gene pairs | 6,356 |
| Variant–ortholog evidence rows | 44,097 |
| Review-candidate contexts | 400 |

Each directional projection received exactly one status. All ten normalized
projected VCFs and all ten exact-equivalence VCFs are compressed, indexed,
reference-consistent, and free of duplicate variant keys. Source checksums were
unchanged after the run.

## Exact matches from MG1655

The 30,973 MG1655-relative variants were attempted against each alternative,
giving 154,865 comparison rows.

| Alternative reference | Accepted exact matches | MG1655 variants matched |
|---|---:|---:|
| 2017-02-2CC | 10,303 | 33.26% |
| DSM 30083 / ATCC 11775 | 7,350 | 23.73% |
| E4742 | 1,236 | 3.99% |
| HT073016 | 770 | 2.49% |
| ATCC 35469 | 550 | 1.78% |

Across all five alternatives, 14,409 distinct MG1655 variants had at least one
accepted exact match. Most were reference-specific: 9,368 matched exactly one
alternative, 4,325 matched two, 675 matched three, 39 matched four, and only two
matched all five.

The complete MG1655 outcome distribution was:

| Outcome | Rows | Share |
|---|---:|---:|
| `nonalignable` | 84,261 | 54.41% |
| `target_reference_match` | 30,699 | 19.82% |
| `accepted_match` | 20,209 | 13.05% |
| `rejected_match` | 10,523 | 6.79% |
| `nonreciprocal_mapping` | 4,313 | 2.79% |
| `target_not_callable` | 3,662 | 2.36% |
| `complex_alignment_context` | 1,037 | 0.67% |
| `callable_no_candidate` | 161 | 0.10% |

The large nonalignable fraction is consistent with Milestone 1: none of the
five alternatives passed the 70% reciprocal whole-genome coverage gate. The
30,699 `target_reference_match` rows are also important. In those cases, the
sample allele called as a difference from MG1655 is already the reference
allele in the alternative genome. That is a reference-allele difference, not
evidence that the allele disappeared.

## What changed in the Project 15 gene conclusions

Each of Project 15's 6,356 comparable one-to-one gene pairs received one
explanation:

| Explanation | Gene pairs | Share |
|---|---:|---:|
| `partially_shared_variants` | 2,303 | 36.23% |
| `different_variants_in_same_gene` | 2,075 | 32.65% |
| `unresolved_liftover` | 1,428 | 22.47% |
| `reference_allele_difference` | 324 | 5.10% |
| `same_exact_variants` | 225 | 3.54% |
| `callability_or_filter_difference` | 1 | 0.02% |

The `same_exact_variants` count needs a careful denominator. Of those 225 gene
pairs, 201 contain zero variants in both genes. Only 24 contain one or more
shared exact variants. The strongest conclusion is therefore that Project 15's
matching raw burdens frequently did **not** arise from completely identical
variant sets.

Only DSM 30083 / ATCC 11775 and 2017-02-2CC contribute large comparable-gene
denominators: 3,206 and 3,137 pairs, respectively. E4742, ATCC 35469, and
HT073016 contribute only 9, 1, and 3 comparable pairs. Those three references
cannot support broad gene-level stability claims under the fixed callability
and one-to-one-orthology rules.

## Consequence stability among accepted gene matches

There were 14,280 accepted exact projections whose target allele also occurred
in the paired one-to-one ortholog:

| Predicted-effect relation | Rows |
|---|---:|
| Same consequence and impact | 14,201 |
| Different consequence, same impact | 2 |
| Different impact | 77 |

Thus, when an exact allele was successfully placed inside the paired ortholog,
the predicted consequence and impact were usually stable. This statement does
not apply to unresolved projections or to accepted projections outside the
paired ortholog. Their target consequences are intentionally left blank.

## Project 14 review candidates

The 80 MG1655 review candidates generated exactly 400 alternative-reference
contexts. Forty-three contexts were accepted exact matches, representing 32
distinct candidates. Twenty-six of those 43 matches also landed in the paired
one-to-one ortholog. Most contexts were nonalignable (274 of 400), so absence of
an exact match must not be interpreted as evidence that a review candidate is
biologically absent.

## Interpretation

Project 15's broad gene-level signal is only partly stable at exact-variant
resolution. The closest alternatives recover substantial exact subsets, but
no alternative provides a complete equivalent coordinate system. Some burden
differences reflect different reference alleles; others reflect different
variants in the same ortholog, filtering or callability, or unresolved genome
alignment.

MG1655 remains the primary interpretation reference. All alternative-reference
results remain sensitivity-analysis evidence because none passed the earlier
whole-genome mapping and reciprocal-coverage eligibility gates.

## Reproducibility

- Python 3.12.14, Minimap2 2.31-r1302, BCFtools 1.24, Samtools 1.24, and
  Matplotlib 3.11.1 were used in the isolated Linux environment.
- Projection and integration outputs were published atomically after their
  cross-artifact checks passed.
- The manifest records parameters, versions, inherited hashes, output hashes,
  status counts, VCF counts, and the unchanged-input check.
- Native reference copies and normalization work stayed under ignored local
  storage and were removed after successful publication.

The dashboard is in
[`results/authentic_analysis/dashboard/variant_liftover_dashboard.png`](../results/authentic_analysis/dashboard/variant_liftover_dashboard.png).
