# Methods

## Study design

Project 16 is a direct extension of Project 15. It compares exact normalized
small variants across six reference coordinate systems and uses those results
to explain the earlier gene-level burden comparisons. No reads are remapped and
no variants are recalled.

MG1655 (`GCF_000005845.2`) is the baseline. The five alternatives are
`GCF_003697165.2`, `GCF_036503815.1`, `GCF_005843885.1`,
`GCF_000026225.1`, and `GCF_002900365.1`. All accepted, evaluated, and rejected
VCFs and callable-region BED files come from Project 13. Gene relationships,
annotation effects, burden tables, and review candidates come from Project 15.

## Input validation

Before alignment, the workflow checks:

- the six-row reference-manifest schema and unique accessions;
- exactly one baseline role matching the requested baseline accession;
- all reference, VCF, index, and callable-BED files and their SHA-256 values;
- accepted and rejected VCFs form an exact disjoint partition of each evaluated
  VCF;
- VCF sample, biallelic allele, coordinate, and duplicate-key rules;
- FASTA identifiers, bases, record lengths, and PAF sequence identities; and
- exact reproduction of Project 15's inherited totals and file hashes.

Source hashes are calculated again after analysis. Any change fails the run.

## Whole-genome alignment

MG1655 is aligned to each alternative and each alternative is independently
aligned back to MG1655, producing ten directional PAF files:

```text
minimap2 -cx asm20 --cs=long --eqx --secondary=yes -N 5 -t 1 \
  TARGET_FASTA SOURCE_FASTA > DIRECTION.paf
```

The `asm20` preset is intended for divergent assembly-to-assembly comparison.
The long `cs` tag carries aligned bases, substitutions, insertions, and
deletions, while `--eqx` makes match and mismatch CIGAR operations explicit.

A block is eligible for projection when it is:

- a primary alignment;
- at least 10,000 aligned bases;
- at least 85% identical; and
- MAPQ 20 or greater.

A competing block is considered ambiguous when its chain score reaches at
least 95% of the best block's score. Secondary alignments are retained so this
ambiguity can be measured rather than hidden.

## Base-level projection

Each accepted source variant is handled independently:

1. Confirm the VCF REF allele exactly matches the source FASTA.
2. Select one eligible, non-ambiguous source-to-target block that contains the
   complete source REF span.
3. Translate every REF base through the long `cs` operations.
4. Reject gaps, discontinuous target spans, unsupported complex alleles, and
   unsafe indel anchors as complex alignment contexts.
5. Reverse-complement alleles for reverse-strand blocks.
6. Re-express the allele against the target reference sequence.
7. Left-align and normalize with `bcftools norm -f`.
8. Project the normalized target allele back through the independently produced
   reverse alignment.
9. Normalize the returned source allele and require its complete key—contig,
   position, REF, and ALT—to equal the original source key.

An exact match is accepted only after both coordinate and normalized-allele
reciprocity succeed.

## Target evidence classification

Successful reciprocal projections are compared with the target callsets:

- `accepted_match`: exact normalized key in the target accepted VCF;
- `rejected_match`: exact normalized key in the target rejected VCF;
- `target_reference_match`: the sample allele equals the target reference;
- `callable_no_candidate`: target REF span is callable but no candidate exists;
- `target_not_callable`: at least one required target REF base is not callable.

Mappings that cannot safely reach evidence classification remain
`ambiguous_mapping`, `nonreciprocal_mapping`, `nonalignable`, or
`complex_alignment_context`.

## Equivalence clusters

Accepted directional exact matches are treated as graph edges connecting
`(accession, normalized_variant_key)` nodes. Connected components form
equivalence clusters. Singleton variants remain explicit one-member clusters;
they are not discarded.

## Project 15 integration

Only Project 15 gene pairs already classified as comparable one-to-one
orthologs are interpreted. The workflow selects each gene's highest-impact
effect for a `(variant, gene)` assignment, joins the exact projection, and
assigns one explanation:

1. `unresolved_liftover` when any relevant baseline projection is ambiguous,
   nonreciprocal, nonalignable, or complex;
2. `same_exact_variants` when the complete baseline and alternative gene sets
   are connected by accepted exact matches;
3. `partially_shared_variants` when at least one, but not all, variants match;
4. `reference_allele_difference` when no exact variant matches but at least one
   sample allele is the alternative reference allele;
5. `callability_or_filter_difference` for rejected or non-callable evidence;
6. `different_variants_in_same_gene` otherwise.

Uncertainty takes precedence so a partial match cannot conceal unresolved
evidence. Target-gene consequences are compared only for accepted exact matches
that fall inside the paired ortholog.

## Output and integrity

The public command writes to a uniquely named sibling staging directory. It
publishes the complete output with one rename only after alignment, projection,
normalization, integration, checksum, VCF, and accounting checks pass. Any
exception removes the staging tree.

Generated compressed VCFs use BGZF and CSI indexes. Every table has an exact
row-count invariant recorded in the manifests. Portable command templates omit
native paths.
