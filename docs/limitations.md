# Limitations

## Reference dependence

Every callset is relative to its own reference. A sample allele can be a variant
against MG1655 and simultaneously be the reference allele in another genome.
`target_reference_match` records this distinction; it does not mean a mutation
was gained or lost.

## Incomplete whole-genome comparability

None of the five alternatives passed the 70% reciprocal-coverage gate in both
directions. The closest alternative reached 68.99% coverage of MG1655 and
61.97% of itself. Exact matches are therefore local evidence within reliable
blocks, not proof that two references provide interchangeable whole-genome
coordinate systems.

The three least comparable references contribute very few Project 15 gene
pairs under the one-to-one orthology and callability rules. Their gene-level
percentages should not be generalized.

## Conservative unresolved states

The workflow does not force variants through:

- equally scoring repeat placements;
- alignment gaps or discontinuous REF spans;
- rearrangement boundaries;
- nonreciprocal mappings;
- missing chromosomes or plasmids;
- circular-origin-spanning alleles; or
- unsupported complex alleles.

These cases remain explicit unresolved outcomes. Unresolved is not equivalent
to biologically absent.

## Variant scope

The workflow handles normalized biallelic SNPs and short insertions/deletions.
It does not analyze structural variants, copy-number changes, large
rearrangements, multiallelic records, symbolic alleles, breakends, or graph
genome paths.

## Orthology and consequences

Project 15 one-to-one orthology is used only to provide a gene context. It does
not guarantee identical gene function. Predicted consequences can differ
because annotations, coding boundaries, strands, or reference alleles differ.
They are not experimental evidence of phenotype or gene damage.

## Interpretation

An accepted exact match means the same normalized allele is supported under the
fixed alignments, reference sequences, VCF partitions, and callability rules.
It does not prove:

- which lineage acquired the allele;
- that the variant is biologically causal;
- antimicrobial resistance;
- phenotype;
- sample identity; or
- clinical significance.

## Out of scope

This project performs no read remapping, variant recalling, pangenome graph
analysis, structural-variant calling, pathway analysis, resistance prediction,
phenotype inference, phylogeny, or clinical interpretation.
