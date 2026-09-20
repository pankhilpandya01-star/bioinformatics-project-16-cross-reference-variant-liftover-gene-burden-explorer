# Milestone 5 — Quality review

## Decision

Milestone 5 passed. The expanded Python 3.12 suite contains 97 test cases,
including parameterized malformed-input cases. All 97 passed with no failures,
errors, or skips.

A separate read-only verification reconfirmed the authentic analysis totals,
all inherited Project 13 and Project 15 checksums, every recorded authentic
output hash, the dashboard hash, and the 20 compressed/indexed directional VCF
artifacts. No authentic result total changed during this review.

## Hardening changes

The review strengthened five boundaries:

1. PAF parsing now rejects non-integer fields, invalid sequence lengths,
   impossible match counts, MAPQ outside 0–255, missing or invalid primary/
   secondary tags, and negative alignment scores.
2. VCF parsing now requires one `#CHROM` header before data, rejects duplicate
   headers, invalid positions, symbolic or multiallelic alleles, and malformed
   rows while still supporting an empty callset with a valid header.
3. BED parsing now reports invalid numeric coordinates and empty contig names
   explicitly.
4. Batch projection rejects duplicate source keys, overlapping accepted and
   rejected evidence, and invalid or unsorted callable intervals before any
   native tool runs.
5. Consequence comparison remains restricted to accepted reciprocal exact
   matches. Coordinate hints from failed reciprocal checks cannot receive a
   target-gene consequence.

Circular-origin-spanning alleles receive an explicit instruction to split and
normalize before liftover. They are not forced through a linear alignment.
Native executables that cannot start are also converted into clear workflow
errors so staged outputs can be removed consistently.

## Test coverage

| Area | Cases reviewed |
|---|---|
| FASTA | duplicate identifiers, invalid bases, sequence before header, empty records |
| VCF | missing or duplicate headers, duplicate keys, invalid positions, equal alleles, symbolic and multiallelic records, empty callsets |
| BED | missing fields, non-integer, negative, empty, and unsorted intervals |
| PAF | malformed fields, count disagreement, invalid MAPQ and tags, CIGAR/`cs` disagreement, unknown contigs, sequence and length disagreement |
| Variant types | SNPs, insertions, deletions, repeat-left-normalized indels, allele swaps |
| Genome structure | forward and reverse strands, multiple contigs, missing plasmids, rearrangement gaps, circular-origin spans |
| Mapping evidence | ambiguous repeats, nonreciprocal mappings, nonalignable loci, non-callable positions, accepted and rejected candidates |
| Project 15 integration | duplicate and missing joins, consequence precedence, all six burden explanations, all 400 review contexts |
| Native tools | failed normalization, missing executables, pinned versions, compressed/indexed VCF validation |
| Atomic publication | controlled, integration, authentic, and quality-review failure cleanup |
| Authentic artifacts | 481,473 attempts, 154,865 baseline rows, 357,581 variants, 6,356 gene pairs, 400 review contexts, hashes and dashboard |
| Publication hygiene | portable commands and scans for machine paths, credentials, personal metadata, and authorship wording |

The machine-readable matrix is in
[`results/quality_review/test_matrix.csv`](../results/quality_review/test_matrix.csv).

## Controlled biological edge cases

The controlled suite still assigns exactly one expected outcome to each test
variant:

- forward SNP, insertion, deletion, and reverse-strand SNP: accepted exact;
- target allele equal to the alternative reference: target-reference match;
- matching weak target call: rejected match;
- callable target with no candidate: callable absent;
- uncovered target locus: not callable;
- equally scoring repeat placements: ambiguous mapping;
- missing source-contig alignment: nonalignable;
- one-way mapping: nonreciprocal;
- allele at an alignment gap: complex context.

A separate homopolymer insertion test confirms that BCFtools left-normalizes an
equivalent repeat-position allele before comparison.

## Tool verification

The live isolated Linux environment reports:

- Python 3.12.14
- Minimap2 2.31-r1302
- BCFtools 1.24
- Samtools 1.24
- Matplotlib 3.11.1
- pytest 9.0.2

These versions match the project environment.

## Atomic and accounting checks

Injected failures during the controlled workflow, Project 15 integration,
authentic analysis, and quality-review publication leave no final output and
no sibling staging directory. Existing output directories are rejected rather
than overwritten.

Authentic-result tests continue to verify:

- all 481,473 directional attempts have unique source/reference keys;
- all 30,973 MG1655 variants have exactly five alternative-reference rows;
- all 357,581 variants occur exactly once in the equivalence-cluster table;
- all 6,356 comparable gene pairs receive one burden explanation;
- all 80 review candidates receive five contexts;
- all projected and exact VCF keys are unique and agree with the manifest; and
- every recorded public-output checksum matches the current file.

## Publication boundary

Native PAF files, copied references, normalization work, and test reports remain
under ignored local storage. No repository was created or published. The public
CLI, continuous integration, final documentation, and publication checklist are
complete; repository creation remains separately authorized work.
