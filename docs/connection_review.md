# Milestone 1 — Connection and whole-genome alignment review

## Purpose

This milestone verifies that Project 16 is directly connected to Project 15
and that the six inherited reference coordinate systems can be compared at
base resolution. No variant was declared equivalent during this milestone.

## Project 15 reproduction

The review independently counted the inherited tables and VCF records rather
than trusting the previous summary alone.

| Metric | Expected | Observed |
|---|---:|---:|
| References | 6 | 6 |
| Accepted variants | 357,581 | 357,581 |
| Annotation effects | 357,959 | 357,959 |
| MG1655-to-alternative gene comparisons | 23,255 | 23,255 |
| Comparable one-to-one gene pairs | 6,356 | 6,356 |
| MG1655 review candidates | 80 | 80 |
| Review-candidate alternative contexts | 400 | 400 |

The six evaluated VCFs were also confirmed to be exact disjoint unions of
their accepted and rejected VCFs.

## Alignment method

Each alternative was aligned to MG1655 and MG1655 was independently aligned
back to that alternative:

```text
minimap2 -cx asm20 --cs=long --eqx --secondary=yes -N 5 -t 1
```

Every PAF record was required to contain parseable `cg` and long `cs` tags.
A primary block qualified when it was at least 10,000 bp, at least 85%
identical, and had MAPQ at least 20. A block was counted as reciprocal only
when an independently generated reverse-direction block covered at least half
of both corresponding intervals.

## Observed reciprocal coverage

| Alternative | MG1655 covered | Alternative covered | Interpretation |
|---|---:|---:|---|
| DSM 30083 / ATCC 11775 | 54.81% | 48.95% | Audit-only |
| 2017-02-2CC | 68.99% | 61.97% | Audit-only; closest to the 70% gate |
| E4742 | 34.05% | 30.94% | Audit-only |
| ATCC 35469 | 23.76% | 23.76% | Audit-only |
| HT073016 | 20.24% | 19.07% | Audit-only |

The ten alignments contain 1,669 PAF records. Of these, 472 passed the block
filters and 439 had a reciprocal counterpart. None of the five alternatives
met 70% coverage in both genomes.

This result is not an implementation failure. It shows why Project 15's
gene-level orthology comparison must not be silently upgraded to a claim of
whole-genome equivalence. Milestone 2 can still measure which individual
variants lie inside the reliable reciprocal blocks, while retaining all other
variants as nonalignable, ambiguous, nonreciprocal, or complex context.

## Reproducibility and integrity

- Python 3.12.14, Minimap2 2.31-r1302, BCFtools 1.24, Samtools 1.24, and
  Matplotlib 3.11.1 were used in an isolated WSL2 environment.
- SHA-256 checksums were recorded for all ten PAF files.
- Source checksums were calculated before and after alignment and were
  identical.
- PAF files and native logs remain under ignored `local/` storage.
- Compact validation, block, coverage, checksum, and manifest evidence is in
  `results/connection_review/`.

## Interpretation boundary

Reciprocal genome alignment establishes where exact coordinate translation
may be attempted. It does not by itself establish that two VCF records are the
same variant, that unmatched variants are biologically absent, or that an
alternative reference is suitable as the sample's primary reference.
