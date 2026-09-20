# Data provenance

## Sample lineage

The sample is public run `SRR13921545`. Project 16 does not download reads,
trim them, map them, or call new variants. It inherits the six reference-
relative Project 13 callsets created from the same complete retained paired-read
cohort used by Projects 12–15.

The inherited chain is:

```text
SRR13921545 paired reads
  → Project 11 trimming and mapping design
  → Project 12 MG1655-relative small variants
  → Project 13 six reference-relative callsets and callable regions
  → Project 14 MG1655 consequence annotations
  → Project 15 six-reference gene orthology and burden comparison
  → Project 16 exact cross-reference variant projection
```

## References

| Role | Accession | Label |
|---|---|---|
| Baseline | `GCF_000005845.2` | MG1655 |
| Alternative | `GCF_003697165.2` | DSM 30083 / ATCC 11775 |
| Alternative | `GCF_036503815.1` | 2017-02-2CC |
| Alternative | `GCF_005843885.1` | E4742 |
| Alternative | `GCF_000026225.1` | ATCC 35469 |
| Alternative | `GCF_002900365.1` | HT073016 |

Exact FASTA, accepted/evaluated/rejected VCF, and callable-BED paths and hashes
are listed in [`data/reference_manifest.csv`](../data/reference_manifest.csv).
The manifest deliberately uses relative paths so it does not publish a native
machine path.

## Inherited counts

| Artifact | Count |
|---|---:|
| References | 6 |
| Accepted variants | 357,581 |
| Annotation effects | 357,959 |
| Project 15 comparable gene pairs | 6,356 |
| MG1655 review candidates | 80 |
| Review-candidate reference contexts | 400 |

Project 16 independently reproduces these totals before creating an alignment.

## Integrity controls

- SHA-256 is checked before and after the workflow.
- Accepted and rejected VCF keys must be disjoint and together equal the
  evaluated VCF keys.
- VCF indexes must be present.
- Every PAF `cs` operation is validated against the source and target FASTA.
- Project 15 manifests and core tables are hashed and compared with the
  inherited records.
- Generated public artifacts have hashes in the authentic run manifest.

The original FASTAs, VCFs, BED files, annotations, and Project 15 tables are
read-only inputs. The workflow creates indexed reference copies under temporary
work storage when a native tool requires a local index.

## Public versus local data

The public repository retains compact results, controlled fixtures, normalized
VCFs, manifests, checksums, tests, documentation, and the dashboard. Native
whole-genome PAF files and copied references remain under ignored `local/`
storage in the portfolio run. A fresh CLI run retains its PAF files in the
chosen output directory for complete auditing.
