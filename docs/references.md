# References

## Software and formats

1. Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*.
   2018;34(18):3094–3100. [Minimap2 repository and documentation](https://github.com/lh3/minimap2).
2. Minimap2 developers. [Assembly alignment, PAF, `asm20`, and liftover examples](https://github.com/lh3/minimap2/blob/master/cookbook.md).
3. Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and
   BCFtools. *GigaScience*. 2021;10(2):giab008.
   [BCFtools manual](https://samtools.github.io/bcftools/bcftools).
4. GA4GH Large Scale Genomics Work Stream.
   [HTS format specifications](https://github.com/samtools/hts-specs), including
   VCF, BED, BGZF, and CSI.
5. Matplotlib development team. [Matplotlib documentation](https://matplotlib.org/stable/).

## Data lineage

6. NCBI Sequence Read Archive. [`SRR13921545`](https://www.ncbi.nlm.nih.gov/sra/SRX10301019).
7. NCBI Assembly. [*Escherichia coli* K-12 MG1655 assembly `GCF_000005845.2`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/).

## Implementation choices

Minimap2's long `cs` representation is used because it carries matching bases,
substitutions, insertions, and deletions in one alignment tag. BCFtools performs
reference checking and left-normalization so repeat-equivalent indels are
compared in one consistent representation. VCF/BED and CSI validation follow
the maintained HTS specifications.

Accessed 2026-09-18.
