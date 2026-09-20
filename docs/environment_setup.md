# Linux and WSL2 environment setup

## Requirements

The workflow is designed for a Linux command line. It was executed in Ubuntu
through WSL2 and can also run on a regular Linux host.

Pinned versions:

- Python 3.12
- Minimap2 2.31
- BCFtools 1.24
- Samtools 1.24
- Matplotlib 3.11.1
- pytest 9.0.2

## Create the environment

Install Micromamba, then run from the repository root:

```bash
micromamba create -f environment.yml
micromamba activate variant-liftover-explorer
python -m pip install --no-build-isolation --no-deps -e .
```

Verify the installation:

```bash
python --version
minimap2 --version
bcftools --version
samtools --version
variant-liftover-explorer --version
```

## Run the offline tests

```bash
python -m pytest -q -p no:cacheprovider
```

The tests use only committed controlled fixtures and compact authentic result
artifacts. Network access is not required.

## Run the workflow

Project 13 and Project 15 must be available as sibling directories, or the
manifest and input-directory arguments must be changed to their actual
locations.

```bash
variant-liftover-explorer \
  --reference-manifest data/reference_manifest.csv \
  --project15-comparison-dir ../bioinformatics-project-15-cross-reference-gene-orthology-consequence-explorer/results/authentic_comparison \
  --project15-orthology-dir ../bioinformatics-project-15-cross-reference-gene-orthology-consequence-explorer/results/orthology_analysis \
  --project15-annotation-dir ../bioinformatics-project-15-cross-reference-gene-orthology-consequence-explorer/local/authentic_annotation \
  --baseline-accession GCF_000005845.2 \
  --sample-id SRR13921545 \
  --alignment-preset asm20 \
  --minimum-alignment-length 10000 \
  --minimum-alignment-identity 0.85 \
  --minimum-alignment-mapq 20 \
  --ambiguity-score-fraction 0.95 \
  --minimum-comparable-coverage 0.70 \
  --threads 1 \
  --output-dir local/full_project16_run
```

The output directory must not exist. Ten PAF files and their native logs are
retained in a full run, so choose a location with sufficient storage.

## WSL2 notes

- Run native tools inside the same Linux environment; do not mix Windows and
  Linux executables in one workflow.
- Prefer a Linux filesystem for the fastest large-file processing. A mounted
  Windows path works but may be slower.
- The portfolio results use one thread for deterministic ordering. Increasing
  threads changes performance, not the evidence rules.
- Set `MPLBACKEND=Agg` on headless hosts.

## GitHub Actions

The repository workflow creates the pinned Micromamba environment on Ubuntu,
installs the package without downloading Python dependencies, verifies the
native tools, and runs the offline suite. Authentic source packages and local
PAFs are not required by CI.
