# Bioinformatics code and data — WithTA

Code and processed data supporting the original Figure 1 and the additional PlasAnn backbone-filtered comparison. This package runs independently of the other TA version.

## Contents

- `code/fig1a.py`, `fig1b.py`, `fig1d.py`: original figure plotting entries.
- `code/extended_fig1.py`: relaxed backbone-filtered article Extended Figure 1 plus strict sensitivity plotting entry.
- `code/pipeline/scripts/`: numbered analysis scripts; see the pipeline README for each purpose, input, output and command.
- `code/pipeline/slurm/`: matching portable Slurm submission scripts.
- `code/pipeline/config/`: analysis categories, frozen database information and optional environment overrides.
- `data/`: the compressed CDS-level input and supporting tables; see `data/README.md`.

Generated images, logs, temporary files and large raw sequence/database inputs are not included. Commands create their output directories automatically. Packaged paths are resolved relative to the code location, so no personal server paths or path edits are needed for the bundled-data workflow.

## Install

```bash
conda env create -f code/pipeline/environment.yml
conda activate nee-bioinformatics
```

The optional AMRFinderPlus analysis uses `code/pipeline/environment-amrfinder.yml`. Environment files omit platform-specific build hashes. Full upstream jobs additionally require the external sequences and exact reference snapshots described below.

## Recalculate the supplied analyses

Run from the package root:

```bash
python code/pipeline/scripts/24_run_backbone_analysis.py
```

This streams `data/19_cds_annotation_membership.csv.gz` without extracting it, reproduces counts before and after backbone exclusion, and recalculates enrichment statistics. It writes CSVs to `results/` and does not create figures.

The article uses relaxed cargo definitions: Figure 1d shows the original four-category comparison, and Extended Figure 1 shows the four-category PlasAnn backbone-filtered comparison. Strict results are SI sensitivity analyses where cited. Both comparisons use two-sided pooled-CDS Fisher tests with BH correction across four tests within each cargo mode.

TA-positive CDSs are not removed as an additional filter. Both groups use the same independent backbone rule. Unmatched, hypothetical and unclassified CDSs remain in the denominators unless assigned to a specified backbone category.

## Draw figures when needed

```bash
python code/fig1a.py
python code/fig1b.py
python code/fig1d.py
python code/extended_fig1.py
```

These separate entries read the supplied tables and export figures under `results/figures/`; they are not run by the statistical workflow. Plotting uses Arial by font name and reports a fallback if Arial is unavailable. No font path needs editing. The plotting entry exports relaxed and strict panels; use the relaxed filtered panel for article Extended Figure 1 and the strict panel only as SI sensitivity output where cited.

## Data flow

```text
Compressed CDS annotations, memberships and original functional-hit flags
  -> 20_filter_and_count_cds.py
     -> counts before filtering
     -> counts after filtering
     -> total / excluded / retained CDS counts
  -> 21_calculate_backbone_enrichment.py
     -> 11_enrichment_relaxed.csv and 13_enrichment_strict.csv (original four-category comparisons)
     -> 21_backbone_enrichment.csv (4-category backbone-filtered comparisons)
  -> Article Figure 1d (raw relaxed) / Extended Figure 1 (filtered relaxed); strict is SI sensitivity only where cited
```

## Full upstream processing

The retained upstream scripts download GBFF records, annotate ISs, define cargo/background, generate FASTAs, and run functional searches. Large GBFF/FASTA files, search results and reference databases are excluded; they are generated or acquired into `work/raw/` and `work/processed/`. NCBI downloads require an email address and optional API key. Never put credentials in a public package.

PlasAnn 1.1.6 (commit `35199e00d90462b9cf84e9f924ef098edd5da122`) uses the frozen 2026-10-19 official database snapshot: 213,884 CDSs / 11,942 accessions. The full-table source and exact SHA256 are in `config/analysis.json` under `code/pipeline/`. The public archive is https://zenodo.org/records/19446428. This snapshot is not labelled as the paper's Database v4.1.5. Verify the recorded SHA256 before using a downloaded database; the live website may change.

Annotation uses identity >=80%, query coverage >=90% and subject coverage >=90%, without a lower-identity fallback. Original functional-hit flags use DIAMOND identity >=90% and query coverage >=60%. Backbone exclusions are Conjugation and Plasmid Maintenance, Replication and Regulation, including the documented category aliases. Statistics use pooled-CDS two-sided Fisher tests and do not adjust for within-plasmid dependence.
