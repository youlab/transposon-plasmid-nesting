# Numbered bioinformatics pipeline

All analysis scripts are in `scripts/`. Numbers describe the analysis steps; 04 and 09 are external-tool jobs under `slurm/`. Scripts 17/18 and 22 are supporting functions used by the subsequent numbered entries.

## Run from processed data

```bash
python code/pipeline/scripts/24_run_backbone_analysis.py
```

This default route reads the compressed CDS table directly, creates before/after counts and calculates statistics. It does not need GBFF, raw FASTA, PlasAnn installation or reference databases, and it does not generate figures. It requires Python and SciPy. Install the full provided environment only if you also need upstream annotation tools.

Scripts resolve the package root from their own location. The commands below are shown from that root; launching a script by its absolute filename from another directory also works. Shared configuration is `config/analysis.json`; optional environment overrides are in `config/environment.example`.

## Script reference

| Step | Script | Purpose | Immediate inputs | Outputs | Example command |
|---|---|---|---|---|---|
| 01 | 01_extract_accessions.py | Select transposase-positive accessions | data/01_transposase_counts.csv | 01_transposase_positive_accessions.txt | python code/pipeline/scripts/01_extract_accessions.py |
| 02 | 02_download_gbff.py | Download GenBank records | Selected accession list; NCBI email | work/raw/gbff/*.gbff | sbatch code/pipeline/slurm/02_download_gbff.sbatch |
| 03 | 03_gbff_to_fna.py | Convert records and deduplicate exact sequences | work/raw/gbff/*.gbff | 03_plasmid_sequence_hashes.csv; deduplicated FNA | python code/pipeline/scripts/03_gbff_to_fna.py |
| 04 | slurm/04_run_isescan.sbatch | Annotate insertion sequences | Deduplicated plasmid FNA | 04_isescan/fna/*.fna.tsv | sbatch code/pipeline/slurm/04_run_isescan.sbatch |
| 05 | 05_compare_transposase_counts.py | Compare annotation counts with ISEScan counts | Sequence map; transposase-count table; ISEScan TSVs | 05_transposase_count_comparison.csv; summary | sbatch code/pipeline/slurm/05_count_qc.sbatch |
| 06 | 06_discordant_breakdown.py | Describe annotation-positive / ISEScan-negative records | Sequence map; GBFF; ISEScan TSVs; source count table | 06_discordant_cds.csv; summary | python code/pipeline/scripts/06_discordant_breakdown.py |
| 07 | 07_extract_cargo.py | Define relaxed or strict cargo/background | GBFF; ISEScan TSVs; sequence map | Cargo/background FASTAs; membership and span tables | sbatch code/pipeline/slurm/07_extract_cargo.sbatch |
| 08 | 08_dedup_proteins.py | Deduplicate protein search queries | 07_cargo_relaxed.faa; 07_background_relaxed.faa | 08_*_relaxed_unique.faa and copy maps | python code/pipeline/scripts/08_dedup_proteins.py |
| 09 | slurm/09_functional_annotation.sbatch | Run original functional-database searches | Deduplicated proteins; CARD/VFDB/BacMet/TADB references | 09_*_relaxed_vs_*.tsv; AMRFinder TSVs | sbatch code/pipeline/slurm/09_functional_annotation.sbatch |
| 10 | 10_function_inventory.py | Summarize functional cargo hits | Protein maps; original DIAMOND hits; cargo table | 10_cargo_function_inventory.csv; top genes; summary | python code/pipeline/scripts/10_function_inventory.py 90 60 |
| 11 | 11_enrichment.py | Original relaxed enrichment from search outputs | Relaxed protein maps and DIAMOND hits | 11_enrichment_relaxed.csv; optional 80/70 sensitivity | python code/pipeline/scripts/11_enrichment.py 90 60 |
| 12 | 12_amrf_enrichment.py | Original AMRFinderPlus sensitivity analysis | AMRFinder TSVs; relaxed protein maps and cargo FASTA | 12_amrfinder_enrichment_summary.txt | python code/pipeline/scripts/12_amrf_enrichment.py |
| 13 | 13_strict_analysis.py | Strict functional-search sensitivity analysis | Strict cargo/background FASTAs; functional databases | 13_enrichment_strict.csv; strict search outputs | sbatch code/pipeline/slurm/13_strict_analysis.sbatch |
| 14 | 14_strict_counts.py | Summarize current strict-set counts | Strict protein maps, cargo table and search outputs | 14_strict_counts.txt | python code/pipeline/scripts/14_strict_counts.py |
| 15 | 15_is_sequence_jumps.py | Identify identical IS copies in distinct outer-flank contexts | ISEScan TSVs; deduplicated FNA; sequence map | 15_is_insertion_contexts.csv; examples; summary | python code/pipeline/scripts/15_is_sequence_jumps.py |
| 16 | 16_final_figures.py | Original count-comparison and enrichment drawing | 05 count comparison; 11 relaxed enrichment table | Figures under results/figures | python code/pipeline/scripts/16_final_figures.py |
| 17 | 17_plasann_cds_helpers.py | Shared CDS/reference functions; prepare reference BLAST index | GBFF and exact PlasAnn Database.csv | Checked reference BLAST cache; helper functions | python code/pipeline/scripts/17_plasann_cds_helpers.py --prepare-db |
| 18 | 18_plasann_category_helpers.py | Normalize categories and check CDS membership sources | PlasAnn outputs; frozen cargo/source files | Category checks and helper functions | Imported by scripts 19/20; see --help for its optional summary command |
| 19 | 19_annotate_plasann_cds.py | Annotate existing CDSs with the agreed 80/90/90 rule | GBFF; sequence map; frozen PlasAnn reference | results/19_plasann_coverage_annotation/plasmids/<accession>/ | sbatch code/pipeline/slurm/19_annotate_plasann_cds.sbatch |
| 20 | 20_filter_and_count_cds.py | Apply the backbone rule and recount before/after CDSs | data/19_cds_annotation_membership.csv.gz; analysis.json | 20_cds_counts_before_filtering.csv; 20_cds_counts_after_filtering.csv; 20_backbone_filter_counts.csv | python code/pipeline/scripts/20_filter_and_count_cds.py |
| 21 | 21_calculate_backbone_enrichment.py | Calculate Fisher and within-mode BH for both gene sets | The two 20_* hit-count tables; analysis.json | 11_enrichment_relaxed.csv; 13_enrichment_strict.csv; 21_backbone_enrichment.csv | python code/pipeline/scripts/21_calculate_backbone_enrichment.py |
| 22 | 22_backbone_plot_helpers.py | Shared Extended Figure 1 drawing functions | Validated fractions passed by script 23 | PDF/SVG/PNG only when plotting is explicitly requested | Imported by script 23 |
| 23 | 23_plot_backbone_enrichment.py | Draw the configured backbone-filtered comparison | 20_cds_counts_after_filtering.csv; 21_backbone_enrichment.csv | Relaxed article Extended Figure 1 panel; strict SI sensitivity panel | python code/extended_fig1.py |
| 24 | 24_run_backbone_analysis.py | Run scripts 20 and 21 without plotting by default | Bundled compressed CDS input; analysis.json | Count and statistical tables under results/ | python code/pipeline/scripts/24_run_backbone_analysis.py |

## Paths and outputs

- `data/` contains supplied processed inputs and current result tables.
- `work/raw/` holds downloaded GBFFs, derived sequence FASTAs and external database snapshots for optional upstream runs.
- `work/processed/` holds generated upstream sequence maps, ISEScan calls, protein maps and search results. Scripts prefer a generated upstream input when available and otherwise use a compatible bundled table.
- `results/` receives the bundled-data recount and statistical outputs.
- `results/figures/` receives images only when a plotting entry is explicitly run.

The raw tables produced by script 21 are the before-filtering four-category outputs: relaxed `11_enrichment_relaxed.csv` supplies article Figure 1d, and strict `13_enrichment_strict.csv` is an SI sensitivity result where cited. Relaxed rows in `21_backbone_enrichment.csv` supply article Extended Figure 1; strict rows are SI sensitivity results where cited. This WithTA package uses four categories for both test families. Dropping TADB as a comparison never changes CDS denominators.

## Full upstream processing

Run steps 01–15 in dependency order, using the matching Slurm jobs for downloads and external-tool searches. Set NCBI_EMAIL for step 02. Restore the archived ISEScan TSVs into `work/processed/04_isescan/` if reusing the supplied IS calls rather than rerunning ISEScan. Cargo/background FASTAs are generated by step 07, not downloaded. The full sequence route requires the large raw files; the bundled-data route above does not.

For a fresh PlasAnn annotation, acquire the full Database.csv into `work/raw/PlasAnn_database/`, verify the frozen SHA256, prepare its BLAST cache with step 17, then run step 19. The official full-table endpoint is:

```text
https://plasann-production.up.railway.app/api/download/Database.csv
```

The archived reference is documented at https://zenodo.org/records/19446428. The frozen CSV has 213,884 CDSs / 11,942 accessions and SHA256 `72b005cf392f6b22f969a1218a2f71c8ad4dfb399efd1327ce9d8c1a695ddd57`. A current website download is acceptable only if this SHA matches. The archive metadata described a v2 database; this package does not call it paper Database v4.1.5.

To rejoin newly produced annotations with the frozen cargo and original search files, script 20 offers `--from-annotation`. This optional path needs the original FASTAs, query-copy maps, DIAMOND TSVs and current annotation directory. It preserves the original accession/coordinate/protein-sequence checks. Normal users can start with the included compressed CDS membership table instead.

## Statistical and annotation rules

The cargo span limit remains 25,000 bp (script 07 default, recorded as `cargo_max_span_bp` in config). This is distinct from PlasAnn's maximum of ten candidate reference targets and from the optional `--limit` number of plasmids used in a trial run.

- Cargo and background receive the same two-category PlasAnn backbone exclusion. Other/unmatched labels remain in the denominators.
- PlasAnn classification uses identity >=80%, query coverage >=90%, subject coverage >=90%, no lower-identity fallback.
- Original functional hits remain CARD/VFDB/BacMet/TADB calls at identity >=90%, query coverage >=60%.
- Counts represent CDS instances, not unique protein sequences. Protein search deduplication is reversed through copy/member mappings.
- Original comparisons use four Fisher tests per cargo mode; backbone-filtered tests use the categories specified by this package's config. Pooled-CDS Fisher is two-sided; BH is performed separately by gene set and cargo mode.
- Scientific probability strings and log10 P/q preserve extremely small results; no true P=0 is inferred from floating-point underflow.
- There is no within-plasmid dependence correction.

## Slurm and software

Slurm jobs locate this package from their own filenames and use the active Python environment or PYTHON_BIN. No personal account, partition or absolute HPC checkout is embedded. Request your site's partition/account on the `sbatch` command if required. The bundled-data job 24 does not submit annotation or draw figures. The step-19 array preserves the original eight-way sharding with 40 CPUs / 64 GB per task.

`environment.yml` lists the analysis, plotting and PlasAnn tools without platform build hashes. `environment-amrfinder.yml` is the optional AMRFinderPlus environment. `config/database_manifest.tsv` records the exact original functional reference snapshots. Scientific category/threshold choices are in `config/analysis.json`; credentials are excluded.

## Checks

```bash
python code/pipeline/scripts/20_filter_and_count_cds.py --self-test
python code/pipeline/scripts/21_calculate_backbone_enrichment.py --self-test
```

The self-tests cover the backbone rule, counts and known Fisher/BH cases. A full upstream sequence/search rerun requires the omitted large raw inputs.

`config/input_manifest.tsv` lists the size and SHA256 of each supplied data file, including the compressed CDS input. It describes the distributed inputs; generated work/results are excluded by `.gitignore`.
