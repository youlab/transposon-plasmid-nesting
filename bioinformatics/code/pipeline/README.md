# NEE sequence-level IS and cargo-function analysis

This folder contains the code used for sequence-level annotation of insertion sequence (IS) elements and for analysis of the functions enriched in putative IS-associated cargo. The code files are listed below in execution order together with their principal outputs.

Input and output directories can be configured with the environment variables documented in `config/environment.example`. Large sequence files and intermediate search results are not included in this repository.

| Step | Code | Analysis | Main output files |
|---:|---|---|---|
| 01 | `scripts/01_extract_accessions.py` | Extract plasmid accessions with transposase-related annotations | `1_Hit_plasmids.txt` |
| 02 | `scripts/02_download_gbff.py`; `slurm/02_download_gbff.sbatch` | Download GenBank records | Not included because the downloaded GBFF records are large |
| 03 | `scripts/03_gbff_to_fna.py` | Convert plasmid records to FASTA and remove exact sequence duplicates | `3_Seqhash_map.csv`; `3_Dedup_summary.txt` (converted FNA files are not included because of their size) |
| 04 | `slurm/04_run_isescan.sbatch` | Annotate IS elements with ISEScan | `4_ISEScan/fna/<md5>.fna.tsv` |
| 05 | `scripts/05_compare_transposase_counts.py`; `slurm/05_count_qc.sbatch` | Compare annotation-based transposase counts with ISEScan-detected IS counts | `4b_Transposase_count_comparison.csv`; `4b_Compare_summary.txt`; `4b_Transposase_count_comparison.pdf`; `4b_Transposase_count_comparison.png` |
| 06 | `scripts/06_discordant_breakdown.py` | Summarize records with transposase annotations but no ISEScan-detected IS | `4c_Discordant_cds.csv`; `4c_Summary.txt` |
| 07 | `scripts/07_extract_cargo.py`; `slurm/07_extract_cargo.sbatch` | Extract putative cargo and within-plasmid background proteins using relaxed and strict definitions | `5_Cargo.faa`; `5_Background.faa`; `5_Cargo_table.csv`; `5_Span_distribution.csv`; `5_Summary.txt`; `5strict_Cargo.faa`; `5strict_Background.faa`; `5strict_Cargo_table.csv`; `5strict_Span_distribution.csv`; `5strict_Summary.txt` |
| 08 | `scripts/08_dedup_proteins.py` | Deduplicate cargo and background protein sequences | `6_Cargo_unique.faa`; `6_Cargo_unique_map.csv`; `6_Background_unique.faa`; `6_Background_unique_map.csv` |
| 09 | `slurm/09_functional_annotation.sbatch` | Search cargo and background proteins against CARD, VFDB, BacMet and TADB, and run AMRFinderPlus | `6_Cargo_vs_{card,vfdb,bacmet,tadb}.tsv`; `6_Background_vs_{card,vfdb,bacmet,tadb}.tsv`; `6_AMRF_Cargo.tsv`; `6_AMRF_Background.tsv` |
| 10 | `scripts/10_function_inventory.py` | Summarize functional annotations among cargo proteins | `6a_Cargo_function_inventory.csv`; `6a_Cargo_top_genes.csv`; `6a_Summary.txt` |
| 11 | `scripts/11_enrichment.py` | Compare functional-database hits between cargo and background proteins | `7_Enrichment_table.csv`; `7_Summary.txt`; `7_Enrichment_table_p80q70.csv`; `7_Summary_p80q70.txt`; `7_Enrichment_figure.pdf`; `7_Enrichment_figure.png` |
| 12 | `scripts/12_amrf_enrichment.py` | Calculate AMRFinderPlus-based antibiotic-resistance enrichment | `6_AMRF_weighted_summary.txt` |
| 13 | `scripts/13_strict_analysis.py`; `slurm/13_strict_analysis.sbatch` | Run the adjacent-same-family cargo sensitivity analysis | `8_Cargo_unique.faa`; `8_Cargo_unique_map.csv`; `8_Background_unique.faa`; `8_Background_unique_map.csv`; `8_{Cargo,Background}_vs_{card,vfdb,bacmet,tadb}.tsv`; `8_Enrichment_table.csv`; `8_Summary.txt` |
| 14 | `scripts/14_strict_counts.py` | Summarize strict-analysis plasmid and protein counts | `8b_Strict_counts.txt` |
| 15 | `scripts/15_is_sequence_jumps.py` | Compare exact IS sequences and their outer 10-bp plasmid flanks | `10_Per_plasmid.csv`; `10_Summary.txt`; `10_Jump_examples.txt` |
| 16 | `scripts/16_final_figures.py` | Generate the final count-comparison and cargo-enrichment figures | `4b_Transposase_count_comparison.pdf`; `4b_Transposase_count_comparison.png`; `7_Enrichment_figure.pdf`; `7_Enrichment_figure.png` |

`environment.yml` and `environment-amrfinder.yml` list the software dependencies. Database versions and file checksums are recorded in `config/database_manifest.tsv`.
