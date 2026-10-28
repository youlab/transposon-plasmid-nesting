# Processed data and figure inputs

This folder contains data used by the scripts, not rendered figures or execution logs. The large CDS table remains compressed. Original raw sequences, databases and temporary search files are not bundled.

| File | Content / use |
|---|---|
| `01_transposase_counts.csv` | Per-replicon transposase-annotation counts; Figure 1a/b and accession selection. |
| `01_plasmid_copy_number_estimates.csv` | Published PIRA copy-number estimates; Figure 1b joins by replicon SeqID. |
| `01_transposase_positive_accessions.txt` | Accessions selected for sequence-level processing. |
| `03_plasmid_sequence_hashes.csv` | Exact sequence deduplication map; 26,693 accessions / 25,801 representatives. |
| `03_deduplication_summary.txt` | Sequence-deduplication totals. |
| `04_isescan_results.tar.gz` | Per-representative ISEScan TSVs; extract only when running upstream steps. |
| `05_transposase_count_comparison.csv` | Per-accession annotation counts versus ISEScan IS counts. |
| `05_transposase_count_summary.txt` | Summary of the IS-count comparison. |
| `06_discordant_cds.csv` | Detailed annotation-positive / ISEScan-negative records. |
| `06_discordant_summary.txt` | Summary of the discordant annotations. |
| `07_cargo_summary_relaxed.txt; 07_cargo_summary_strict.txt` | Current cargo/background counts from the supplied CDS membership input. |
| `09_amrfinder_cargo.tsv; 09_amrfinder_background.tsv` | Original AMRFinderPlus outputs; not backbone-filtered reruns. |
| `10_cargo_function_inventory.csv; 10_cargo_top_genes.csv; 10_cargo_function_summary.txt` | Original functional cargo inventory and reported gene/carrier summaries. |
| `11_enrichment_relaxed.csv` | Current original four-category relaxed enrichment; Figure 1d source. |
| `13_enrichment_strict.csv` | Current original four-category strict sensitivity result. |
| `11_enrichment_relaxed_summary.txt; 13_enrichment_strict_summary.txt` | Readable summaries of those current enrichment results. |
| `11_enrichment_relaxed_p80q70.csv; 11_enrichment_relaxed_p80q70_summary.txt` | Original 80% identity / 70% query-coverage sensitivity output; a separate input/threshold analysis, not regenerated from the supplied 90/60 hit flags. |
| `12_amrfinder_enrichment_summary.txt` | Original AMRFinderPlus sensitivity summary; not rerun after backbone filtering. |
| `15_is_insertion_contexts.csv; 15_is_insertion_summary.txt; 15_is_insertion_examples.txt` | Per-plasmid identical-IS / distinct-context analysis supporting Figure 1c and Table S4. |
| `19_cds_annotation_membership.csv.gz` | One row per eligible CDS per cargo mode: accession, coordinates, PlasAnn category, cargo/background group, backbone/retention flags and all four original functional-hit flags. Contains both retained and excluded rows; read directly as gzip. |
| `20_cds_counts_before_filtering.csv` | Integer counts and denominators for the original four categories, for both groups and cargo modes. |
| `20_cds_counts_after_filtering.csv` | Integer retained counts and denominators for the 4 backbone-filtered comparison categories. |
| `20_backbone_filter_counts.csv` | Raw, excluded and retained CDS counts, separately by group and cargo mode. |
| `21_backbone_enrichment.csv` | Two-sided pooled-CDS Fisher and BH(4) results; relaxed rows supply article Extended Figure 1, while strict rows are SI sensitivity results where cited. |

## Scientific scope

The supplied CDS input includes 3,323,101 eligible CDSs per cargo mode (6,646,202 rows across relaxed and strict), not 6,646,202 distinct proteins. Corrected raw cargo/background totals are 511,814 / 2,811,287 (relaxed) and 240,858 / 3,082,243 (strict). After backbone exclusion, totals are 469,783 / 2,484,691 and 221,029 / 2,733,445, respectively. Filtered odds ratios (CARD, VFDB, BacMet, TADB) are 14.444048 / 3.467556 / 1.679686 / 2.601322 for relaxed and 8.389165 / 4.100302 / 1.673068 / 2.019852 for strict. The strict comparison is an SI sensitivity analysis only where cited; article Extended Figure 1 uses the relaxed results.

Original comparisons always use four functional categories. The WithTA backbone-filtered comparison uses 4. TA-positive CDSs remain in the denominators unless excluded by the independent backbone rule.

The supplied membership master and six derived CSVs are the corrected current WithTA outputs. The 07 extractor checks strict adjacent pairs and the circular last/first pair, and includes CDS coordinates equal to either wrap boundary. The correction changes cargo/background group labels only; CDS coordinates, PlasAnn categories and all four functional-hit flags are unchanged. Separate 80/70 and AMRFinder sensitivity outputs retain their own analysis scope; their source searches cannot be recreated by changing the supplied 90/60 Boolean hit flags. Relative to the previous local 9/8 strict output, 41 additional wrap-qualified CDSs are included; none has a hit in any of the four functional databases. Five same-locus records differ only in coordinates and remain present.

The two public Figure 1a/b source tables are included once per independent package; duplicate gzip copies were removed. Joining them on SeqID gives PCN-bin fractions 1,923/2,322, 4,198/6,396, 166/2,124 and 6/166; 11,010 records join, with two >=1000 outside the plotted bins.
