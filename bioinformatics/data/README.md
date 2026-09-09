# Essential data supporting the bioinformatics analyses

Ha et al., *Transposon–plasmid nesting shapes population-level gene-dosage
responses to fluctuating environments* (Nature Ecology & Evolution, manuscript
NATECOLEVOL-25061819).

This folder contains the essential data **used or generated** by the
bioinformatics analyses reported in the manuscript (Fig. 1a–d, Supplementary
Fig. 1, Tables S4–S5, Supplemental Note 1). File names and step numbers match
the analysis code in `../code/pipeline/scripts/`.



## Contents

### Figure 1a/b source data (external inputs)

| File | Content | Supports |
|---|---|---|
| `1_Transposase_counts.csv.gz` | per-replicon transposase-annotation counts for the full curated plasmid dataset (63,425 replicons; original file name `transposase-counts.csv`, analysis input from R. Maddamsetti, [darwinian-circuit](https://github.com/rohanmaddamsetti/darwinian-circuit)) | Fig. 1a densities (53.2 plasmid / 10.6 chromosome annotations per Mbp), Fig. 1b hit counts |
| `1_PIRA_PCN_estimates.csv.gz` | published pseuPIRA plasmid copy-number estimates (16,861 rows: 12,006 plasmids from 4,644 genomes; original file `41467_2025_61205_MOESM3_ESM.csv`, Supplementary Information of Maddamsetti et al., *Nat. Commun.* 16:6023 (2025), DOI 10.1038/s41467-025-61205-2; downloaded 2026-08-31) | Fig. 1b PCN bins |

Note on Fig. 1b: joining these two tables on `SeqID` (replicon accession;
`AnnotationAccession` is the per-assembly GCF and must not be used as the join
key) reproduces the PCN-binned fractions with slightly different denominators
(1,923/2,322; 4,198/6,396; 166/2,124; 6/166 plasmids per bin = 82.8 / 65.6 /
7.8 / 3.6%) than the current figure labels (1,825/2,192; 3,972/6,012;
155/2,051; 6/163 = 83.3 / 66.1 / 7.6 / 3.7%), because the figure was drawn
with an earlier-coverage version of the PIRA table that has not been made
public. The published table included here is the traceable public source; the
qualitative pattern is unchanged. (11,010 plasmids join in total; 2 with
PIRACopyNumber ≥ 1000 fall outside the four displayed bins.)

### Dataset definition and deduplication

| File | Content | Supports |
|---|---|---|
| `1_Hit_plasmids.txt` | 26,693 plasmid accessions with ≥1 transposase-annotation hit (selected from `transposase-counts.csv`) | sequence-level dataset, SI Note 1 |
| `3_Seqhash_map.csv` | exact-sequence MD5 → accession(s) mapping (dedup 26,693 → 25,801 unique sequences) | denominators used throughout |
| `3_Dedup_summary.txt` | deduplication audit (892 exact duplicates in 503 multi-copy groups) | Methods |

### Sequence-level confirmation of the annotation screen

| File | Content | Supports |
|---|---|---|
| `4_ISEScan_fna_tsv.tar.gz` | 25,801 per-plasmid ISEScan v1.7.3 result tables (`.fna.tsv`: IS family, boundaries, complete/partial), including header-only tables for the 1,404 no-hit sequences | all sequence-level analyses |
| `4b_Transposase_count_comparison.csv` | per-accession annotation-based transposase counts vs ISEScan IS-copy counts | Supplementary Fig. 1 (94.5%, Spearman ρ = 0.92, log–log slope 0.96) |
| `4b_Compare_summary.txt` | summary statistics for the comparison | Supplementary Fig. 1 legend |
| `4c_Discordant_cds.csv`, `4c_Summary.txt` | breakdown of the 1,404 annotation-positive / ISEScan-negative plasmids (Tn3-family 39.4%, other TEs 50.5%, regex FP 2.5%, annotation drift 7.6%) | SI Note 1 |

### Cargo definition

| File | Content | Supports |
|---|---|---|
| `5_Cargo_table.csv.gz` | per-CDS cargo assignments under the primary (relaxed) definition: same-family IS pair span ≤ 25 kb, IS-overlapping CDS excluded; records IS family and flank context | all cargo numbers |
| `5strict_Cargo_table.csv.gz` | same, strict adjacent-only sensitivity definition | sensitivity |
| `5_Summary.txt`, `5strict_Summary.txt` | plasmid/protein totals (17,227/25,801 = 66.77%; 511,814 cargo proteins relaxed) and IS-pair span statistics validating the 25 kb threshold | SI Note 1, Table S5 legend |

### Functional annotation and enrichment

| File | Content | Supports |
|---|---|---|
| `6_AMRF_Cargo.tsv`, `6_AMRF_Background.tsv` | raw AMRFinderPlus outputs (NCBI Reference Gene Catalog 2026-01-21.1) | orthogonal AMR validation |
| `6_AMRF_weighted_summary.txt` | copy-weighted AMRFinderPlus enrichment (OR = 14.35) | SI Note 1 |
| `6a_Cargo_function_inventory.csv` | per-hit inventory of functional cargo (database, gene, identity, coverage, copies, plasmid, IS family) at the primary threshold | reported category counts (6,953 / 2,494 / 3,394 / 4,704 plasmids) |
| `6a_Cargo_top_genes.csv`, `6a_Summary.txt` | top cargo genes per category; inventory summary | text |
| `7_Enrichment_table.csv`, `7_Summary.txt` | **primary enrichment** (DIAMOND pident ≥ 90 / qcov ≥ 60; copy-weighted; Fisher exact, BH-FDR): cargo 7.29 / 2.70 / 2.25 / 2.44% vs background 0.52 / 0.80 / 1.30 / 0.94%; OR 14.90 / 3.44 / 1.74 / 2.62 | Fig. 1d, Table S5 |
| `7_Enrichment_table_p80q70.csv`, `7_Summary_p80q70.txt` | 80/70 threshold sensitivity | SI Note 1 |
| `8_Enrichment_table.csv`, `8_Summary.txt`, `8b_Strict_counts.txt` | strict-definition sensitivity (OR 8.62 / 4.07 / 1.73 / 2.05) and counts at both thresholds | SI Note 1 |

### Identical IS copies at distinct insertion contexts

| File | Content | Supports |
|---|---|---|
| `10_Summary.txt` | funnel: 19,137 (≥2 same-family copies) → 9,729 (exact-identical pair) → **9,244 plasmids (35.8% of 25,801)** with the same exact IS sequence at ≥2 distinct 10-bp outer-flank contexts | Fig. 1c, Table S4 |
| `10_Per_plasmid.csv` | per-plasmid results | Table S4 |
| `10_Jump_examples.txt` | representative example records (accession, IS family, distinct flanks) | sequence-level interpretation |

---

