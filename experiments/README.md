# Experimental data and figure-generation workflow

This directory contains the experimental data, processed tables, and plotting notebook for the experimental analyses in Ha et al., *Transposon-plasmid nesting enables fast response to fluctuating environments*.

The notebook supports the experimental analyses in **Figs. 3, 4, and 5b**, together with the associated Supplementary and Extended Data figures listed below. Figure assembly and panel lettering were performed separately from numerical analysis and plotting.

## Main analysis notebook

[`Fig3_Fig4_Fig5_SuppFig7_8_plotting_code.ipynb`](Fig3_Fig4_Fig5_SuppFig7_8_plotting_code.ipynb)

Run the notebook from the `experiments/` directory, executing analysis sections in order. Some later sections depend on dataframes or processed CSV files generated earlier. The notebook filename and some internal headings retain numbering from an earlier revision; use the **final figure mapping** below.

The analysis uses Python 3 and packages including `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `altair`, `statannotations`, and `marimo`, together with the Excel-reading dependencies needed by pandas.

## Experimental design

The primary fluctuating-selection experiments used four plasmid backbones (**p15A, pBR, CloDF, and pUC**) and five independent biological replicates per condition. Each experiment comprised three consecutive 7-day cycles:

- 5 daily passages without tetracycline;
- 2 daily passages with tetracycline;
- a 30,000-fold dilution at each daily passage.

The main comparisons are:

1. **Transposase on versus off:** matched hosts containing the same chromosomal `tetA/gfp` transposable element, differing in the presence of the Tn5 transposase gene.
2. **Plasmid maintenance:** kanamycin maintained (`+Kan`) versus omitted (`-Kan`) during fluctuating tetracycline selection.
3. **Measured outputs:** population-averaged plasmid copy number (PCN), total transposable-element copy number (TCN), GFP/OD600, and derived loss/gain rates.

The earlier chromosome-free transposase-off control shown in **Extended Data Fig. 5** is distinct from the matched-host transposase-off comparison in Fig. 3.

## Directory contents

| Path | Contents |
| --- | --- |
| `qPCR_withTransposase/` | Raw daily qPCR data for transposase-on experiments (`+Kan` and `-Kan`), including calibration data |
| `qPCR_withoutTransposase/` | Raw daily qPCR data for transposase-off experiments (`+Kan` and `-Kan`), including calibration/control data |
| `GFP_withTransposase/` | Daily plate-reader data and plate-layout information for transposase-on experiments |
| `GFP_withoutTransposase/` | Daily plate-reader data and plate-layout information for transposase-off experiments |
| `9.12.26_Yuanchi_growthcurves.xlsx` | Short-term growth-curve source data for Figs. 3b and 4b |
| `PlotCalculatedRates/` | Processed GFP and qPCR tables used in response-rate analyses |
| `Figures_final/` | Saved PDF panels used during manuscript figure assembly |

Root-level processed CSV files retained for downstream analysis and plotting include:

```text
withTransposase_qPCR.csv
withTransposase_GFP-OD.csv
qPCR_no-transposase_with-Kan.csv
qPCR_no-transposase_no-Kan.csv
noTransposase_GFP-OD_r.csv
withJumping_gain_qPCR.csv
withJumping_loss_qPCR.csv
noJumping_gain_qPCR_r.csv
noJumping_loss_qPCR_r.csv
```

These intermediate files allow the downstream plotting sections to be run without reprocessing every raw instrument file.

## Final manuscript figure mapping

### Figure 3 — Transposase-on versus transposase-off comparisons

| Panel | Analysis |
| --- | --- |
| **3b** | Short-term growth curves for transposase-on and transposase-off strains |
| **3c–d** | Total TCN trajectories and loss/gain rates under fluctuating tetracycline selection |
| **3e–f** | GFP/OD600 trajectories and loss/gain rates |

Representative saved outputs are in `Figures_final/Figure3/`, including:

```text
Figure3B_+transposase.pdf
Figure3B_-transposase.pdf
Fig3C_transposase_withK_TCN.pdf
Fig3C_No_transposase_withK_TCN.pdf
Fig3E_transposase.pdf
Fig3E_noTransposase.pdf
Figure3DF.pdf
```

### Figure 4 — With and without kanamycin maintenance

| Panel | Analysis |
| --- | --- |
| **4b** | Short-term growth curves under `+Kan +Tet` and `-Kan +Tet` |
| **4c–d** | PCN trajectories and loss/gain rates |
| **4e–f** | Total TCN trajectories and loss/gain rates |
| **4g–h** | GFP/OD600 trajectories and loss/gain rates |

Representative saved outputs are in `Figures_final/Figure4/`, including:

```text
Figure4B_+Kan_+Tet.pdf
Figure4B_-Kan_+Tet.pdf
Fig4C_withK_PCN.pdf
Fig4C_noK_PCN.pdf
Fig4E_withK_TCN.pdf
Fig4E_noK_TCN.pdf
Fig4G_withKan_noChr.pdf
Fig4G_noKan.pdf
Figure4DFH.pdf
```

### Figure 5b — TCN response rates versus measured mean PCN

Each point pairs one biological replicate's TCN response rate with the arithmetic mean of its available daily PCN measurements under the same plasmid and experimental condition. **Figure 5b replots replicate-level data used in Figs. 3 and 4 and Extended Data Fig. 6; it is not an independent validation dataset.**

Saved outputs are in `Figures_final/Figure5B/`:

```text
F5B_gain_noT.pdf
F5B_loss_noT.pdf
F5B_gain_withT.pdf
F5B_loss_withT.pdf
```

### Supplementary and Extended Data figures

| Final figure | Experimental analysis |
| --- | --- |
| **Supplementary Fig. 5** | GFP/OD600 versus qPCR-measured total TCN |
| **Supplementary Fig. 6** | Gel-based plasmid-form assessment (panel a) and population-level TCN/PCN trajectories (panel b); the gel and final composite figure require separate figure assembly |
| **Supplementary Fig. 7** | TCN/PCN responses with and without kanamycin maintenance in **transposase-on** strains |
| **Supplementary Fig. 8** | qPCR assay characterization and primer-efficiency calibration |
| **Extended Data Fig. 5** | Earlier **chromosome-free**, transposase-off control |
| **Extended Data Fig. 6** | PCN, total TCN, and GFP/OD600 responses with and without kanamycin maintenance in **matched transposase-off** strains |

The primary plotting notebook contains the analyses for the main experimental figures and the associated time-course/rate and reporter comparisons. Not every panel in the final assembled Supplementary and Extended Data figures is generated directly by that notebook.

#### Historical numbering in notebook headings and saved files

Do **not** interpret older notebook section labels or filenames as final manuscript figure numbers:

| Notebook label or saved-output prefix | Final manuscript designation |
| --- | --- |
| Notebook section `Supplementary Figure 4` | **Supplementary Fig. 5** (GFP/OD600 versus qPCR TCN) |
| Saved files with prefix `SF7` | **Supplementary Fig. 7** (transposase-on, ±Kan, TCN/PCN) |
| Saved files with prefix `SF8` | **Extended Data Fig. 6** (transposase-off, ±Kan) |

These mappings update documentation only; they do not change the measurements or analysis.

## Response-rate calculations

For each biological replicate and each 7-day cycle, loss and gain rates are calculated from the logarithmic change of the measured quantity `X` (PCN, total TCN, TCN/PCN, or GFP/OD600), normalized by the number of generation equivalents implied by the serial dilution:

```text
loss rate = -ln(X_day5 / X_day0) / [5 * log2(30000)]
gain rate =  ln(X_day7 / X_day5) / [2 * log2(30000)]
```

Rates from the three cycles are averaged within each biological replicate. The response-rate plots show individual biological replicates, with mean ± standard deviation summaries where indicated. Pairwise condition comparisons use Mann–Whitney U tests.

Total TCN includes the chromosomal contribution in the matched hosts. Consequently, `TCN/PCN` is a **population-level ratio**, not a direct measurement of the number of transposon insertions in each plasmid.

## Reproducibility notes

- Run the plotting notebook from `experiments/`, so that relative data paths resolve.
- Execute the upstream processing sections before downstream plotting sections, or use the included processed CSV files where appropriate.
- The raw qPCR and plate-reader records and relevant processed tables are retained alongside the plotted outputs.
- Saved output filenames may retain historical figure numbers; use the final mapping above.
- Final panel lettering, combination of plots and gel images, and manuscript figure layout were performed separately from the analysis notebook.
