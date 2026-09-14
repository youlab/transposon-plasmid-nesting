# Figures 3, 4, 5B, Supplementary Figures 4, 7, and 8 — Plotting Notebook

## Overview

This notebook contains the analysis and plotting workflow used to generate the manuscript figures associated with the plasmid copy-number, transposon copy-number, GFP/OD, and transposon-to-plasmid (`T/P`) measurements.

The notebook is organized into major sections corresponding to the experimental condition and analysis stage:

1. **qPCR without Transposase**
2. **GFP without Transposase**
3. **qPCR with Transposase**
4. **GFP with Transposase**
5. **Plot Calculated Rates**
6. **Supplementary Figure 4**
7. **Rate Plots for Figures 3, 4, Supplementary Figures 7 and 8**
8. **Figure 5B**
9. **Figure 3B**
10.**Figure 4B**

The workflow combines raw qPCR and plate-reader measurements, applies qPCR calibration, calculates plasmid/transposon copy numbers and `T/P`, calculates gain/loss rates, and generates the final PDF figures.

> **Important:** Run each major section from top to bottom. Later sections depend on dataframes, processed CSV files, and functions created earlier in the notebook.

---

# 1. qPCR WITHOUT TRANSPOSASE

This section processes the qPCR time-course data collected without transposase.

## 1.1 Imports and global parameters

The main experimental parameters are defined here:

- Up/selection days: `0, 7, 14, 21`
- Down/no-selection days: `5, 12, 19`
- Treatment dilution: `30000`
- Loss-treatment duration: `5` days
- Selection-treatment duration: `2` days
- Plasmid/PCN categories: p15A, pBR, CloDF, pUC

Treatment windows are represented as:

```text
Day 5–7
Day 12–14
Day 19–21
```

These windows are shaded on the time-course plots.

## 1.2 qPCR standard-curve calibration

The notebook first reads the qPCR calibration file:

```text
qPCR_withTransposase/2024-1001_yuanchi.xls
```

Calibration data are separated by:

- `cmR` → chromosome target
- `tetA` → transposon target
- `kanR` → plasmid target

Standard-curve plots are generated for both pUC and p15A calibration material.

The calibration function calculates:

```text
slope
amplification factor = 10^(-1/slope)
amplification efficiency = (amplification factor - 1) × 100
```

These calibration parameters are subsequently applied to the daily qPCR measurements.

## 1.3 Import and process the 22-day qPCR time course

Daily qPCR files are imported for:

- `Kan`
- `noKan`

covering Days 0–21.

For each day, the notebook:

1. Assigns the appropriate calibration slope.
2. Assigns amplification factors and efficiencies.
3. Removes `H2O` and `Blank` wells.
4. Organizes measurements by qPCR replicate, experimental replicate, and PCN/plasmid.
5. Calculates plasmid copy number.
6. Calculates transposon copy number.
7. Calculates the transposon-to-plasmid ratio (`T/P`).
8. Groups the resulting measurements by day, experimental replicate, and PCN.

Calibration/testing strains are removed before downstream analysis.
The `p15A++` condition is kept separate because it is p15A with transposase.

## 1.4 Add chromosome controls

Chromosome control measurements are imported from:

```text
qPCR_withoutTransposase/noTransposon_kan_qPCR_timeCourse.csv
qPCR_withoutTransposase/noTransposon_noKan_qPCR_timeCourse.csv
```

These chromosome measurements are combined with the plasmid measurements where appropriate.

## 1.5 qPCR time-course plots

The processed qPCR data are plotted across the 22-day experiment.

The plots show individual experimental replicate trajectories using:

- x-axis: day
- y-axis: log plasmid copy number, log transposon copy number, or `T/P`
- color: plasmid
- markers: replicate measurements
- shaded treatment windows

### Figures generated

| Figure | Panel | Output |
|---|---|---|
| **Figure 3** | 3C | `Figures_final/Figure3/Fig3C_No_transposase_withK_TCN.pdf` |
| **Supplementary Figure 8** | 8A | `Figures_final/supplementary/SF8A_No_transposase_withK_PCN.pdf` |
| **Supplementary Figure 8** | 8A | `Figures_final/supplementary/SF8A_No_transposase_noK_PCN.pdf` |
| **Supplementary Figure 8** | 8C | `Figures_final/supplementary/SF8C-No_transposase_K_TCN.pdf` |
| **Supplementary Figure 8** | 8C | `Figures_final/supplementary/SF8C-No_transposase_noK_TCN.pdf` |

The section also creates intermediate CSV files for later rate calculations:

```text
qPCR_no-transposase_with-Kan.csv
qPCR_no-transposase_no-Kan.csv
qPCR_transposase_with-Kan_p15A.csv
qPCR_transposase_no-Kan_p15A.csv
```

---

# 2. GFP WITHOUT TRANSPOSASE

This section processes the GFP/OD time-course data collected without transposase.

## 2.1 Import daily plate-reader data

The notebook imports the daily GFP/OD measurements from the `GFP_withoutTransposase` directory.

The data include:

- OD600
- GFP
- plasmid identity
- experimental replicate
- day
- Kan/noKan condition

Plate-layout metadata are read from the corresponding metadata CSV files.

## 2.2 Prepare plate data

Helper functions convert the plate-reader files into a consistent long-format dataframe.

The workflow:

1. Reads plate-layout metadata.
2. Converts well positions into sample identifiers.
3. Extracts OD600 measurements.
4. Extracts GFP measurements.
5. Combines GFP and OD measurements.
6. Associates each measurement with its plasmid and experimental replicate.
7. Labels Kan/noKan conditions.
8. Adds treatment/selection information.

The `p15A++` condition is kept separate because it is p15A with transposase.

## 2.3 Generate GFP/OD time courses

GFP/OD is plotted as:

```text
log10(GFP/OD)
```

Individual experimental replicate trajectories are retained rather than averaging the lines together.

Treatment windows are shaded at:

```text
Days 5–7
Days 12–14
Days 19–21
```

### Figures generated

| Figure | Panel | Output |
|---|---|---|
| **Figure 3** | 3E | `Figures_final/Figure3/Fig3E_noTransposase.pdf` |
| **Supplementary Figure 8** | 8E | `Figures_final/supplementary/SF8E_noTransposase_Kan.pdf` |
| **Supplementary Figure 8** | 8E | `Figures_final/supplementary/SF8E_noTransposase_noKan.pdf` |

The processed GFP data are also saved as:

```text
noTransposase_GFP-OD_r.csv
```

---

# 3. qPCR WITH TRANSPOSASE

This section processes the qPCR measurements from the transposase-containing experiments.

## 3.1 qPCR calibration

The same general calibration strategy is used:

- chromosome → `cmR`
- transposon → `tetA`
- plasmid → `kanR`

The calibration curves are visualized before processing the daily data.

## 3.2 Import daily qPCR measurements

qPCR files are imported for Days 0–21 under both:

- `Kan`
- `noKan`

conditions.

The notebook also imports the separately generated p15A transposase data:

```text
qPCR_transposase_with-Kan_p15A.csv
qPCR_transposase_no-Kan_p15A.csv
```

## 3.3 Calculate copy numbers

The `processOneDay()` function applies the calibration parameters to each daily dataset.

The processed data contain measurements for:

- plasmid copy number
- transposon copy number
- `T/P`
- experimental replicate
- day
- PCN/plasmid
- Kan/noKan condition


## 3.4 Combine p15A and the remaining plasmids

The separately processed p15A data are combined with the pBR, CloDF, and pUC datasets.

This creates the primary qPCR dataframe used for the transposase-containing condition.

### Figures generated

| Figure | Panel | Output |
|---|---|---|
| **Figure 3** | 3C | `Figures_final/Figure3/Fig3C_transposase_withK_TCN.pdf` |
| **Figure 4** | 4C | `Figures_final/Figure4/Fig4C_withK_PCN.pdf` |
| **Figure 4** | 4E | `Figures_final/Figure4/Fig4E_withK_TCN.pdf` |
| **Supplementary Figure 7** | 7A | `Figures_final/supplementary/SF7A_withK_TP.pdf` |
| **Figure 4** | 4C | `Figures_final/Figure4/Fig4C_noK_PCN.pdf` |
| **Figure 4** | 4E | `Figures_final/Figure4/Fig4E_noK_TCN.pdf` |
| **Supplementary Figure 7** | 7A | `Figures_final/supplementary/SF7A_noK_TP.pdf` |

An intermediate processed dataset is saved as:

```text
withTransposase_qPCR.csv
```

---

# 4. GFP DATA WITH TRANSPOSASE

This section processes the GFP/OD data corresponding to the transposase-containing experiments.

## 4.1 Import and organize GFP data

Daily GFP/OD measurements are imported from the `GFP_withTransposase` datasets.

Plate metadata are used to identify:

- plasmid
- experimental replicate
- day
- Kan/noKan condition

A chromosome control is added from the previously generated no-transposase dataset.

## 4.2 Handle the transposase-containing p15A condition

The p15A transposase-containing dataset is incorporated into the combined dataframe.

## 4.3 Generate GFP/OD time-course plots

The primary measurement plotted is:

```text
log10(GFP/OD)
```

The plots retain individual experimental replicate trajectories.

### Figures generated

| Figure | Panel | Output |
|---|---|---|
| **Figure 3** | 3E | `Figures_final/Figure3/Fig3E_transposase.pdf` |
| **Figure 4** | 4G | `Figures_final/Figure4/Fig4G_withKan_noChr.pdf` |
| **Figure 4** | 4G | `Figures_final/Figure4/Fig4G_noKan.pdf` |

The combined processed GFP data are also written to:

```text
withTransposase_GFP-OD.csv
```

---

# 5. PLOT CALCULATED RATES

This section collects the calculated rate datasets produced by the preceding analyses and prepares them for the final statistical plots.

## 5.1 Import rate datasets

The notebook imports GFP rate data for:

- no transposase / no jumping
- transposase / jumping
- gain
- loss

and corresponding qPCR rate data.

Examples include:

```text
noJumping_gain_r.csv
noJumping_loss_r.csv
withJumping_gain.csv
withJumping_loss.csv
```

The qPCR rate datasets are similarly imported.

---

# 6. SUPPLEMENTARY FIGURE 4

Supplementary Figure 4 examines the relationship between GFP/OD and transposon copy number.

Two datasets are compared:

1. **No jumping / no transposase**
2. **Jumping / transposase**

The notebook creates scatter plots of:

```text
x = log10(GFP/OD)
y = log(transposon copy)
```

Both full ranges and restricted positive-transposon-copy versions are generated.

The final combined plots are constructed by adding the corresponding no-jumping and jumping charts.

### Figure generated

**Supplementary Figure 4**

The notebook constructs the figure interactively through the Altair objects:

```text
SF4A
SF4B
J
J_B
SF4A_revision
combined_r
```

---

# 7. COMBINED RATE PLOTS — FIGURES 3, 4, SUPPLEMENTARY FIGURES 7 AND 8

This is the main statistical plotting section of the notebook.

## 7.1 Statistical plotting functions

Two main plotting functions are defined:

```text
plotStatQPCR()
plotStatGFP()
```

Both functions:

- separate measurements by plasmid
- separate Kan and noKan conditions
- plot individual replicate values
- calculate the mean
- show standard deviation
- overlay the mean point
- use the plasma color scheme for plasmids
- use silver for chromosome where applicable
- compare paired plasmid conditions
- add statistical significance annotations

The statistical comparison is performed using the **Mann–Whitney test**, with significance displayed using star notation.

## 7.2 Figure 3D/F

### Figure 3D

qPCR:

- rate measurements
- loss
- gain

### Figure 3F

GFP:
- rate measurements
- loss
- gain

The combined four-panel figure is also saved as:

```text
Figures_final/Figure3/Figure3DF.pdf
```

### Figure 3 rate outputs

```text
Fig3F_GFP_Loss.pdf
Fig3D_qPCR_Loss.pdf
Fig3F_GFP_Gain.pdf
Fig3D_qPCR_Gain.pdf
Figure3DF.pdf
```

## 7.3 Figure 4D/F/H

The three rows correspond to:

### Figure 4D — plasmid copy-number rate

```text
rate_log(plasmid copy)
```

### Figure 4F — transposon copy-number rate

```text
rate_log(transposon copy)
```

### Figure 4H — GFP rate

```text
rate_log(GFP/OD)
```

Both gain and loss are shown.

The complete combined figure is saved as:

```text
Figures_final/Figure4/Figure4DFH.pdf
```

### Figure 4 rate outputs

```text
Fig4H_GFP_Gain.pdf
Fig4F_qPCR_Gain.pdf
Fig4D_qPCR_PCN_Gain.pdf

Fig4H_GFP_Loss.pdf
Fig4F_qPCR_Loss.pdf
Fig4D_qPCR_PCN_Loss.pdf

Figure4DFH.pdf
```

## 7.4 Supplementary Figure 8B/D/F

These plots repeat the rate analysis for the **no-transposase** condition.

### Supplementary Figure 8B

Plasmid copy-number gain/loss rates.

### Supplementary Figure 8D

Transposon copy-number gain/loss rates.

### Supplementary Figure 8F

GFP gain/loss rates.

The combined figure is saved as:

```text
Figures_final/supplementary/SF8BDF.pdf
```

### Outputs

```text
SF8B_plasmid_qPCR_Loss.pdf
SF8B_plasmid_qPCR_Gain.pdf

SF8D_qPCR_Loss.pdf
SF8D_qPCR_Gain.pdf

SF8F_GFP_Loss.pdf
SF8F_GFP_Gain.pdf

SF8BDF.pdf
```

## 7.5 Supplementary Figure 7B

Supplementary Figure 7B focuses on the transposon-to-plasmid ratio:

```text
rate_log(T/P)
```

Gain and loss are shown separately.

The combined figure is saved as:

```text
Figures_final/supplementary/SF7B.pdf
```

### Outputs

```text
SF7B_plasmid_TP_Loss.pdf
SF7B_plasmid_TP_Gain.pdf
SF7B.pdf
```

---

# 8. FIGURE 5B

The section first reconstructs the data needed for the four Figure 5B conditions.

## 8.1 Data groups

The GFP and qPCR datasets are divided into:

```text
GFP_gain_noT
GFP_gain_withT
GFP_loss_noT
GFP_loss_withT

qPCR_gain_noT
qPCR_gain_withT
qPCR_loss_noT
qPCR_loss_withT
```

## 8.2 Replicate-level average plasmid copy number

For the qPCR data, plasmid copy number is averaged within each:

```text
PCN + experimental replicate
```

combination across days.

This creates:

```text
replicate_averagePCN
```

This variable is then assigned back to the gain/loss qPCR datasets so that each rate measurement has a corresponding replicate-level average plasmid copy number.

## 8.3 Figure 5B plotting style

The final Figure 5B plots use:

- x-axis: replicate-average plasmid copy number
- logarithmic x-axis
- y-axis: transposon-copy gain/loss rate
- plasmid-specific plasma colors
- circular markers
- open circles for `-Kan`
- filled circles for `+Kan`
- individual replicate points
- no grid
- fixed x-axis major ticks

The four conditions are:

### Figure 5B — gain, no transposase

```text
F5B_gain_noT.pdf
```

### Figure 5B — loss, no transposase

The corresponding code is present, but the final `savefig()` line is currently commented out in the notebook.

Intended filename:

```text
F5B_loss_noT.pdf
```

### Figure 5B — gain, transposase

```text
F5B_gain_withT.pdf
```

### Figure 5B — loss, transposase

```text
F5B_loss_withT.pdf
```

---

# 9. FIGURE 3B, 4B

The growth curve data is plotted for strains in the appropriate antibiotic conditions.

### Figure 3B +transposase
```text
Figure3B_+transposase.pdf
```

### Figure 3B -transposase
```text
Figure3B_-transposase.pdf
```

### Figure 4B +Kan, +Tet
```text
Figure4B_+Kan_+Tet.pdf
```

### Figure 4B -Kan, +Tet
```text
Figure4B_-Kan_+Tet.pdf
```

---

# 10. FIGURE OUTPUT SUMMARY

## Main Figures

### Figure 3

| Panel | Description | Output |
|---|---|---|
| 3C | qPCR transposon copy, no transposase | `Fig3C_No_transposase_withK_TCN.pdf` |
| 3C | qPCR transposon copy, transposase | `Fig3C_transposase_withK_TCN.pdf` |
| 3D | qPCR loss/gain rates | `Fig3D_qPCR_Loss.pdf`, `Fig3D_qPCR_Gain.pdf` |
| 3E | GFP/OD, no transposase | `Fig3E_noTransposase.pdf` |
| 3E | GFP/OD, transposase | `Fig3E_transposase.pdf` |
| 3F | GFP loss/gain rates | `Fig3F_GFP_Loss.pdf`, `Fig3F_GFP_Gain.pdf` |
| 3D/F | Combined rate figure | `Figure3DF.pdf` |

### Figure 4

| Panel | Description | Output |
|---|---|---|
| 4C | Plasmid copy, +Kan | `Fig4C_withK_PCN.pdf` |
| 4C | Plasmid copy, −Kan | `Fig4C_noK_PCN.pdf` |
| 4D | Plasmid copy-number gain/loss rates | `Fig4D_qPCR_PCN_Gain.pdf`, `Fig4D_qPCR_PCN_Loss.pdf` |
| 4E | Transposon copy, +Kan | `Fig4E_withK_TCN.pdf` |
| 4E | Transposon copy, −Kan | `Fig4E_noK_TCN.pdf` |
| 4F | Transposon-copy gain/loss rates | `Fig4F_qPCR_Gain.pdf`, `Fig4F_qPCR_Loss.pdf` |
| 4G | GFP/OD, +Kan | `Fig4G_withKan_noChr.pdf` |
| 4G | GFP/OD, -Kan | `Fig4G_noKan.pdf` |
| 4H | GFP gain/loss rates | `Fig4H_GFP_Gain.pdf`, `Fig4H_GFP_Loss.pdf` |
| 4D/F/H | Combined rate figure | `Figure4DFH.pdf` |

### Figure 5B

| Condition | Output |
|---|---|
| Gain, −T | `F5B_gain_noT.pdf` |
| Loss, −T | `F5B_loss_noT.pdf` |
| Gain, +T | `F5B_gain_withT.pdf` |
| Loss, +T | `F5B_loss_withT.pdf` |

### Figure 3B

| Condition | Output |
|---|---|
| +T | `Figure3B_+transposase.pdf` |
| −T | `Figure3B_-transposase.pdf` |

### Figure 4B

| Condition | Output |
|---|---|
| +Kan, +Tet | `Figure4B_+Kan_+Tet.pdf` |
| -Kan, +Tet | `Figure4B_-Kan_+Tet.pdf` |


---

# 11. SUPPLEMENTARY FIGURE OUTPUT SUMMARY

## Supplementary Figure 4

Relationship between GFP/OD and transposon copy number for jumping and no-jumping datasets.


## Supplementary Figure 7

### 7A

Time-course measurements of `T/P`.

```text
SF7A_withK_TP.pdf
SF7A_noK_TP.pdf
```

### 7B

Gain/loss rate analysis of `T/P`.

```text
SF7B_plasmid_TP_Loss.pdf
SF7B_plasmid_TP_Gain.pdf
SF7B.pdf
```

## Supplementary Figure 8

### 8A

Plasmid copy-number time courses without transposase.

```text
SF8A_No_transposase_withK_PCN.pdf
SF8A_No_transposase_noK_PCN.pdf
```

### 8B

Plasmid copy-number gain/loss rates without transposase.

```text
SF8B_plasmid_qPCR_Loss.pdf
SF8B_plasmid_qPCR_Gain.pdf
```

### 8C

Transposon copy-number time courses without transposase.

```text
SF8C-No_transposase_K_TCN.pdf
SF8C-No_transposase_noK_TCN.pdf
```

### 8D

Transposon-copy gain/loss rates without transposase.

```text
SF8D_qPCR_Loss.pdf
SF8D_qPCR_Gain.pdf
```

### 8E

GFP/OD time courses without transposase.

```text
SF8E_noTransposase_Kan.pdf
SF8E_noTransposase_noKan.pdf
```

### 8F

GFP gain/loss rates without transposase.

```text
SF8F_GFP_Loss.pdf
SF8F_GFP_Gain.pdf
SF8BDF.pdf
```

---

# 12. INPUT DATA ORGANIZATION

The notebook expects the input data to be available using the relative paths embedded in the code.

## qPCR — no transposase

```text
qPCR_withoutTransposase/
├── daily qPCR .xls files for Days 0–21
├── noTransposon_kan_qPCR_timeCourse.csv
├── noTransposon_noKan_qPCR_timeCourse.csv
└── 2024-1001_yuanchi.xls
```

## qPCR — transposase

```text
qPCR_withTransposase/
├── daily qPCR .xls files for Days 0–21
├── noTransposon_kan_qPCR_timeCourse.csv
├── noTransposon_noKan_qPCR_timeCourse.csv
└── 2024-1001_yuanchi.xls
```

## GFP — no transposase

```text
GFP_withoutTransposase/
├── daily Day 0–21 .xlsx files
├── plate-layout-metadata.csv
├── plate-layout-metadata_Day0.csv
├── plate-layout-metadata_tet10.csv
└── noJumping_rawData_GFP.csv
```

## GFP — transposase

```text
GFP_withTransposase/
├── daily Day 0–21 .xlsx files
├── plate-layout-metadata.csv
├── plate-layout-metadata_Day0.csv
├── plate-layout-metadata_tet10.csv
└── noJumping_rawData_GFP.csv
```

## Previously processed datasets

The notebook also uses processed CSV files from previous analysis stages, including:

```text
PlotCalculatedRates/
noJumping_gain_r.csv
noJumping_loss_r.csv
noJumping_gain_qPCR_r.csv
noJumping_loss_qPCR_r.csv
withJumping_gain.csv
withJumping_loss.csv
withJumping_gain_qPCR.csv
withJumping_loss_qPCR.csv
```

---

# 13. KEY DATAFRAME NAMING

The main dataframes used throughout the notebook are:

| Dataframe | Purpose |
|---|---|
| `completeDF_tet10_kan` | Processed qPCR data with Kan |
| `completeDF_tet10_noKan` | Processed qPCR data without Kan |
| `tet10_kan_DF_Use_grouped` | Grouped qPCR data with Kan |
| `tet10_noKan_DF_Use_grouped` | Grouped qPCR data without Kan |
| `combined_df` | Combined GFP rate dataset |
| `combined_qPCR_df` | Combined qPCR rate dataset |
| `combined_scatter_DF` | Combined loss-rate dataframe |
| `combined_scatter_gain_DF` | Combined gain-rate dataframe |
| `combined_df_Fig5` | GFP data prepared for Figure 5B |
| `combined_qPCR_df_Fig5` | qPCR data prepared for Figure 5B |

---


# 14. OUTPUT DIRECTORY STRUCTURE

The final figures are organized as:

```text
Figures_final/
├── Figure3/
│   ├── Fig3C_No_transposase_withK_TCN.pdf
│   ├── Fig3C_transposase_withK_TCN.pdf
│   ├── Fig3E_noTransposase.pdf
│   ├── Fig3E_transposase.pdf
│   ├── Figure3B_-transposase.pdf
│   ├── Figure3B_+transposase.pdf
│   └── Figure3DF.pdf
│
├── Figure4/
│   ├── Fig4C_withK_PCN.pdf
│   ├── Fig4C_noK_PCN.pdf
│   ├── Fig4E_withK_TCN.pdf
│   ├── Fig4E_noK_TCN.pdf
│   ├── Fig4G_withKan_noChr.pdf
│   ├── Fig4G_noKan.pdf
│   ├── Figure4B_-Kan_+Tet.pdf
│   ├── Figure4B_+Kan_+Tet.pdf
│   └── Figure4DFH.pdf
│
├── Figure5B/
│   ├── F5B_gain_noT.pdf
│   ├── F5B_loss_noT.pdf
│   ├── F5B_gain_withT.pdf
│   └── F5B_loss_withT.pdf
│
└── supplementary/
    ├── SF7A_withK_TP.pdf
    ├── SF7A_noK_TP.pdf
    ├── SF7B.pdf
    ├── SF8A_No_transposase_withK_PCN.pdf
    ├── SF8A_No_transposase_noK_PCN.pdf
    ├── SF8BDF.pdf
    ├── SF8C-No_transposase_K_TCN.pdf
    ├── SF8C-No_transposase_noK_TCN.pdf
    ├── SF8E_noTransposase_Kan.pdf
    └── SF8E_noTransposase_noKan.pdf
```
