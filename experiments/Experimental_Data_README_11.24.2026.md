# Experimental data and figure analysis

This archive contains qPCR, GFP/OD600, and growth-curve measurements; processed tables; and the analysis notebook used for experimental figures in *Transposon nesting within plasmids enables fast response to fluctuating environments*.

## Analysis notebook

[`Fig3-5_ED5-6_SF5-7_plotting_code.ipynb`](Fig3-5_ED5-6_SF5-7_plotting_code.ipynb)

Run the notebook from the directory containing this README so that its relative input paths resolve. Its sections include calculations and plots for Figs. 3, 4, and 5b, Extended Data Figs. 5–6, and associated Supplementary Figures. Final multi-panel figure assembly was performed separately. Historical names may remain in intermediate CSVs and saved PDF panels.

## Files and directories

| Path | Contents |
| --- | --- |
| `qPCR_withTransposase/` | Transposase-on qPCR instrument files, calibration, and associated data |
| `qPCR_withoutTransposase/` | Transposase-off qPCR instrument files, calibration, and chromosome-free control data |
| `GFP_withTransposase/`, `GFP_withoutTransposase/` | Daily fluorescence and OD600 records and plate metadata |
| `9.12.26_Yuanchi_growthcurves.xlsx` | Short-term growth-curve measurements |
| Root-level `.csv` files | Processed qPCR, fluorescence, and replicate-level rate tables |
| `Figures_final/` | Individual saved PDF panels, organized by figure |

The transposase-on calibration is `qPCR_withTransposase/2024-1001_yuanchi.xls`. The transposase-off calibration is `qPCR_withoutTransposase/20241125-calibration_tetAFAM_dxsNED_kanRVIC.xls`. In transposase-off instrument exports, the target named `cmR` represents the `dxs` assay; each host's qPCR calculations use its corresponding calibration factors.

## Experimental conditions and figure mapping

Primary experiments used p15A, pBR, CloDF, and pUC plasmids, five biological replicates per condition, and three 7-day cycles (five passages without tetracycline followed by two with tetracycline). Each daily passage used a 30,000-fold dilution.

| Figure | Data or comparison |
| --- | --- |
| Fig. 3 | Matched transposase-on and transposase-off hosts: total TCN, GFP/OD600, growth curves, and response rates |
| Fig. 4 | Transposase-on +Kan versus −Kan: PCN, total TCN, GFP/OD600, growth curves, and rates |
| Fig. 5b | Replicate-level TCN response rates versus measured mean PCN |
| Extended Data Fig. 5 | Earlier transposase-off, chromosome-free host control |
| Extended Data Fig. 6 | Matched transposase-off +Kan versus −Kan: PCN, total TCN, GFP/OD600, and rates |
| Supplementary Figs. 5–7 | GFP/OD600–TCN comparison, plasmid-associated TCN/PCN, and maintenance-condition comparisons |

For each biological replicate, response rates are computed for each cycle from the natural-log change in the measured quantity divided by `5 × log2(30000)` (loss) or `2 × log2(30000)` (gain), and then averaged across cycles. Mean PCN for Fig. 5b is the arithmetic mean of available daily PCN values for the corresponding replicate and condition. The matched hosts retain one chromosomal transposable-element copy in total TCN measurements; plasmid-associated TCN/PCN analyses account for that chromosomal contribution.
