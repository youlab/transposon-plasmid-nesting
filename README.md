# Transposon–plasmid nesting: analysis code and data

Code, processed data, and figure-generation workflows supporting *Transposon-plasmid nesting enables fast response to fluctuating environments* (Ha et al.).

- **[Bioinformatics](bioinformatics/README.md)** — annotation-based plasmid survey, insertion-sequence analysis, putative cargo enrichment, and PlasAnn backbone-filtered comparison (Fig. 1 and Extended Data Fig. 1).
- **[Experiments](experiments/README.md)** — raw and processed qPCR, fluorescence and growth data; figure-generation notebook for Figs. 3, 4 and 5b and related Supplementary/Extended Data analyses.
- **[Simulations](simulations/README_simulations.md)** — shared coarse-grained TCN-state model and figure-specific notebooks for Figs. 2 and 5a, Extended Data Figs. 2–4 and Supplementary simulation figures.

Each component has its own instructions, input-data descriptions, and figure-generation entry points. The simulation notebooks load `simulation_core.ipynb` from their common directory. The bioinformatics component includes compressed processed data for recalculating the reported backbone-filtered comparisons; large upstream sequence databases and intermediate files are not bundled.
