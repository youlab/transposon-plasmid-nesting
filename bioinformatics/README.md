# Code and data for Figure 1

Code and data supporting Figure 1 of Ha et al., *Nature Ecology & Evolution*
(manuscript NATECOLEVOL-25061819).

## Contents

- `data/` — the complete essential-data package for the bioinformatics
  analyses (all key result tables supporting Fig. 1a–d, Supplementary Fig. 1,
  Tables S4–S5 and Supplemental Note 1; see `data/README.md` for the
  file-by-file documentation). The two Figure 1a/b source tables are also
  provided as plain CSVs for direct use by the scripts:
  - `data/transposase-counts.csv` — per-replicon transposase-annotation counts
    (63,425 replicons; analysis input from R. Maddamsetti,
    github.com/rohanmaddamsetti/darwinian-circuit).
  - `data/PIRA-PCN-estimates.csv` — published pseuPIRA plasmid copy-number
    estimates (16,861 rows), Supplementary Information of Maddamsetti et al.,
    Nat. Commun. 16:6023 (2025), DOI 10.1038/s41467-025-61205-2
    (file 41467_2025_61205_MOESM3_ESM.csv, downloaded 2026-08-31).
  (The same two tables also appear in gzipped form as `1_*.csv.gz` inside the
  essential-data set; the plain CSVs are byte-identical in content.)
- `code/fig1a.py` — Fig. 1a from `transposase-counts.csv` alone.
- `code/fig1b.py` — data-driven Fig. 1b from the two public source tables.
- `code/fig1d.py` — Fig. 1d enrichment bars from
  `data/7_Enrichment_table.csv`.
- `code/pipeline/` — reference copy of the full analysis pipeline
  (`scripts/`: scripts 01–16, `slurm/` jobs, conda environments) that generated
  every result table in `data/`; see `code/pipeline/README.md`. These scripts
  assume the original project layout (configurable via the
  `NEE_PROJECT_ROOT`/`NEE_RESULTS_DIR` environment variables) and are provided
  for provenance; the standalone scripts above are the ones to run against
  this package.

## How to run

```
cd code
python3 fig1a.py
python3 fig1b.py
python3 fig1d.py
```

Requires python3 with matplotlib. Figure files are not stored in this
package; running the scripts regenerates them next to the scripts as vector
PDFs with editable TrueType text, plus 300-dpi PNGs.

## Expected results

- Fig. 1a densities: 53.2 (plasmid) / 10.6 (chromosome) annotations per Mbp.
- Fig. 1b recomputed fractions: 1923/2322 = 82.8%, 4198/6396 = 65.6%,
  166/2124 = 7.8%, 6/166 = 3.6% for PCN bins <1, 1-10, 10-100, 100-1000.
- Fig. 1d enrichment bars: cargo 7.29/2.70/2.25/2.44% vs background
  0.52/0.80/1.30/0.94% (CARD/VFDB/BacMet/TADB).

## Version note

The manuscript's Fig. 1b labels (1825/2192, 3972/6012, 155/2051, 6/163)
were computed with an earlier-coverage PIRA table that is not public;
recomputation from the published table gives slightly different denominators,
with the qualitative pattern unchanged. 11,010 plasmids join in total; 2
plasmids with PIRACopyNumber >= 1000 fall outside the four displayed bins.

## Scope note

The `code/` scripts regenerate Figures 1a, 1b and 1d directly from `data/`.
Figure 1c is a hand-drawn schematic assembled in Illustrator (not
script-generated); its headline number — 9,244/25,801 = 35.8% of plasmids
with the same exact IS sequence at distinct insertion contexts — comes from
`data/10_Summary.txt` (per-plasmid records in `data/10_Per_plasmid.csv`).
The file-to-claim mapping for all tables is documented in `data/README.md`.
