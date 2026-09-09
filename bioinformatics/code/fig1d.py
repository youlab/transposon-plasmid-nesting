#!/usr/bin/env python3
"""
fig1d.py — Fig. 1d: enrichment of functional categories among
cargo proteins vs background, drawn from the shipped result table.

Standalone port of the Fig. 1d panel in
pipeline/scripts/16_final_figures.py (fig_11), reading
../data/7_Enrichment_table.csv relative to this script. Values are not
hardcoded; expected bars: cargo 7.29/2.70/2.25/2.44% vs background
0.52/0.80/1.30/0.94% (see data/README.md).

Run:  python3 fig1d.py   (requires python3 with matplotlib)
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FixedLocator, FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = os.path.join(HERE, "..", "data", "7_Enrichment_table.csv")

FONT_FILE = "/System/Library/Fonts/Supplemental/Arial.ttf"
if os.path.exists(FONT_FILE):
    font_manager.fontManager.addfont(FONT_FILE)
    plt.rcParams["font.family"] = font_manager.FontProperties(fname=FONT_FILE).get_name()
else:
    plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["pdf.fonttype"] = 42  # TrueType: text editable in Illustrator


def fig1d():
    with open(TABLE, newline="") as f:
        rows = [r for r in csv.reader(f)][1:]
    labels = [r[1].replace(" (", "\n(") for r in rows]
    cargo = [float(r[4]) for r in rows]
    bg = [float(r[7]) for r in rows]
    print("Fig. 1d bars:")
    for lab, c, b in zip(labels, cargo, bg):
        print(f"  {lab.replace(chr(10), ' ')}: cargo {c*100:.2f}%  background {b*100:.2f}%")

    x = range(len(labels))
    w, off = 0.30, 0.19
    fig, ax = plt.subplots(figsize=(7, 4.4))
    ax.bar([i - off for i in x], cargo, w, label="Cargos",
           color="#D25D5D", edgecolor="black", linewidth=0.8)
    ax.bar([i + off for i in x], bg, w, label="Background",
           color="#EEEEEE", edgecolor="black", linewidth=0.8)
    ax.set_ylim(0, 0.08)
    ax.yaxis.set_major_locator(FixedLocator([0, 0.02, 0.04, 0.06, 0.08]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v * 100:g}%"))
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Fraction of proteins with database hit")
    ax.legend(frameon=False, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    base = os.path.join(HERE, "Fig1d_enrichment_bars")
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=300)
    plt.close(fig)
    print("figure -> " + base + ".pdf")


if __name__ == "__main__":
    fig1d()
