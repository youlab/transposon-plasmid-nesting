#!/usr/bin/env python3
"""
fig1a.py — regenerate Figure 1a as an editable vector PDF.

Figure 1a shows transposase-annotation density per Mbp for plasmids and
chromosomes, computed from transposase-counts.csv.

Run:  python3 fig1a.py   (reads ../data/transposase-counts.csv; writes figures next to this script)
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE = os.path.dirname(os.path.abspath(__file__))
TC_CSV = os.path.join(HERE, "..", "data", "transposase-counts.csv")
OUT = HERE
os.makedirs(OUT, exist_ok=True)

# ---- adjustable style ----
FS_BASE = 9    # tick labels / in-plot text
FS_LABEL = 10  # axis titles
BAR_COLOR = "#7F7F7F"

FONT_FILE = "/System/Library/Fonts/Supplemental/Arial.ttf"
if os.path.exists(FONT_FILE):
    font_manager.fontManager.addfont(FONT_FILE)
    plt.rcParams["font.family"] = font_manager.FontProperties(fname=FONT_FILE).get_name()
else:
    plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = FS_BASE
plt.rcParams["axes.labelsize"] = FS_LABEL
plt.rcParams["pdf.fonttype"] = 42  # TrueType: text editable in Illustrator
plt.rcParams["svg.fonttype"] = "none"


def save(fig, name):
    base = os.path.join(OUT, name)
    fig.savefig(base + ".pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("figure -> " + base + ".pdf")


# ---------- Fig 1a ----------
def fig1a():
    agg = {}
    with open(TC_CSV, newline="") as f:
        for row in csv.DictReader(f):
            if row["SeqType"] not in ("plasmid", "chromosome"):
                continue
            d = agg.setdefault(row["SeqType"], [0, 0])
            d[0] += int(row["transposase_count"])
            d[1] += int(row["SeqLength"])
    dens = {k: v[0] / v[1] * 1e6 for k, v in agg.items()}
    print("1a densities per Mbp:", {k: round(v, 3) for k, v in dens.items()})

    labels = ["chromosome", "plasmid"]  # bottom -> top
    vals = [dens[k] for k in labels]
    fig, ax = plt.subplots(figsize=(3.6, 1.5))
    ax.barh(labels, vals, color=BAR_COLOR, height=0.62)
    ax.set_xlim(0, 60)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_xlabel("Transposase-related annotations per Mbp")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3)
    fig.tight_layout()
    save(fig, "Fig1a_density_per_Mbp")



if __name__ == "__main__":
    fig1a()
