#!/usr/bin/env python3
"""
16_final_figures.py
Publication-styled figures:
  A) Step 05 scatter — log ticks only 10^0..10^3 (no minor ticks);
     single color, no outlier highlighting; axis titles
     "Transposase counts based on CDS annotation" /
     "IS element counts based on ISEScan"; no title,
     Spearman rho computed from the input table as in-plot text.
  B) Step 11 enrichment bars — y axis 0-8% with % ticks; no q labels
     (a sentence in text covers it); colors #D25D5D / #EEEEEE with black
     edges; wider bar separation; legend "Cargos" / "Background".
Reads existing result tables; overwrites the two figure files (+ .png).
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import LogLocator, NullLocator, FixedLocator, FuncFormatter
from scipy.stats import spearmanr

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
RES = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))

font_file = os.environ.get("NEE_FONT_FILE")
if font_file:
    if not os.path.isfile(font_file):
        raise FileNotFoundError(f"NEE_FONT_FILE does not exist: {font_file}")
    font_manager.fontManager.addfont(font_file)
    default_family = font_manager.FontProperties(fname=font_file).get_name()
else:
    default_family = os.environ.get("NEE_FONT_FAMILY", "DejaVu Sans")
plt.rcParams["font.family"] = default_family
plt.rcParams["axes.unicode_minus"] = False


# ---------- A) Step 05 scatter ----------

def fig_05():
    xs, ys = [], []
    with open(os.path.join(RES, "4b_Transposase_count_comparison.csv")) as f:
        for r in csv.reader(f):
            if r[0] == "accession" or r[2] == "":
                continue
            rohan, iscan = int(r[1]), int(r[2])
            if iscan >= 1:
                xs.append(rohan)
                ys.append(iscan)
    fig, ax = plt.subplots(figsize=(4.8, 4.8))
    ax.scatter(xs, ys, s=4, alpha=0.15, edgecolors="none", color="#4C72B0")
    lim = [0.8, max(max(xs), max(ys)) * 1.2]
    ax.plot(lim, lim, color="grey", lw=0.8, ls="--", label="y = x")
    # log-log linear fit (power law y = a * x^b) on the displayed points
    import numpy as np
    lx, ly = np.log10(xs), np.log10(ys)
    b, a = np.polyfit(lx, ly, 1)
    xfit = np.logspace(np.log10(lim[0]), np.log10(lim[1]), 100)
    ax.plot(xfit, 10 ** (a + b * np.log10(xfit)), color="#C44E52", lw=1.2,
            label="fitted line")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(LogLocator(base=10, numticks=4))
        axis.set_minor_locator(NullLocator())
    ax.set_xlabel("Transposase counts based on CDS annotation")
    ax.set_ylabel("IS element counts based on ISEScan")
    rho, _ = spearmanr(xs, ys)
    ax.text(0.05, 0.95, f"Spearman $\\rho$ = {rho:.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    base = os.path.join(RES, "4b_Transposase_count_comparison")
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=150)
    plt.close(fig)
    print("figure -> " + base + ".pdf")


# ---------- B) Step 11 enrichment bars ----------

def fig_11():
    rows = [r for r in csv.reader(
        open(os.path.join(RES, "7_Enrichment_table.csv")))][1:]
    labels = [r[1].replace(" (", "\n(") for r in rows]
    cargo = [float(r[4]) for r in rows]
    bg = [float(r[7]) for r in rows]

    x = range(len(labels))
    w, off = 0.30, 0.19
    fig, ax = plt.subplots(figsize=(7, 4.4))
    ax.bar([i - off for i in x], cargo, w, label="Cargos",
           color="#D25D5D", edgecolor="black", linewidth=0.8)
    ax.bar([i + off for i in x], bg, w, label="Background",
           color="#EEEEEE", edgecolor="black", linewidth=0.8)
    ax.set_ylim(0, 0.08)
    ax.yaxis.set_major_locator(FixedLocator([0, 0.02, 0.04, 0.06, 0.08]))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, _: f"{v * 100:g}%"))
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Fraction of proteins with database hit")
    ax.legend(frameon=False, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    base = os.path.join(RES, "7_Enrichment_figure")
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=150)
    plt.close(fig)
    print("figure -> " + base + ".pdf")


if __name__ == "__main__":
    fig_05()
    fig_11()
