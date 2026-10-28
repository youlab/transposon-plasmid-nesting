#!/usr/bin/env python3
"""Purpose: create QC and relaxed-enrichment figures from existing tables.
Inputs: 05_transposase_count_comparison.csv and 11_enrichment_relaxed.csv (from work/processed/ or data/).
Outputs: PDF and PNG figures under results/figures/.
Usage: python code/pipeline/scripts/16_final_figures.py
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
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))))
DATA_DIR = os.path.abspath(os.environ.get(
    "NEE_DATA_DIR", os.path.join(ROOT, "work", "raw")))
RESULTS_DIR = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "work", "processed")))
FIGURES_DIR = os.path.join(ROOT, "results", "figures")
STATIC_DIR = os.path.join(ROOT, "data")


def input_path(name, generated=False):
    """Prefer current work products, then raw inputs or bundled static tables."""
    directories = (RESULTS_DIR, DATA_DIR, STATIC_DIR) if generated else (DATA_DIR, STATIC_DIR, RESULTS_DIR)
    candidates = [os.path.join(directory, name) for directory in directories]
    return next((path for path in candidates if os.path.exists(path)), candidates[0])

RES = RESULTS_DIR

font_file = os.environ.get("NEE_FONT_FILE")
if font_file:
    if not os.path.isfile(font_file):
        raise FileNotFoundError(f"NEE_FONT_FILE does not exist: {font_file}")
    font_manager.fontManager.addfont(font_file)
    default_family = font_manager.FontProperties(fname=font_file).get_name()
else:
    default_family = os.environ.get("NEE_FONT_FAMILY", "Arial")
plt.rcParams["font.family"] = default_family
plt.rcParams["axes.unicode_minus"] = False


# ---------- A) Step 05 scatter ----------

def fig_05():
    xs, ys = [], []
    with open(input_path("05_transposase_count_comparison.csv", generated=True)) as f:
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
    os.makedirs(FIGURES_DIR, exist_ok=True)
    base = os.path.join(FIGURES_DIR, "05_transposase_count_comparison")
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=150)
    plt.close(fig)
    print("figure -> " + base + ".pdf")


# ---------- B) Step 11 enrichment bars ----------

def fig_11():
    rows = [r for r in csv.reader(
        open(input_path("11_enrichment_relaxed.csv", generated=True)))][1:]
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
    base = os.path.join(FIGURES_DIR, "11_enrichment_relaxed")
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=150)
    plt.close(fig)
    print("figure -> " + base + ".pdf")


if __name__ == "__main__":
    fig_05()
    fig_11()
