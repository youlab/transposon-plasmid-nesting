#!/usr/bin/env python3
"""
fig1b.py — data-driven Fig. 1b from the two PUBLIC source tables.

Inputs (the plain CSVs shipped in this package, resolved relative to this file):
  - ../data/transposase-counts.csv    (63,425 replicons)
  - ../data/PIRA-PCN-estimates.csv    (16,861 rows; Nat Commun 16:6023 (2025) supp.)

Join on SeqID (the replicon accession, e.g. NZ_CP107404.1), NOT on
AnnotationAccession — in both tables AnnotationAccession is the GCF *assembly*
accession shared by every replicon of a genome, so joining on it scrambles the
hit lookup. SeqType == plasmid on both sides. PCN bins: <1, [1,10),
[10,100), [100,1000); plasmids with PIRACopyNumber >= 1000 (2 of them) fall
outside the plotted bins and are counted separately. A "hit" = transposase_count >= 1.
95% CIs: same normal approximation as Rohan's transposase-plasmid-figure.R
(+/-1.96*se, truncated to [0,1]).

Expected values for the included inputs: joined plasmids 11,010,
bin totals 2322/6396/2124/166 (+2 with PCN>=1000), fractions 82.8/65.6/7.8/3.6%.

Layout identical to fig1a.py (compressed vertically, editable TrueType text).

Run:  python3 fig1b.py   (writes figures next to this script)
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
PIRA_CSV = os.path.join(HERE, "..", "data", "PIRA-PCN-estimates.csv")
OUT = HERE

FS_BASE = 9
FS_LABEL = 10
TITLE_1B = "Fraction of plasmids containing at least one\ntransposase-related annotation"

FONT_FILE = "/System/Library/Fonts/Supplemental/Arial.ttf"
if os.path.exists(FONT_FILE):
    font_manager.fontManager.addfont(FONT_FILE)
    plt.rcParams["font.family"] = font_manager.FontProperties(fname=FONT_FILE).get_name()
else:
    plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = FS_BASE
plt.rcParams["axes.labelsize"] = FS_LABEL
plt.rcParams["pdf.fonttype"] = 42


def load_hits():
    hit = {}
    with open(TC_CSV, newline="") as f:
        for row in csv.DictReader(f):
            hit[row["SeqID"]] = int(row["transposase_count"]) >= 1
    return hit


def bin_pcn(pcn):
    if pcn < 1:
        return "<1"
    if pcn < 10:
        return "1-10"
    if pcn < 100:
        return "10-100"
    if pcn < 1000:
        return "100-1000"
    return ">=1000"


def compute_rows():
    hit = load_hits()
    bins = ["<1", "1-10", "10-100", "100-1000"]
    tot = {b: 0 for b in bins}
    pos = {b: 0 for b in bins}
    n_join = n_ge1000 = 0
    with open(PIRA_CSV, newline="") as f:
        for row in csv.DictReader(f):
            if row["SeqType"] != "plasmid":
                continue
            acc = row["SeqID"]
            if acc not in hit:
                continue
            n_join += 1
            b = bin_pcn(float(row["PIRACopyNumber"]))
            if b == ">=1000":
                n_ge1000 += 1
                continue
            tot[b] += 1
            pos[b] += hit[acc]
    print(f"joined plasmids: {n_join} (excluded PCN>=1000: {n_ge1000})")
    return [(b, pos[b], tot[b]) for b in bins]


def fig1b(rows):
    ys, ps, los, his, texts = [], [], [], [], []
    for i, (lab, h, t) in enumerate(rows):
        p = h / t
        se = math.sqrt(p * (1 - p) / t)
        lo, hi = max(0.0, p - 1.96 * se), min(1.0, p + 1.96 * se)
        ys.append(i); ps.append(p); los.append(p - lo); his.append(hi - p)
        texts.append(f"{h}/{t}")
        print(f"1b {lab}: {h}/{t} = {p*100:.2f}%  CI [{lo:.4f}, {hi:.4f}]")

    fig, ax = plt.subplots(figsize=(5.2, 1.7))  # compressed vertically
    ax.errorbar(ps, ys, xerr=[los, his], fmt="o", ms=3.5, color="black",
                elinewidth=0.7, capsize=0, zorder=3)
    for y, p, hi, t in zip(ys, ps, his, texts):
        ax.text(p + hi + 0.03, y, t, va="center", ha="left", fontsize=FS_BASE)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_ylim(-0.55, 3.55)
    ax.set_xlim(0, 1.06)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xlabel("proportion of plasmids")
    ax.set_ylabel("PCN bin")
    if TITLE_1B:
        ax.set_title(TITLE_1B, fontsize=FS_LABEL, loc="left", pad=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3)
    fig.tight_layout()
    base = os.path.join(OUT, "Fig1b_PCN_bin_fractions_publicPCN")
    fig.savefig(base + ".pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("figure -> " + base + ".pdf")


if __name__ == "__main__":
    fig1b(compute_rows())
