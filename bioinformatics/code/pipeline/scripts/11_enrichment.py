#!/usr/bin/env python3
"""
11_enrichment.py
Enrichment of functional categories in transposon cargo vs background
(proteins on the same IS-carrying plasmids), with Fisher exact tests and
Benjamini-Hochberg FDR across categories, plus the summary figure.

Categories: CARD = antibiotic resistance, VFDB = virulence,
            BacMet = metal/biocide resistance, TADB = toxin-antitoxin.

Counts are copy-weighted: each unique protein's hit status maps back to all
its copies via 6_{Cargo,Background}_unique_map.csv (protein md5 dedup was
only a compute saver and has no statistical effect).

Usage:  python scripts/11_enrichment.py [min_pident] [min_qcov] [tag]
        default 90 60 (main); e.g. python scripts/11_enrichment.py 80 70 _p80q70

Input:  2_Result/6_{Cargo,Background}_vs_{card,vfdb,bacmet,tadb}.tsv
        2_Result/6_{Cargo,Background}_unique_map.csv
Output: 2_Result/7_Enrichment_table{tag}.csv
        2_Result/7_Enrichment_figure{tag}.pdf/.png
        2_Result/7_Summary{tag}.txt
"""
import csv
import os
import sys

from scipy.stats import fisher_exact
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
RES = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))
MIN_PIDENT = float(sys.argv[1]) if len(sys.argv) > 1 else 90.0
MIN_QCOV = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
TAG = sys.argv[3] if len(sys.argv) > 3 else ""

CATS = {"card": "Antibiotic resistance\n(CARD)",
        "vfdb": "Virulence\n(VFDB)",
        "bacmet": "Metal/biocide resistance\n(BacMet)",
        "tadb": "Toxin-antitoxin\n(TADB)"}


def load_map(tag):
    m = {}
    with open(os.path.join(RES, f"6_{tag}_unique_map.csv")) as f:
        for row in csv.DictReader(f):
            m[row["md5"]] = int(row["n_copies"])
    return m


def load_hits(setname, cat, umap):
    """Weighted hit count for one set x category."""
    n = 0
    with open(os.path.join(RES, f"6_{setname}_vs_{cat}.tsv"),
              errors="replace") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 7:
                continue
            if float(row[2]) >= MIN_PIDENT and float(row[4]) >= MIN_QCOV:
                n += umap.get(row[0], 1)
    return n


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted q-values (same order as input)."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        val = min(prev, pvals[i] * m / (m - rank + 1))
        q[i] = val
        prev = val
    return q


def main():
    maps = {s: load_map(s) for s in ("Cargo", "Background")}
    totals = {s: sum(maps[s].values()) for s in maps}

    rows = []
    for cat, label in CATS.items():
        a = load_hits("Cargo", cat, maps["Cargo"])
        c = load_hits("Background", cat, maps["Background"])
        b = totals["Cargo"] - a
        d = totals["Background"] - c
        odds, p = fisher_exact([[a, b], [c, d]])
        fc = a / totals["Cargo"]
        fb = c / totals["Background"]
        rows.append([cat, label.replace("\n", " "), a, totals["Cargo"],
                     f"{fc:.4f}", c, totals["Background"], f"{fb:.4f}",
                     f"{odds:.3f}", f"{p:.3e}"])

    qvals = bh_fdr([float(r[9]) for r in rows])
    for r, q in zip(rows, qvals):
        r.append(f"{q:.3e}")

    out_csv = os.path.join(RES, f"7_Enrichment_table{TAG}.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["category", "label", "cargo_hits", "cargo_total",
                    "cargo_frac", "bg_hits", "bg_total", "bg_frac",
                    "odds_ratio", "p_fisher", "q_fdr_bh"])
        w.writerows(rows)

    lines = [
        f"# Step 7: cargo vs background enrichment "
        f"(pident>={MIN_PIDENT:g}, qcov>={MIN_QCOV:g})",
        f"cargo_total: {totals['Cargo']}, background_total: "
        f"{totals['Background']} (copy-weighted)",
    ]
    for r in rows:
        lines.append(f"{r[1]}: cargo {r[4]} vs bg {r[7]} "
                     f"| OR={r[8]} | P={r[9]} | q={r[10]}")
    summary = "\n".join(lines)
    print(summary)
    with open(os.path.join(RES, f"7_Summary{TAG}.txt"), "w") as f:
        f.write(summary + "\n")

    # figure
    labels = list(CATS.values())
    cargo_frac = [float(r[4]) for r in rows]
    bg_frac = [float(r[7]) for r in rows]
    x = range(len(labels))
    w = 0.35
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.bar([i - w / 2 for i in x], cargo_frac, w,
           label="Transposon cargo", color="#C44E52")
    ax.bar([i + w / 2 for i in x], bg_frac, w,
           label="Background (same plasmids)", color="#8C8C8C")
    ymax = max(max(cargo_frac), max(bg_frac))
    for i, r in enumerate(rows):
        ax.text(i, max(cargo_frac[i], bg_frac[i]) + ymax * 0.02,
                f"OR={r[8]}\nq={r[10]}", ha="center", fontsize=7.5)
    ax.set_ylim(0, ymax * 1.25)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Fraction of proteins with database hit")
    ax.set_title("Transposon cargo vs background "
                 f"(pident>={MIN_PIDENT:g}%, qcov>={MIN_QCOV:g}%)",
                 fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_pdf = os.path.join(RES, f"7_Enrichment_figure{TAG}.pdf")
    fig.savefig(out_pdf)
    out_png = os.path.join(RES, f"7_Enrichment_figure{TAG}.png")
    fig.savefig(out_png, dpi=150)
    print(f"figure -> {out_pdf} (+ .png)")


if __name__ == "__main__":
    main()
