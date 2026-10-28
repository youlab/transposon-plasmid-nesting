#!/usr/bin/env python3
"""Plot functional-category enrichment for cargo proteins and background.

Input: ``data/11_enrichment_relaxed.csv`` in this package. Fractions are
computed from integer hit and denominator counts. Output: PDF and PNG in
``results/figures`` by default.

Example: ``python code/fig1d.py --output-dir results/figures`` (from anywhere).
"""

import argparse
import csv
import math
import warnings
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = PACKAGE_ROOT / "data"
DEFAULT_OUTPUT_DIR = PACKAGE_ROOT / "results" / "figures"


def _configure_font():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    try:
        font_manager.findfont("Arial", fallback_to_default=False)
        plt.rcParams["font.family"] = "Arial"
    except ValueError:
        warnings.warn("Arial was not found; using matplotlib's default sans-serif font.", RuntimeWarning)
        plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams.update({"axes.unicode_minus": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
    return plt


def _fraction(row, hit_field, total_field):
    hits = int(row[hit_field])
    total = int(row[total_field])
    if total <= 0 or hits < 0 or hits > total:
        raise ValueError(f"Invalid hit/denominator counts: {hit_field}={hits}, {total_field}={total}")
    return hits / total


def fig1d(input_dir=DEFAULT_INPUT_DIR, output_dir=DEFAULT_OUTPUT_DIR):
    plt = _configure_font()
    from matplotlib.ticker import FixedLocator, FuncFormatter

    table = Path(input_dir) / "11_enrichment_relaxed.csv"
    with table.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"label", "cargo_hits", "cargo_total", "bg_hits", "bg_total"}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"{table} must contain columns: {', '.join(sorted(required))}")
    labels = [row["label"].replace(" (", "\n(") for row in rows]
    cargo = [_fraction(row, "cargo_hits", "cargo_total") for row in rows]
    background = [_fraction(row, "bg_hits", "bg_total") for row in rows]
    print("Fig. 1d bars:")
    for label, cargo_fraction, background_fraction in zip(labels, cargo, background):
        print(f"  {label.replace(chr(10), ' ')}: cargo {cargo_fraction*100:.2f}%  background {background_fraction*100:.2f}%")

    x = range(len(labels))
    width, offset = 0.30, 0.19
    fig, ax = plt.subplots(figsize=(7, 4.4))
    ax.bar([i - offset for i in x], cargo, width, label="Cargos",
           color="#D25D5D", edgecolor="black", linewidth=0.8)
    ax.bar([i + offset for i in x], background, width, label="Background",
           color="#EEEEEE", edgecolor="black", linewidth=0.8)
    upper = max(0.08, math.ceil(max(cargo + background) * 1.05 / 0.02) * 0.02)
    ticks = [index * 0.02 for index in range(int(round(upper / 0.02)) + 1)]
    ax.set_ylim(0, upper)
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value * 100:g}%"))
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Fraction of proteins with database hit")
    ax.legend(frameon=False, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    base = output_dir / "Fig1d_enrichment_bars"
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".png"), dpi=300)
    plt.close(fig)
    print(f"figure -> {base}.pdf")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args(argv)
    fig1d(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
