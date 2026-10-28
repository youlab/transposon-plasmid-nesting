#!/usr/bin/env python3
"""Plot the fraction of plasmids with a transposase hit by PCN bin.

Inputs: ``data/01_transposase_counts.csv`` and
``data/01_plasmid_copy_number_estimates.csv`` in this package. Output: PDF and
PNG in ``results/figures`` by default.

Example: ``python code/fig1b.py --output-dir results/figures`` (from anywhere).
"""

import argparse
import csv
import math
import warnings
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = PACKAGE_ROOT / "data"
DEFAULT_OUTPUT_DIR = PACKAGE_ROOT / "results" / "figures"
FS_BASE = 9
FS_LABEL = 10
TITLE_1B = "Fraction of plasmids containing at least one\ntransposase-related annotation"


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
    plt.rcParams.update({
        "axes.unicode_minus": False,
        "font.size": FS_BASE,
        "axes.labelsize": FS_LABEL,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })
    return plt


def load_hits(input_dir=DEFAULT_INPUT_DIR):
    hit = {}
    path = Path(input_dir) / "01_transposase_counts.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
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


def compute_rows(input_dir=DEFAULT_INPUT_DIR):
    input_dir = Path(input_dir)
    hit = load_hits(input_dir)
    bins = ["<1", "1-10", "10-100", "100-1000"]
    totals = {label: 0 for label in bins}
    positives = {label: 0 for label in bins}
    n_joined = n_ge1000 = 0
    with (input_dir / "01_plasmid_copy_number_estimates.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        for row in csv.DictReader(handle):
            if row["SeqType"] != "plasmid":
                continue
            accession = row["SeqID"]
            if accession not in hit:
                continue
            n_joined += 1
            label = bin_pcn(float(row["PIRACopyNumber"]))
            if label == ">=1000":
                n_ge1000 += 1
                continue
            totals[label] += 1
            positives[label] += hit[accession]
    print(f"joined plasmids: {n_joined} (excluded PCN>=1000: {n_ge1000})")
    return [(label, positives[label], totals[label]) for label in bins]


def fig1b(rows, output_dir=DEFAULT_OUTPUT_DIR):
    plt = _configure_font()
    ys, proportions, lower_errors, upper_errors, texts = [], [], [], [], []
    for index, (label, hits, total) in enumerate(rows):
        proportion = hits / total
        standard_error = math.sqrt(proportion * (1 - proportion) / total)
        lower = max(0.0, proportion - 1.96 * standard_error)
        upper = min(1.0, proportion + 1.96 * standard_error)
        ys.append(index)
        proportions.append(proportion)
        lower_errors.append(proportion - lower)
        upper_errors.append(upper - proportion)
        texts.append(f"{hits}/{total}")
        print(f"1b {label}: {hits}/{total} = {proportion*100:.2f}%  CI [{lower:.4f}, {upper:.4f}]")

    fig, ax = plt.subplots(figsize=(5.2, 1.7))
    ax.errorbar(proportions, ys, xerr=[lower_errors, upper_errors], fmt="o", ms=3.5,
                color="black", elinewidth=0.7, capsize=0, zorder=3)
    for y, proportion, upper_error, text in zip(ys, proportions, upper_errors, texts):
        ax.text(proportion + upper_error + 0.03, y, text, va="center", ha="left", fontsize=FS_BASE)
    ax.set_yticks(ys)
    ax.set_yticklabels([row[0] for row in rows])
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
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    base = output_dir / "Fig1b_PCN_bin_fractions_publicPCN"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"figure -> {base}.pdf")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args(argv)
    fig1b(compute_rows(args.input_dir), args.output_dir)


if __name__ == "__main__":
    main()
