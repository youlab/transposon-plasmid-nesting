#!/usr/bin/env python3
"""Plot transposase-annotation density for plasmids and chromosomes.

Input: ``data/01_transposase_counts.csv`` in this package. Output: editable PDF
and PNG in ``results/figures`` by default.

Example: ``python code/fig1a.py --output-dir results/figures`` (from anywhere).
"""

import argparse
import csv
import warnings
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = PACKAGE_ROOT / "data"
DEFAULT_OUTPUT_DIR = PACKAGE_ROOT / "results" / "figures"
FS_BASE = 9
FS_LABEL = 10
BAR_COLOR = "#7F7F7F"


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


def fig1a(input_dir=DEFAULT_INPUT_DIR, output_dir=DEFAULT_OUTPUT_DIR):
    plt = _configure_font()
    input_path = Path(input_dir) / "01_transposase_counts.csv"
    output_dir = Path(output_dir)
    aggregate = {}
    with input_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["SeqType"] not in ("plasmid", "chromosome"):
                continue
            totals = aggregate.setdefault(row["SeqType"], [0, 0])
            totals[0] += int(row["transposase_count"])
            totals[1] += int(row["SeqLength"])
    densities = {kind: hits / length * 1e6 for kind, (hits, length) in aggregate.items()}
    print("1a densities per Mbp:", {kind: round(value, 3) for kind, value in densities.items()})

    labels = ["chromosome", "plasmid"]
    values = [densities[label] for label in labels]
    fig, ax = plt.subplots(figsize=(3.6, 1.5))
    ax.barh(labels, values, color=BAR_COLOR, height=0.62)
    ax.set_xlim(0, 60)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_xlabel("Transposase-related annotations per Mbp")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3)
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    base = output_dir / "Fig1a_density_per_Mbp"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"figure -> {base}.pdf")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args(argv)
    fig1a(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
