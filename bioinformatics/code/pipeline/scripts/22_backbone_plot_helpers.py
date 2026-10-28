#!/usr/bin/env python3
"""Shared plotting helpers for Extended Data Figure 1.

Purpose: render the current relaxed or strict backbone-filtered hit fractions.
Input: validated statistics/count rows prepared by 23_plot_backbone_enrichment.py.
Output: PDF, SVG and 600-dpi PNG in the caller-selected output directory.
Example: ``python code/extended_fig1.py --help`` (from anywhere).
"""

import warnings
from pathlib import Path

MODES = ("relaxed", "strict")
GROUPS = ("cargo", "background")
DATABASES = (
    ("card", "Antibiotic resistance\n(CARD)"),
    ("vfdb", "Virulence\n(VFDB)"),
    ("bacmet", "Metal/biocide resistance\n(BacMet)"),
    ("tadb", "Toxin–antitoxin systems\n(TADB)"),
)
PLOTTED_DATABASES = DATABASES


def _load_plotting():
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.ticker import FixedLocator, FuncFormatter
    return mpl, plt, font_manager, FixedLocator, FuncFormatter


def _nice_percent_axis(max_fraction):
    """Return a shared upper fraction and readable percentage ticks."""
    import math

    target = max(2.0, max_fraction * 100.0 * 1.12)
    raw_step = target / 5.0
    scale = 10 ** math.floor(math.log(raw_step, 10.0))
    step = next((factor * scale for factor in (1, 2, 2.5, 5, 10) if factor * scale >= raw_step), 10 * scale)
    upper_pct = math.ceil(target / step) * step
    tick_count = int(round(upper_pct / step))
    ticks = [index * step / 100.0 for index in range(tick_count + 1)]
    return upper_pct / 100.0, ticks


def _configure_matplotlib():
    mpl, _, font_manager, _, _ = _load_plotting()
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
        font = "Arial"
    except ValueError:
        warnings.warn("Arial was not found; using matplotlib's default sans-serif font.", RuntimeWarning)
        font = "DejaVu Sans"
    mpl.rcParams.update({
        "font.family": font,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
        "font.size": 9,
        "savefig.facecolor": "white",
    })


def _plot_mode(mode, hit_map, y_max, y_ticks, out_dir):
    _, plt, _, FixedLocator, FuncFormatter = _load_plotting()
    out_dir = Path(out_dir)
    x = range(len(PLOTTED_DATABASES))
    width, offset = 0.30, 0.19
    fig, ax = plt.subplots(figsize=(7, 4.4))
    for group, shift, color, label in (
        ("cargo", -offset, "#D25D5D", "Cargos"),
        ("background", offset, "#EEEEEE", "Background"),
    ):
        fractions = [hit_map[(mode, group, database)] for database, _ in PLOTTED_DATABASES]
        ax.bar(
            [index + shift for index in x], fractions, width,
            label=label, color=color, edgecolor="black", linewidth=0.8,
        )

    ax.set_ylim(0, y_max)
    ax.yaxis.set_major_locator(FixedLocator(y_ticks))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value * 100:g}%"))
    ax.set_xticks(list(x), [label for _, label in PLOTTED_DATABASES], fontsize=10)
    ax.set_ylabel("Fraction", fontsize=10)
    ax.tick_params(axis="both", labelsize=9)
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    out_dir.mkdir(parents=True, exist_ok=True)
    base = out_dir / f"ExtendedFig1_{mode}"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    return base
