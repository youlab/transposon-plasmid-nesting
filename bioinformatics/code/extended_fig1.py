#!/usr/bin/env python3
"""Run the shared Extended Data Figure 1 plotting wrapper.

Inputs: package ``data/20_cds_counts_after_filtering.csv`` and
``data/21_backbone_enrichment.csv`` plus ``code/pipeline/config/analysis.json``.
Output: current relaxed and strict PDF/SVG/PNG panels in
``results/figures`` by default.

Example: ``python code/extended_fig1.py --output-dir results/figures``
(from anywhere).
"""

import importlib.util
import sys
from pathlib import Path


WRAPPER_PATH = Path(__file__).resolve().parent / "pipeline" / "scripts" / "23_plot_backbone_enrichment.py"
WRAPPER_SPEC = importlib.util.spec_from_file_location("plot_backbone_enrichment", WRAPPER_PATH)
if WRAPPER_SPEC is None or WRAPPER_SPEC.loader is None:
    raise ImportError(f"Could not load plotting wrapper from {WRAPPER_PATH}")
WRAPPER = importlib.util.module_from_spec(WRAPPER_SPEC)
sys.modules[WRAPPER_SPEC.name] = WRAPPER
WRAPPER_SPEC.loader.exec_module(WRAPPER)


def main(argv=None):
    return WRAPPER.main(argv)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        sys.exit(f"ERROR: {exc}")
