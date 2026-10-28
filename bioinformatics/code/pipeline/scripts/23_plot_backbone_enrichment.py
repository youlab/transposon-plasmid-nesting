#!/usr/bin/env python3
"""Plot current PlasAnn backbone-filtered hit fractions.

Inputs: ``data/21_backbone_enrichment.csv``,
``data/20_cds_counts_after_filtering.csv``, and
``code/pipeline/config/analysis.json`` in this package. Output: relaxed and
strict Extended Data Figure 1 panels (PDF, SVG and 600-dpi PNG) in
``results/figures`` by default.

Example: ``python code/extended_fig1.py --output-dir results/figures``
(from anywhere).
"""

import argparse
import csv
import importlib.util
import json
import math
import sys
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT_DIR = PACKAGE_ROOT / "data"
DEFAULT_OUTPUT_DIR = PACKAGE_ROOT / "results" / "figures"
CONFIG_PATH = PACKAGE_ROOT / "code" / "pipeline" / "config" / "analysis.json"
HELPERS_PATH = Path(__file__).with_name("22_backbone_plot_helpers.py")
HELPERS_SPEC = importlib.util.spec_from_file_location("backbone_plot_helpers", HELPERS_PATH)
if HELPERS_SPEC is None or HELPERS_SPEC.loader is None:
    raise ImportError(f"Could not load plotting helpers from {HELPERS_PATH}")
HELPERS = importlib.util.module_from_spec(HELPERS_SPEC)
sys.modules[HELPERS_SPEC.name] = HELPERS
HELPERS_SPEC.loader.exec_module(HELPERS)

STAT_FIELDS = (
    "mode", "database", "label", "cargo_hit_cds", "cargo_other_retained_cds",
    "cargo_retained_cds_total", "cargo_hit_fraction", "background_hit_cds",
    "background_other_retained_cds", "background_retained_cds_total",
    "background_hit_fraction", "odds_ratio_cargo_vs_background",
    "p_value_two_sided", "log10_p", "q_value_BH", "log10_q", "n_tests_BH",
)
COUNT_FIELDS = ("mode", "group", "database", "hit_cds", "cds_total", "hit_fraction", "hit_carriers")
MODES = HELPERS.MODES
GROUPS = HELPERS.GROUPS
DATABASE_LABELS = dict(HELPERS.DATABASES)


def _read_rows(path, required_fields):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not set(required_fields) <= set(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: {', '.join(required_fields)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} contains no data rows")
    return rows


def _integer(row, name, *, minimum=0):
    try:
        value = int(row[name])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid integer {name}={row.get(name)!r}") from exc
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}, got {value}")
    return value


def _number(row, name):
    try:
        value = float(row[name])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid number {name}={row.get(name)!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value}")
    return value


def _config_categories(config_path=CONFIG_PATH):
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise FileNotFoundError(f"Missing package analysis config: {config_path}") from exc
    categories = config.get("backbone_categories")
    variant = config.get("variant")
    if not isinstance(categories, list) or not categories or any(db not in DATABASE_LABELS for db in categories):
        raise ValueError(f"Invalid backbone_categories in {config_path}")
    expected_count = {"WithTA": 4, "WithoutTA": 3}.get(variant)
    if expected_count is None or len(categories) != expected_count or len(set(categories)) != len(categories):
        raise ValueError(f"Unexpected backbone category list for variant {variant!r} in {config_path}")
    return categories


def _validate_inputs(input_dir, categories):
    input_dir = Path(input_dir)
    count_rows = _read_rows(input_dir / "20_cds_counts_after_filtering.csv", COUNT_FIELDS)
    count_map = {}
    for row in count_rows:
        key = (row["mode"], row["group"], row["database"])
        if key in count_map:
            raise ValueError(f"Duplicate CDS count row: {key}")
        if key[0] not in MODES or key[1] not in GROUPS or key[2] not in categories:
            raise ValueError(f"Unexpected mode/group/database: {key}")
        hit = _integer(row, "hit_cds")
        total = _integer(row, "cds_total", minimum=1)
        carriers = _integer(row, "hit_carriers")
        fraction = _number(row, "hit_fraction")
        if hit > total or carriers > hit:
            raise ValueError(f"Hit/carrier counts exceed CDS counts for {key}")
        if not 0.0 <= fraction <= 1.0 or not math.isclose(
            fraction, hit / total, rel_tol=1e-8, abs_tol=1e-12
        ):
            raise ValueError(f"Fraction does not agree with hit/denominator counts for {key}")
        count_map[key] = {**row, "_fraction": fraction}

    expected_count_keys = {
        (mode, group, database)
        for mode in MODES for group in GROUPS for database in categories
    }
    if set(count_map) != expected_count_keys:
        missing = sorted(expected_count_keys - set(count_map))
        extra = sorted(set(count_map) - expected_count_keys)
        raise ValueError(f"CDS counts do not cover configured combinations; missing={missing}, extra={extra}")

    for mode in MODES:
        for group in GROUPS:
            denominators = {count_map[(mode, group, db)]["cds_total"] for db in categories}
            if len(denominators) != 1:
                raise ValueError(f"Inconsistent retained-CDS denominators for {mode}/{group}")

    stats_rows = _read_rows(input_dir / "21_backbone_enrichment.csv", STAT_FIELDS)
    stats_map = {}
    for row in stats_rows:
        key = (row["mode"], row["database"])
        if key in stats_map:
            raise ValueError(f"Duplicate enrichment row: {key}")
        if key[0] not in MODES or key[1] not in categories:
            raise ValueError(f"Unexpected enrichment mode/database: {key}")
        for group, hit_field, other_field, total_field, fraction_field in (
            ("cargo", "cargo_hit_cds", "cargo_other_retained_cds", "cargo_retained_cds_total", "cargo_hit_fraction"),
            ("background", "background_hit_cds", "background_other_retained_cds", "background_retained_cds_total", "background_hit_fraction"),
        ):
            count = count_map[(key[0], group, key[1])]
            hit = _integer(row, hit_field)
            other = _integer(row, other_field)
            total = _integer(row, total_field, minimum=1)
            fraction = _number(row, fraction_field)
            if hit != _integer(count, "hit_cds") or total != _integer(count, "cds_total"):
                raise ValueError(f"Enrichment counts disagree with CDS count table for {key}/{group}")
            if other != total - hit:
                raise ValueError(f"Enrichment complement count is inconsistent for {key}/{group}")
            if not math.isclose(fraction, count["_fraction"], rel_tol=1e-8, abs_tol=1e-12):
                raise ValueError(f"Enrichment fraction disagrees with CDS count table for {key}/{group}")
        for probability in ("p_value_two_sided", "q_value_BH"):
            value = _number(row, probability)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{probability} must be between 0 and 1 for {key}")
        for probability_log in ("log10_p", "log10_q"):
            if _number(row, probability_log) > 0.0:
                raise ValueError(f"{probability_log} must be non-positive for {key}")
        stats_map[key] = row

    expected_stats_keys = {(mode, db) for mode in MODES for db in categories}
    if set(stats_map) != expected_stats_keys:
        missing = sorted(expected_stats_keys - set(stats_map))
        extra = sorted(set(stats_map) - expected_stats_keys)
        raise ValueError(f"Enrichment rows do not cover configured Fisher families; missing={missing}, extra={extra}")

    n_tests = len(categories)
    for mode in MODES:
        if any(_integer(stats_map[(mode, db)], "n_tests_BH", minimum=1) != n_tests for db in categories):
            raise ValueError(f"BH test count must be {n_tests} for {mode}")
        ordered = sorted(categories, key=lambda db: _number(stats_map[(mode, db)], "log10_p"))
        running = 0.0
        expected_log10_q = {}
        for rank in range(n_tests, 0, -1):
            database = ordered[rank - 1]
            log10_p = _number(stats_map[(mode, database)], "log10_p")
            running = min(running, log10_p + math.log10(n_tests / rank))
            expected_log10_q[database] = min(0.0, running)
        for database in categories:
            observed = _number(stats_map[(mode, database)], "log10_q")
            if not math.isclose(observed, expected_log10_q[database], rel_tol=0.0, abs_tol=1e-6):
                raise ValueError(f"BH log10_q does not match the configured {n_tests}-test family for {mode}/{database}")

    fraction_map = {
        key: row["_fraction"] for key, row in count_map.items()
    }
    return fraction_map


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--config", type=Path, default=CONFIG_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args(argv)

    categories = _config_categories(args.config)
    HELPERS.PLOTTED_DATABASES = tuple((db, DATABASE_LABELS[db]) for db in categories)
    fractions = _validate_inputs(args.input_dir, categories)
    max_fraction = max(fractions.values())
    y_max, y_ticks = HELPERS._nice_percent_axis(max_fraction)
    HELPERS._configure_matplotlib()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for mode in MODES:
        base = HELPERS._plot_mode(mode, fractions, y_max, y_ticks, args.output_dir)
        print(f"figure -> {base}.pdf")
    print(f"Generated current relaxed and strict figures in {args.output_dir}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        sys.exit(f"ERROR: {exc}")
