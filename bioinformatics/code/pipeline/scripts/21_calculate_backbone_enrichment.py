#!/usr/bin/env python3
"""Calculate pooled-CDS Fisher tests and within-mode BH correction.

Inputs are 20_cds_counts_before_filtering.csv and
20_cds_counts_after_filtering.csv under --input-dir (results/ by default,
falling back to the bundled data directory). Outputs are the legacy relaxed
and strict tables plus 21_backbone_enrichment.csv under results/.
"""
import argparse
import csv
import json
import math
from pathlib import Path

from scipy.special import logsumexp
from scipy.stats import fisher_exact, hypergeom

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = ROOT / "code" / "pipeline" / "config" / "analysis.json"
DATABASES = ("card", "vfdb", "bacmet", "tadb")
LABELS = {
    "card": "Antibiotic resistance (CARD)",
    "vfdb": "Virulence (VFDB)",
    "bacmet": "Metal/biocide resistance (BacMet)",
    "tadb": "Toxin-antitoxin (TADB)",
}
COUNT_FIELDS = ("mode", "group", "database", "hit_cds", "cds_total",
                "hit_fraction", "hit_carriers")


def _read_counts(path):
    rows = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != COUNT_FIELDS:
            raise ValueError(f"unexpected count columns in {path}: {reader.fieldnames}")
        for line_number, row in enumerate(reader, 2):
            mode, group, database = (row[key].strip().casefold()
                                     for key in ("mode", "group", "database"))
            key = (mode, group, database)
            if key in rows:
                raise ValueError(f"duplicate count row in {path}: {key}")
            if mode not in {"relaxed", "strict"} or group not in {"cargo", "background"}:
                raise ValueError(f"unexpected mode/group at {path}:{line_number}")
            if database not in DATABASES:
                raise ValueError(f"unexpected database at {path}:{line_number}: {database}")
            hit, total = int(row["hit_cds"]), int(row["cds_total"])
            if total < 0 or hit < 0 or hit > total:
                raise ValueError(f"invalid hit/total counts at {path}:{line_number}")
            rows[key] = {"hit": hit, "total": total}
    return rows


def _fisher_log_p(table):
    a, b, c, d = (int(value) for row in table for value in row)
    total, successes, draws = a + b + c + d, a + c, a + b
    low, high = max(0, draws - (total - successes)), min(draws, successes)
    support = range(low, high + 1)
    log_probs = hypergeom.logpmf(support, total, successes, draws)
    observed = float(hypergeom.logpmf(a, total, successes, draws))
    return min(0.0, float(logsumexp(log_probs[log_probs <= observed + 1e-12])))


def _bh_adjust_log_probabilities(log10_ps):
    order = sorted(range(len(log10_ps)), key=log10_ps.__getitem__)
    q_values, running, count = [0.0] * len(order), 0.0, len(order)
    for rank in range(count, 0, -1):
        index = order[rank - 1]
        running = min(running, log10_ps[index] + math.log(count / rank, 10.0))
        q_values[index] = min(0.0, running)
    return q_values


def _probability_text(log10_value):
    if log10_value == -math.inf:
        return "0.000000e+0"
    exponent = math.floor(log10_value)
    return f"{10 ** (log10_value - exponent):.6f}e{exponent:+d}"


def _check_complete(rows, databases):
    expected = {(mode, group, database)
                for mode in ("relaxed", "strict")
                for group in ("cargo", "background")
                for database in databases}
    if set(rows) != expected:
        raise ValueError(f"count table rows do not match expected mode/group/database set; missing={sorted(expected - set(rows))}, extra={sorted(set(rows) - expected)}")
    for mode in ("relaxed", "strict"):
        for group in ("cargo", "background"):
            totals = {rows[(mode, group, database)]["total"] for database in databases}
            if len(totals) != 1:
                raise ValueError(f"inconsistent {mode}/{group} denominators across databases: {totals}")


def _build_stats(rows, databases, mode):
    stats, log10_ps = [], []
    for database in databases:
        cargo = rows[(mode, "cargo", database)]
        background = rows[(mode, "background", database)]
        cargo_total, bg_total = cargo["total"], background["total"]
        if cargo_total < 1 or bg_total < 1:
            raise ValueError(f"empty cargo/background denominator for {mode}/{database}")
        cargo_hit, bg_hit = cargo["hit"], background["hit"]
        table = [[cargo_hit, cargo_total - cargo_hit],
                 [bg_hit, bg_total - bg_hit]]
        odds_ratio = float(fisher_exact(table, alternative="two-sided").statistic)
        log10_p = _fisher_log_p(table) / math.log(10.0)
        log10_ps.append(log10_p)
        stats.append({
            "database": database, "label": LABELS[database],
            "cargo_hit_cds": cargo_hit,
            "cargo_other_retained_cds": cargo_total - cargo_hit,
            "cargo_retained_cds_total": cargo_total,
            "cargo_hit_fraction": format(cargo_hit / cargo_total, ".17g"),
            "background_hit_cds": bg_hit,
            "background_other_retained_cds": bg_total - bg_hit,
            "background_retained_cds_total": bg_total,
            "background_hit_fraction": format(bg_hit / bg_total, ".17g"),
            "odds_ratio_cargo_vs_background": format(odds_ratio, ".17g"),
            "log10_p": format(log10_p, ".17g"),
        })
    for row, log10_q in zip(stats, _bh_adjust_log_probabilities(log10_ps)):
        row["p_value_two_sided"] = _probability_text(float(row["log10_p"]))
        row["q_value_BH"] = _probability_text(log10_q)
        row["log10_q"] = format(log10_q, ".17g")
        row["n_tests_BH"] = len(databases)
    return stats


def _write_csv(path, fields, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _legacy_rows(stats):
    return [{
        "category": row["database"], "label": row["label"],
        "cargo_hits": row["cargo_hit_cds"],
        "cargo_total": row["cargo_retained_cds_total"],
        "cargo_frac": row["cargo_hit_fraction"],
        "bg_hits": row["background_hit_cds"],
        "bg_total": row["background_retained_cds_total"],
        "bg_frac": row["background_hit_fraction"],
        "odds_ratio": row["odds_ratio_cargo_vs_background"],
        "p_fisher": row["p_value_two_sided"],
        "q_fdr_bh": row["q_value_BH"],
        "log10_p": row["log10_p"], "log10_q": row["log10_q"],
        "n_tests_BH": row["n_tests_BH"],
    } for row in stats]


def self_test():
    table = [[1, 9], [11, 3]]
    p = math.exp(_fisher_log_p(table))
    assert math.isclose(p, 0.0027594561852200836, rel_tol=1e-10)
    assert math.isclose(float(fisher_exact(table, alternative="two-sided").statistic), 1 / 33,
                        rel_tol=1e-12)
    three = _bh_adjust_log_probabilities([math.log10(0.01), math.log10(0.03), math.log10(0.04)])
    four = _bh_adjust_log_probabilities([math.log10(value) for value in (0.01, 0.02, 0.03, 0.04)])
    assert len(three) == 3 and all(math.isclose(10 ** value, expected, rel_tol=1e-10)
                                    for value, expected in zip(three, (0.03, 0.04, 0.04)))
    assert len(four) == 4 and all(math.isclose(10 ** value, 0.04, rel_tol=1e-10)
                                  for value in four)
    print("self-test passed: two-sided Fisher plus three- and four-test BH families")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path,
                        help="directory containing the two 20_* count tables")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return

    config = json.loads(args.config.read_text(encoding="utf-8"))
    original_categories = [str(value).casefold() for value in config.get("original_categories", [])]
    backbone_categories = [str(value).casefold() for value in config.get("backbone_categories", [])]
    if set(original_categories) != set(DATABASES):
        raise ValueError("config original_categories must contain card, vfdb, bacmet and tadb")
    if not backbone_categories or any(value not in original_categories for value in backbone_categories):
        raise ValueError("config backbone_categories must be a nonempty subset of original_categories")

    before_name = Path(config.get("before_counts", "20_cds_counts_before_filtering.csv")).name
    after_name = Path(config.get("after_counts", "20_cds_counts_after_filtering.csv")).name
    if args.input_dir is not None:
        input_dir = args.input_dir
    else:
        input_dir = next((folder for folder in (ROOT / "results", ROOT / "data")
                          if (folder / before_name).is_file() and (folder / after_name).is_file()),
                         ROOT / "results")
    before = _read_counts(input_dir / before_name)
    after = _read_counts(input_dir / after_name)
    _check_complete(before, DATABASES)
    _check_complete(after, backbone_categories)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    legacy_fields = ["category", "label", "cargo_hits", "cargo_total", "cargo_frac",
                     "bg_hits", "bg_total", "bg_frac", "odds_ratio", "p_fisher",
                     "q_fdr_bh", "log10_p", "log10_q", "n_tests_BH"]
    for mode, filename in (("relaxed", "11_enrichment_relaxed.csv"),
                           ("strict", "13_enrichment_strict.csv")):
        stats = _build_stats(before, DATABASES, mode)
        _write_csv(output_dir / filename, legacy_fields, _legacy_rows(stats))

    filtered_rows = []
    for mode in ("relaxed", "strict"):
        for row in _build_stats(after, backbone_categories, mode):
            filtered_rows.append({"mode": mode, **row})
    filtered_fields = ["mode", "database", "label", "cargo_hit_cds",
                       "cargo_other_retained_cds", "cargo_retained_cds_total",
                       "cargo_hit_fraction", "background_hit_cds",
                       "background_other_retained_cds", "background_retained_cds_total",
                       "background_hit_fraction", "odds_ratio_cargo_vs_background",
                       "p_value_two_sided", "log10_p", "q_value_BH", "log10_q", "n_tests_BH"]
    output_name = Path(config.get("after_statistics", "21_backbone_enrichment.csv")).name
    _write_csv(output_dir / output_name, filtered_fields, filtered_rows)
    print(f"Wrote raw four-category and filtered {len(backbone_categories)}-category Fisher/BH tables to {output_dir}")


if __name__ == "__main__":
    main()
