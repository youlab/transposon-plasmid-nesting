#!/usr/bin/env python3
"""Purpose: report strict-mode counts at the 90/60 and 80/70 thresholds.
Inputs: 03_plasmid_sequence_hashes.csv, 07_cargo_strict.faa, 13 strict maps and hit tables.
Output: work/processed/14_strict_counts.txt.
Usage: python code/pipeline/scripts/14_strict_counts.py
"""
import csv
import hashlib
import os
from collections import defaultdict

from scipy.stats import fisher_exact

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
CATS = ["card", "vfdb", "bacmet", "tadb"]
LABELS = {"card": "CARD", "vfdb": "VFDB", "bacmet": "BacMet", "tadb": "TADB"}
THRESHOLDS = [(90.0, 60.0), (80.0, 70.0)]


def count_csv_rows(path):
    with open(path) as f:
        return sum(1 for _ in csv.DictReader(f))


def count_fasta_accessions(path):
    accessions = set()
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                accessions.add(line[1:].split("|", 1)[0])
    return len(accessions)


def load_map(path):
    m = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            m[row["md5"]] = int(row["n_copies"])
    return m


def load_hits(path, umap, pid, qcov):
    w, hits = 0, set()
    with open(path, errors="replace") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 7:
                continue
            if float(row[2]) >= pid and float(row[4]) >= qcov:
                w += umap.get(row[0], 1)
                hits.add(row[0])
    return w, hits


def plasmid_counts(hit_sets):
    """unique plasmids (representative accessions) with >=1 cargo hit per cat."""
    accs = defaultdict(set)
    with open(input_path("07_cargo_strict.faa", generated=True)) as f:
        for line in f:
            if line.startswith(">"):
                acc = line[1:].split("|")[0]
            else:
                h = hashlib.md5(line.strip().encode()).hexdigest()
                for cat, s in hit_sets.items():
                    if h in s:
                        accs[cat].add(acc)
    return {cat: len(v) for cat, v in accs.items()}


def main():
    total_plasmids = count_csv_rows(input_path("03_plasmid_sequence_hashes.csv", generated=True))
    strict_cargo_plasmids = count_fasta_accessions(
        input_path("07_cargo_strict.faa", generated=True))
    mc = load_map(os.path.join(RES, "13_cargo_strict_unique_map.csv"))
    mb = load_map(os.path.join(RES, "13_background_strict_unique_map.csv"))
    tc, tb = sum(mc.values()), sum(mb.values())

    lines = [
        f"Strict mode: {strict_cargo_plasmids} plasmids with cargo "
        f"({strict_cargo_plasmids / total_plasmids:.1%} of {total_plasmids}); "
        f"{tc} cargo proteins; {tb} background proteins", ""]
    for pid, qcov in THRESHOLDS:
        lines.append(f"=== pident >= {pid:g}%, qcov >= {qcov:g}% ===")
        hit_sets = {}
        for cat in CATS:
            cw, hs = load_hits(os.path.join(RES, f"13_cargo_strict_vs_{cat}.tsv"),
                               mc, pid, qcov)
            bw, _ = load_hits(os.path.join(RES, f"13_background_strict_vs_{cat}.tsv"),
                              mb, pid, qcov)
            hit_sets[cat] = hs
            odds, p = fisher_exact([[cw, tc - cw], [bw, tb - bw]])
            hit_sets[cat] = hs
            lines.append(
                f"{LABELS[cat]:7s} cargo {cw / tc:.2%} ({cw:,} of {tc:,}) | "
                f"bg {bw / tb:.2%} ({bw:,} of {tb:,}) | "
                f"OR={odds:.2f} | P={p:.1e}")
        plc = plasmid_counts(hit_sets)
        for cat in CATS:
            lines.append(f"  {LABELS[cat]:7s} plasmids with cargo hits: "
                         f"{plc.get(cat, 0):,} "
                         f"({plc.get(cat, 0) / total_plasmids:.1%} of "
                         f"{total_plasmids})")
        lines.append("")
    out = "\n".join(lines)
    print(out)
    with open(os.path.join(RES, "14_strict_counts.txt"), "w") as f:
        f.write(out + "\n")


if __name__ == "__main__":
    main()
