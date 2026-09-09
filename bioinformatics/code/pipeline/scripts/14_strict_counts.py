#!/usr/bin/env python3
"""
14_strict_counts.py
Mirror of the relaxed results table for STRICT mode: per category, both
thresholds (90/60 and 80/70) — weighted cargo/background hits, Fisher OR,
and the number of unique plasmids carrying >=1 cargo hit (instance-level
scan of 5strict_Cargo.faa). No new BLAST: re-filters existing 8_*_vs_*.tsv.

Usage:  python scripts/14_strict_counts.py
Input:  2_Result/3_Seqhash_map.csv, 2_Result/5strict_Cargo.faa,
        2_Result/8_{Cargo,Background}_unique_map.csv,
        2_Result/8_{Cargo,Background}_vs_{card,vfdb,bacmet,tadb}.tsv
Output: 2_Result/8b_Strict_counts.txt (+ stdout)
"""
import csv
import hashlib
import os
from collections import defaultdict

from scipy.stats import fisher_exact

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
RES = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))
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
    with open(os.path.join(RES, "5strict_Cargo.faa")) as f:
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
    total_plasmids = count_csv_rows(os.path.join(RES, "3_Seqhash_map.csv"))
    strict_cargo_plasmids = count_fasta_accessions(
        os.path.join(RES, "5strict_Cargo.faa"))
    mc = load_map(os.path.join(RES, "8_Cargo_unique_map.csv"))
    mb = load_map(os.path.join(RES, "8_Background_unique_map.csv"))
    tc, tb = sum(mc.values()), sum(mb.values())

    lines = [
        f"Strict mode: {strict_cargo_plasmids} plasmids with cargo "
        f"({strict_cargo_plasmids / total_plasmids:.1%} of {total_plasmids}); "
        f"{tc} cargo proteins; {tb} background proteins", ""]
    for pid, qcov in THRESHOLDS:
        lines.append(f"=== pident >= {pid:g}%, qcov >= {qcov:g}% ===")
        hit_sets = {}
        for cat in CATS:
            cw, hs = load_hits(os.path.join(RES, f"8_Cargo_vs_{cat}.tsv"),
                               mc, pid, qcov)
            bw, _ = load_hits(os.path.join(RES, f"8_Background_vs_{cat}.tsv"),
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
    with open(os.path.join(RES, "8b_Strict_counts.txt"), "w") as f:
        f.write(out + "\n")


if __name__ == "__main__":
    main()
