#!/usr/bin/env python3
"""
Compare transposase counts from CDS product annotations with IS-copy counts
from ISEScan. This sequence-annotation QC does not infer insertion contexts;
exact insertion-context analysis is performed in Step 15.

Inputs:
    1_Data/transposase-counts.csv
    2_Result/3_Seqhash_map.csv
    2_Result/4_ISEScan/fna/<md5>.fna.tsv

Outputs:
    2_Result/4b_Transposase_count_comparison.csv
    2_Result/4b_Transposase_count_comparison.pdf/.png
    2_Result/4b_Compare_summary.txt

Usage:
    python scripts/05_compare_transposase_counts.py
    python scripts/05_compare_transposase_counts.py \
        transposase-counts.csv 3_Seqhash_map.csv isescan_tsv_dir output_prefix
"""
import csv
import glob
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
DATA_DIR = os.path.abspath(os.environ.get(
    "NEE_DATA_DIR", os.path.join(ROOT, "1_Data")))
RESULTS_DIR = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))

ROHAN_CSV = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    DATA_DIR, "transposase-counts.csv")
MAP_CSV = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    RESULTS_DIR, "3_Seqhash_map.csv")
TSV_DIR = sys.argv[3] if len(sys.argv) > 3 else os.path.join(
    RESULTS_DIR, "4_ISEScan", "fna")
OUT_PREFIX = sys.argv[4] if len(sys.argv) > 4 else os.path.join(
    RESULTS_DIR, "4b_")
def parse_isescan_tsv(path):
    """Return total, complete, and partial IS-copy counts for one plasmid."""
    total = complete = partial = 0
    with open(path) as f:
        for line in f:
            if line.startswith("seqID") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 22:
                continue
            total += 1
            call_type = parts[21].strip().lower()
            if call_type == "c":
                complete += 1
            elif call_type == "p":
                partial += 1
    return total, complete, partial


def main():
    os.makedirs(os.path.dirname(OUT_PREFIX), exist_ok=True)

    is_counts = {}
    for path in glob.glob(os.path.join(TSV_DIR, "*.fna.tsv")):
        md5 = os.path.basename(path)[:-len(".fna.tsv")]
        total, complete, partial = parse_isescan_tsv(path)
        is_counts[md5] = (total, complete, partial)

    acc2md5 = {}
    with open(MAP_CSV) as f:
        for row in csv.DictReader(f):
            for accession in row["all_accessions"].split(";"):
                acc2md5[accession] = row["md5"]

    rows = []
    missing = 0
    with open(ROHAN_CSV) as f:
        records = csv.reader(f)
        next(records)
        for record in records:
            if record[2] != "plasmid" or int(record[6]) == 0:
                continue
            accession = record[1]
            annotation_count = int(record[6])
            md5 = acc2md5.get(accession)
            if md5 is None or md5 not in is_counts:
                missing += 1
                rows.append([accession, annotation_count, "", "", ""])
                continue
            total, complete, partial = is_counts[md5]
            rows.append([accession, annotation_count, total, complete, partial])

    comparison_csv = OUT_PREFIX + "Transposase_count_comparison.csv"
    with open(comparison_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "accession", "annotation_transposase_count", "isescan_is_copies",
            "isescan_complete", "isescan_partial"])
        writer.writerows(rows)

    plotted = [(int(r[1]), int(r[2])) for r in rows if r[2] not in ("", 0)]
    if not plotted:
        raise RuntimeError("No paired ISEScan-positive plasmids were found")
    x = [v[0] for v in plotted]
    y = [v[1] for v in plotted]
    rho, pvalue = spearmanr(x, y)

    summary = "\n".join([
        "# Step 4b: CDS-annotation versus ISEScan count QC",
        f"annotation_positive_accessions: {len(rows)}",
        f"accessions_with_isescan_tsv: {len(rows) - missing}",
        f"missing_isescan_tsv: {missing}",
        f"plotted_isescan_positive_accessions: {len(plotted)}",
        f"spearman_rho_plotted_points: {rho:.4f}",
        f"spearman_pvalue_plotted_points: {pvalue:.3e}",
        f"comparison_csv: {comparison_csv}",
        "note: insertion-context inference is performed only by Step 15",
    ])
    with open(OUT_PREFIX + "Compare_summary.txt", "w") as f:
        f.write(summary + "\n")

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, s=4, alpha=0.15, edgecolors="none", color="#4C72B0")
    limit = [0.8, max(max(x), max(y)) * 1.2]
    ax.plot(limit, limit, color="grey", linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(limit)
    ax.set_ylim(limit)
    ax.set_xlabel("Transposase counts based on CDS annotation")
    ax.set_ylabel("IS element counts based on ISEScan")
    ax.text(0.05, 0.95, f"Spearman $\\rho$ = {rho:.2f}",
            transform=ax.transAxes, ha="left", va="top")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    figure_base = OUT_PREFIX + "Transposase_count_comparison"
    fig.savefig(figure_base + ".pdf")
    fig.savefig(figure_base + ".png", dpi=150)
    plt.close(fig)

    print(summary)


if __name__ == "__main__":
    main()
