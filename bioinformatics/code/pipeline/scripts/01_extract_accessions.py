#!/usr/bin/env python3
"""Purpose: list plasmid accessions with annotated transposases.
Input: data/01_transposase_counts.csv (or NEE_DATA_DIR/NEE_RESULTS_DIR).
Output: work/processed/01_transposase_positive_accessions.txt.
Usage: python code/pipeline/scripts/01_extract_accessions.py [input_csv] [output_txt]
"""
import csv
import os
import sys

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

IN_CSV = sys.argv[1] if len(sys.argv) > 1 else input_path("01_transposase_counts.csv")
OUT_TXT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    RESULTS_DIR, "01_transposase_positive_accessions.txt")


def main():
    with open(IN_CSV) as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    # columns: AnnotationAccession, SeqID, SeqType, SeqLength, CDS_count,
    #          CDS_length, transposase_count, " transposase_length" (note space)
    seen = set()
    accessions = []
    for r in data:
        if r[2] == "plasmid" and int(r[6]) > 0:
            if r[1] not in seen:
                seen.add(r[1])
                accessions.append(r[1])

    os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(accessions) + "\n")
    print(f"{len(accessions)} accessions written -> {OUT_TXT}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
