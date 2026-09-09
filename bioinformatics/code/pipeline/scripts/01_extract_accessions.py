#!/usr/bin/env python3
"""
01_extract_accessions.py
Extract plasmid accessions with >=1 transposase-annotation hit from Rohan's
transposase-counts.csv (the Fig. 1 screen of Ha et al.).

Default paths (relative to the repository root):
    Input:  1_Data/transposase-counts.csv
    Output: 2_Result/1_Hit_plasmids.txt

Usage:
    python scripts/01_extract_accessions.py [input_csv] [output_txt]
"""
import csv
import os
import sys

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
DATA_DIR = os.path.abspath(os.environ.get(
    "NEE_DATA_DIR", os.path.join(ROOT, "1_Data")))
RESULTS_DIR = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))
IN_CSV = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    DATA_DIR, "transposase-counts.csv")
OUT_TXT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    RESULTS_DIR, "1_Hit_plasmids.txt")


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
