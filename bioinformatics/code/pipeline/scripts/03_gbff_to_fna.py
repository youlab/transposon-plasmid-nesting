#!/usr/bin/env python3
"""
03_gbff_to_fna.py
Convert downloaded .gbff files to .fna for ISEScan, and deduplicate by exact
sequence md5 (identical plasmid sequences sampled from many genomes are kept
once; a mapping table records all accessions sharing a sequence).

Defaults (relative to the repository root):
    Input:  1_Data/gbff/*.gbff
    Output: 1_Data/fna/<md5>.fna         (header: >representative_acc|<md5>)
            2_Result/3_Seqhash_map.csv  (md5, representative_acc, n_accessions, all_accessions)
            2_Result/3_Dedup_summary.txt

Usage:  python scripts/03_gbff_to_fna.py

Notes:
    - fna files are written on first sight of each unique sequence (sequences
      are not held in memory).
    - gbff files that fail to parse, or have empty/N-only sequence, are
      skipped and listed in the summary.
"""
import datetime
import glob
import hashlib
import os

from Bio import SeqIO

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
DATA_DIR = os.path.abspath(os.environ.get(
    "NEE_DATA_DIR", os.path.join(ROOT, "1_Data")))
RESULTS_DIR = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))
GBFF_DIR = os.path.join(DATA_DIR, "gbff")
FNA_DIR = os.path.join(DATA_DIR, "fna")
MAP_CSV = os.path.join(RESULTS_DIR, "3_Seqhash_map.csv")
SUMMARY_TXT = os.path.join(RESULTS_DIR, "3_Dedup_summary.txt")


def main():
    os.makedirs(FNA_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(MAP_CSV), exist_ok=True)

    by_hash = {}   # md5 -> [accessions]  (first = representative)
    skipped = []   # (accession, reason)
    files = sorted(glob.glob(os.path.join(GBFF_DIR, "*.gbff")))

    for i, path in enumerate(files):
        acc = os.path.basename(path)[:-len(".gbff")]
        try:
            rec = SeqIO.read(path, "genbank")
            seq = str(rec.seq).upper()
        except Exception as e:
            skipped.append((acc, f"parse error: {e}"))
            continue
        if not seq or set(seq) <= {"N"}:
            skipped.append((acc, "empty or N-only sequence"))
            continue

        h = hashlib.md5(seq.encode()).hexdigest()
        if h not in by_hash:
            by_hash[h] = []
            with open(os.path.join(FNA_DIR, h + ".fna"), "w") as fo:
                fo.write(f">{acc}|{h}\n")
                for j in range(0, len(seq), 70):
                    fo.write(seq[j:j + 70] + "\n")
        by_hash[h].append(acc)

        if (i + 1) % 2000 == 0:
            print(f"{i + 1}/{len(files)} parsed", flush=True)

    with open(MAP_CSV, "w") as f:
        f.write("md5,representative_acc,n_accessions,all_accessions\n")
        for h, accs in sorted(by_hash.items()):
            f.write(f"{h},{accs[0]},{len(accs)},{';'.join(accs)}\n")

    n_parsed = sum(len(v) for v in by_hash.values())
    n_unique = len(by_hash)
    lines = [
        "# Step 3: gbff -> fna + exact-md5 dedup summary",
        f"date: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"gbff_files_input: {len(files)}",
        f"accessions_before_dedup: {n_parsed}",
        f"accessions_skipped: {len(skipped)}",
        f"unique_sequences_after_dedup: {n_unique}",
        f"exact_duplicates_removed: {n_parsed - n_unique}",
    ]
    if skipped:
        lines.append("skipped_accessions:")
        lines.extend(f"  {acc}\t{reason}" for acc, reason in skipped)
    with open(SUMMARY_TXT, "w") as f:
        f.write("\n".join(lines) + "\n")

    print("\n".join(lines[:7]))
    print(f"summary -> {SUMMARY_TXT}")


if __name__ == "__main__":
    main()
