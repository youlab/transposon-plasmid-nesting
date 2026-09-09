#!/usr/bin/env python3
"""
08_dedup_proteins.py
Deduplicate identical protein sequences in 5_Cargo.faa / 5_Background.faa
(exact md5 on the amino-acid sequence) before DIAMOND. Hit status of each
unique sequence maps back to all its copies via the map csv, so downstream
counts are unchanged — this is purely a compute saver.

Usage:  python scripts/08_dedup_proteins.py

Input:  2_Result/5_Cargo.faa, 2_Result/5_Background.faa
Output: 2_Result/6_Cargo_unique.faa, 6_Cargo_unique_map.csv
        2_Result/6_Background_unique.faa, 6_Background_unique_map.csv
        (unique fasta id = md5; map csv: md5, n_copies, representative_header)
"""
import hashlib
import os

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
RES = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))


def read_fasta(path):
    """Yield (header_without_>, sequence)."""
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def dedup(tag):
    src = os.path.join(RES, f"5_{tag}.faa")
    # first pass: count copies per md5, keep first header as representative
    by_hash = {}  # md5 -> [count, first_header]
    n_total = 0
    for header, seq in read_fasta(src):
        n_total += 1
        h = hashlib.md5(seq.encode()).hexdigest()
        if h in by_hash:
            by_hash[h][0] += 1
        else:
            by_hash[h] = [1, header]
    # second pass: write unique fasta in input order (id = md5)
    out_faa = os.path.join(RES, f"6_{tag}_unique.faa")
    seen = set()
    with open(out_faa, "w") as fo:
        for header, seq in read_fasta(src):
            h = hashlib.md5(seq.encode()).hexdigest()
            if h in seen:
                continue
            seen.add(h)
            fo.write(f">{h}\n{seq}\n")
    out_map = os.path.join(RES, f"6_{tag}_unique_map.csv")
    with open(out_map, "w") as fm:
        fm.write("md5,n_copies,representative_header\n")
        for h, (n, hdr) in sorted(by_hash.items()):
            fm.write(f"{h},{n},{hdr}\n")
    print(f"{tag}: {n_total} proteins -> {len(by_hash)} unique "
          f"({n_total - len(by_hash)} duplicates removed, "
          f"{(n_total - len(by_hash)) / max(n_total, 1):.1%})")


def main():
    dedup("Cargo")
    dedup("Background")


if __name__ == "__main__":
    main()
