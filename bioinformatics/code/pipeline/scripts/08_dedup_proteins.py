#!/usr/bin/env python3
"""Purpose: deduplicate relaxed cargo and background proteins by exact sequence.
Inputs: 07_cargo_relaxed.faa and 07_background_relaxed.faa.
Outputs: 08_{cargo,background}_relaxed_unique.faa and corresponding map CSVs.
Usage: python code/pipeline/scripts/08_dedup_proteins.py
"""
import hashlib
import os

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
    mode = "relaxed"
    src = input_path(f"07_{tag.lower()}_{mode}.faa", generated=True)
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
    out_faa = os.path.join(RES, f"08_{tag.lower()}_{mode}_unique.faa")
    seen = set()
    with open(out_faa, "w") as fo:
        for header, seq in read_fasta(src):
            h = hashlib.md5(seq.encode()).hexdigest()
            if h in seen:
                continue
            seen.add(h)
            fo.write(f">{h}\n{seq}\n")
    out_map = os.path.join(RES, f"08_{tag.lower()}_{mode}_unique_map.csv")
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
