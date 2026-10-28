#!/usr/bin/env python3
"""Purpose: calculate copy-weighted AMRFinder ARG enrichment.
Inputs: 09_amrfinder_{cargo,background}.tsv, 08_*_relaxed_unique_map.csv, and 07_cargo_relaxed.faa.
Outputs: 12_amrfinder_enrichment.csv and 12_amrfinder_enrichment_summary.txt.
Usage: python code/pipeline/scripts/12_amrf_enrichment.py
"""
import csv
import hashlib
import os

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


def load_copy_map(path):
    with open(path) as f:
        return {row["md5"]: int(row["n_copies"])
                for row in csv.DictReader(f)}


def load_amrf_queries(path):
    """Return unique query protein identifiers reported by AMRFinderPlus."""
    queries = set()
    with open(path, errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            query = line.rstrip("\n").split("\t", 1)[0].strip()
            if query.lower() in {"protein identifier", "protein_identifier"}:
                continue
            queries.add(query)
    return queries


def cargo_plasmid_count(hit_queries):
    accessions = set()
    accession = None
    sequence = []
    with open(input_path("07_cargo_relaxed.faa", generated=True)) as f:
        for line in f:
            if line.startswith(">"):
                if accession is not None:
                    digest = hashlib.md5("".join(sequence).encode()).hexdigest()
                    if digest in hit_queries:
                        accessions.add(accession)
                accession = line[1:].split("|", 1)[0]
                sequence = []
            else:
                sequence.append(line.strip())
        if accession is not None:
            digest = hashlib.md5("".join(sequence).encode()).hexdigest()
            if digest in hit_queries:
                accessions.add(accession)
    return len(accessions)


def main():
    maps = {
        name: load_copy_map(input_path(f"08_{name.lower()}_relaxed_unique_map.csv", generated=True))
        for name in ("Cargo", "Background")
    }
    hits = {
        name: load_amrf_queries(input_path(f"09_amrfinder_{name.lower()}.tsv", generated=True))
        for name in ("Cargo", "Background")
    }
    totals = {name: sum(maps[name].values()) for name in maps}
    weighted = {
        name: sum(maps[name].get(query, 0) for query in hits[name])
        for name in maps
    }
    odds_ratio, pvalue = fisher_exact([
        [weighted["Cargo"], totals["Cargo"] - weighted["Cargo"]],
        [weighted["Background"],
         totals["Background"] - weighted["Background"]],
    ])
    n_plasmids = cargo_plasmid_count(hits["Cargo"])

    out_csv = os.path.join(RES, "12_amrfinder_enrichment.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "cargo_hits", "cargo_total", "cargo_fraction",
            "background_hits", "background_total", "background_fraction",
            "odds_ratio", "p_fisher", "cargo_plasmids_with_hit"])
        writer.writerow([
            weighted["Cargo"], totals["Cargo"],
            weighted["Cargo"] / totals["Cargo"],
            weighted["Background"], totals["Background"],
            weighted["Background"] / totals["Background"],
            odds_ratio, pvalue, n_plasmids])

    summary = "\n".join([
        "# AMRFinderPlus copy-weighted ARG enrichment",
        f"cargo: {weighted['Cargo']:,}/{totals['Cargo']:,} "
        f"({weighted['Cargo'] / totals['Cargo']:.2%})",
        f"background: {weighted['Background']:,}/{totals['Background']:,} "
        f"({weighted['Background'] / totals['Background']:.2%})",
        f"odds_ratio: {odds_ratio:.2f}",
        f"fisher_pvalue: {pvalue:.3e}",
        f"cargo_plasmids_with_hit: {n_plasmids:,}",
    ])
    with open(os.path.join(RES, "12_amrfinder_enrichment_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
