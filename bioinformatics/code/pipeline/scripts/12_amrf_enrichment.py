#!/usr/bin/env python3
"""
Calculate copy-weighted ARG enrichment from AMRFinderPlus output.

Inputs:
    2_Result/6_AMRF_Cargo.tsv
    2_Result/6_AMRF_Background.tsv
    2_Result/6_Cargo_unique_map.csv
    2_Result/6_Background_unique_map.csv
    2_Result/5_Cargo.faa

Outputs:
    2_Result/6_AMRF_enrichment.csv
    2_Result/6_AMRF_weighted_summary.txt

Usage:
    python scripts/12_amrf_enrichment.py
"""
import csv
import hashlib
import os

from scipy.stats import fisher_exact


ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
RES = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))


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
    with open(os.path.join(RES, "5_Cargo.faa")) as f:
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
        name: load_copy_map(os.path.join(RES, f"6_{name}_unique_map.csv"))
        for name in ("Cargo", "Background")
    }
    hits = {
        name: load_amrf_queries(os.path.join(RES, f"6_AMRF_{name}.tsv"))
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

    out_csv = os.path.join(RES, "6_AMRF_enrichment.csv")
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
    with open(os.path.join(RES, "6_AMRF_weighted_summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
