#!/usr/bin/env python3
"""
10_function_inventory.py
Functional inventory of transposon cargo. Enrichment relative to the
within-plasmid background is calculated in Step 11.

Per category (CARD=antibiotic resistance, VFDB=virulence, BacMet=metal/
biocide resistance, TADB=toxin-antitoxin):
  - which cargo proteins hit (DIAMOND, pident>=P and qcov>=Q),
  - weighted back to full copy counts via 6_Cargo_unique_map.csv,
  - joined to 5_Cargo_table.csv for plasmid / IS-family context,
  - top genes by copy count, plasmids carrying >=1 hit,
  - cargo annotated fraction (any DB) and hypothetical fraction.

Usage:  python scripts/10_function_inventory.py [min_pident] [min_qcov] [tag]
        default 90 60, no tag (main outputs); tag renames outputs to 6a<tag>_*

Input:  2_Result/6_Cargo_vs_{card,vfdb,bacmet,tadb}.tsv
        2_Result/6_Cargo_unique_map.csv
        2_Result/5_Cargo_table.csv
Output: 2_Result/6a_Cargo_function_inventory.csv
        2_Result/6a_Cargo_top_genes.csv
        2_Result/6a_Summary.txt
"""
import csv
import hashlib
import os
import re
import sys
from collections import Counter, defaultdict

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
RES = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))
MIN_PIDENT = float(sys.argv[1]) if len(sys.argv) > 1 else 90.0
MIN_QCOV = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
TAG = sys.argv[3] if len(sys.argv) > 3 else ""
OUT_PFX = f"6a{TAG}_"

CATS = {"card": "Antibiotic resistance (CARD)",
        "vfdb": "Virulence (VFDB)",
        "bacmet": "Metal/biocide resistance (BacMet)",
        "tadb": "Toxin-antitoxin (TADB)"}


def gene_of(cat, stitle):
    """Extract a short gene name from a subject title."""
    try:
        if cat == "card":      # gb|ACC|ARO:xxx|Gene [org]
            return stitle.split("|")[3].split(" [")[0]
        if cat == "vfdb":      # VFGx(gb|ACC) (gene) desc [...]
            m = re.match(r"\S+\([^()]*\)\s*\(([^()]+)\)", stitle)
            if m:
                return m.group(1)
            m = re.search(r"\(([^()]+)\)", stitle)
            return m.group(1) if m else stitle[:30]
        if cat == "bacmet":    # BAC0001|abeM|tr|...
            return stitle.split("|")[1]
        if cat == "tadb":      # TADB|T6317 ...| Gene (plasmid) [org]
            s = stitle.rsplit("|", 1)[1].strip()
            g = re.split(r"[\[(]", s)[0].strip()
            if not g or "locus_tag" in g or g.lower().startswith("locus"):
                m = re.match(r"TADB\|(\S+)", stitle)
                return m.group(1) if m else g[:40]
            return g[:60]
    except Exception:
        pass
    return stitle[:40]


def load_hits(cat):
    """md5 -> (pident, qcov, gene, stitle) for hits passing thresholds."""
    hits = {}
    with open(os.path.join(RES, f"6_Cargo_vs_{cat}.tsv"),
              errors="replace") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 7:
                continue
            pident, qcov = float(row[2]), float(row[4])
            if pident >= MIN_PIDENT and qcov >= MIN_QCOV:
                hits[row[0]] = (pident, qcov,
                                gene_of(cat, row[7] if len(row) > 7 else ""),
                                row[7] if len(row) > 7 else "")
    return hits


def main():
    # md5 -> [n_copies, representative_header]
    umap = {}
    with open(os.path.join(RES, "6_Cargo_unique_map.csv")) as f:
        for row in csv.DictReader(f):
            umap[row["md5"]] = [int(row["n_copies"]),
                                row["representative_header"]]

    # (acc, locus) -> cargo-table row (is_family, flanks, span, product)
    cargo_ctx = {}
    n_cargo_cds = n_hypo = 0
    with open(os.path.join(RES, "5_Cargo_table.csv")) as f:
        for r in csv.reader(f):
            if r[0] == "accession":
                continue
            n_cargo_cds += 1
            if "hypothetical" in r[3].lower():
                n_hypo += 1
            cargo_ctx[(r[0], r[2])] = r

    hits = {cat: load_hits(cat) for cat in CATS}

    # instance-level scan of cargo faa for exact plasmid sets per category
    # (recompute md5 per protein instance; every instance = one cargo CDS)
    plasmids_of = defaultdict(set)   # cat -> {acc}
    weighted = Counter()             # cat -> copies
    annotated_md5 = set()
    with open(os.path.join(RES, "5_Cargo.faa")) as f:
        for line in f:
            if line.startswith(">"):
                acc = line[1:].split("|")[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
                # protein records are single-line sequences in 5_Cargo.faa
                h = hashlib.md5("".join(seq_lines).encode()).hexdigest()
                seq_lines = []
                for cat, hmap in hits.items():
                    if h in hmap:
                        plasmids_of[cat].add(acc)
                        weighted[cat] += 1
                        annotated_md5.add(h)

    # inventory rows (one per md5 x category)
    inv_rows = []
    for cat, hmap in hits.items():
        for md5, (pident, qcov, gene, stitle) in sorted(hmap.items()):
            n_copies, rep = umap.get(md5, [0, "||"])
            parts = rep.split("|")
            racc = parts[0] if len(parts) > 0 else ""
            rlocus = parts[1] if len(parts) > 1 else ""
            rprod = parts[2] if len(parts) > 2 else ""
            ctx = cargo_ctx.get((racc, rlocus),
                                ["", "", "", "", "", "", "", "", "", "", ""])
            inv_rows.append([md5, cat, gene, pident, qcov, n_copies,
                             racc, rlocus, rprod,
                             ctx[4], ctx[5], ctx[6], stitle[:80]])

    with open(os.path.join(RES, OUT_PFX + "Cargo_function_inventory.csv"),
              "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["md5", "category", "gene", "pident", "qcov", "n_copies",
                    "rep_accession", "rep_locus", "rep_product",
                    "is_family", "left_flank", "right_flank", "subject_title"])
        w.writerows(inv_rows)

    # top genes per category (weighted by copies)
    top_rows = []
    gene_copies = {cat: Counter() for cat in CATS}
    for cat, hmap in hits.items():
        for md5, (_, _, gene, _) in hmap.items():
            gene_copies[cat][gene] += umap.get(md5, [1, ""])[0]
    for cat in CATS:
        for gene, n in gene_copies[cat].most_common(25):
            top_rows.append([cat, gene, n,
                             len(plasmids_of[cat]), CATS[cat]])
    with open(os.path.join(RES, OUT_PFX + "Cargo_top_genes.csv"), "w",
              newline="") as f:
        w = csv.writer(f)
        w.writerow(["category", "gene", "weighted_copies",
                    "plasmids_with_category", "category_label"])
        w.writerows(top_rows)

    n_cargo_proteins = sum(v[0] for v in umap.values())
    lines = [
        "# Step 6a: cargo functional inventory "
        f"(pident>={MIN_PIDENT:g}, qcov>={MIN_QCOV:g})",
        f"cargo_proteins_total: {n_cargo_proteins} "
        f"(unique: {len(umap)}; cargo CDS: {n_cargo_cds})",
        f"cargo_annotated_any_db: "
        f"{sum(umap[m][0] for m in annotated_md5)} copies "
        f"({sum(umap[m][0] for m in annotated_md5) / max(n_cargo_proteins, 1):.1%}), "
        f"{len(annotated_md5)} unique",
        f"cargo_hypothetical_product: {n_hypo}/{n_cargo_cds} "
        f"({n_hypo / max(n_cargo_cds, 1):.1%} of cargo CDS)",
        "",
    ]
    for cat, label in CATS.items():
        hmap = hits[cat]
        lines.append(f"## {label}")
        lines.append(f"  unique proteins with hit: {len(hmap)}")
        lines.append(f"  weighted copies: {weighted[cat]} "
                     f"({weighted[cat] / max(n_cargo_proteins, 1):.2%} of cargo)")
        lines.append(f"  plasmids carrying >=1: {len(plasmids_of[cat])}")
        lines.append(f"  top genes: " + "; ".join(
            f"{g}({n})" for g, n in gene_copies[cat].most_common(10)))
        lines.append("")
    summary = "\n".join(lines)
    with open(os.path.join(RES, OUT_PFX + "Summary.txt"), "w") as f:
        f.write(summary + "\n")
    print(summary)


if __name__ == "__main__":
    main()
