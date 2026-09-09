#!/usr/bin/env python3
"""
07_extract_cargo.py
Define transposon cargo from ISEScan results (TSV) and RefSeq CDS annotations
(GBFF).

Cargo definition (composite-transposon logic, relaxed adjacency):
  - Two SAME-FAMILY IS copies on one plasmid with 0 < span <= MAX_SPAN form a
    candidate pair; adjacency NOT required (nested IS inside the span are
    allowed). Both arcs of the circular plasmid are considered (linear arc
    and wrap-around arc).
  - Cargo = CDS fully inside the UNION of qualifying arcs; CDS overlapping
    any IS copy (transposases) are EXCLUDED from both cargo and background.
  - Each CDS is written once (flags dedup).
  - --strict reproduces the original adjacent-only rule (+ last/first wrap)
    as a sensitivity check (outputs get the 5strict_ prefix).
  - Partial IS copies (tsv type=p) are valid flanks; the cargo table records
    each flank's type.
Scope note: plasmids with ZERO IS copies (the 1,404 no-hit set) contribute
nothing — no cargo by construction, and their CDS are NOT background either
(background = non-cargo, non-IS CDS on IS-carrying plasmids only).

Usage:
    python scripts/07_extract_cargo.py [max_span_bp] [--strict] [--limit N]
    max_span_bp: default 25000; 0 = no span limit

Defaults (relative to the repository root):
    Input:  2_Result/3_Seqhash_map.csv, 1_Data/gbff/*.gbff,
            2_Result/4_ISEScan/fna/*.fna.tsv
    Output: 2_Result/5_Cargo.faa, 5_Background.faa, 5_Cargo_table.csv,
            5_Span_distribution.csv/.pdf, 5_Summary.txt
            (prefix 5strict_ when --strict)
"""
import csv
import glob
import os
import sys

from Bio import SeqIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
DATA_DIR = os.path.abspath(os.environ.get(
    "NEE_DATA_DIR", os.path.join(ROOT, "1_Data")))
RESULTS_DIR = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))
MAP_CSV = os.path.join(RESULTS_DIR, "3_Seqhash_map.csv")
GBFF_DIR = os.path.join(DATA_DIR, "gbff")
TSV_DIR = os.path.join(RESULTS_DIR, "4_ISEScan", "fna")

MAX_SPAN = 25000
STRICT = False
LIMIT = None
args = sys.argv[1:]
pos = [a for a in args if not a.startswith("--")]
if pos:
    MAX_SPAN = int(pos[0])
STRICT = "--strict" in args
if "--limit" in args:
    LIMIT = int(args[args.index("--limit") + 1])

PREFIX = "5strict_" if STRICT else "5_"
if len(pos) > 1:  # optional output-prefix override (for test runs)
    PREFIX = pos[1]
OUT_CARGO = os.path.join(RESULTS_DIR, PREFIX + "Cargo.faa")
OUT_BG = os.path.join(RESULTS_DIR, PREFIX + "Background.faa")
OUT_TABLE = os.path.join(RESULTS_DIR, PREFIX + "Cargo_table.csv")
OUT_SPAN = os.path.join(RESULTS_DIR, PREFIX + "Span_distribution.csv")
OUT_SPANPDF = os.path.join(RESULTS_DIR, PREFIX + "Span_distribution.pdf")
OUT_SUM = os.path.join(RESULTS_DIR, PREFIX + "Summary.txt")

SPAN_DIST_MAX = 200000   # spans beyond this are counted, not listed


def parse_isescan_tsv(path):
    """Return list of dicts (start, end, family, type) sorted by start."""
    elements = []
    with open(path) as f:
        for line in f:
            if line.startswith("seqID") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 22:
                continue
            elements.append({
                "start": int(parts[3]), "end": int(parts[4]),
                "family": parts[1], "type": parts[21].strip().lower(),
            })
    elements.sort(key=lambda e: e["start"])
    return elements


def load_cds(gbff_path):
    """Return (plasmid_length, list of CDS dicts)."""
    rec = SeqIO.read(gbff_path, "genbank")
    cds_list = []
    for feat in rec.features:
        if feat.type != "CDS" or "translation" not in feat.qualifiers:
            continue
        cds_list.append({
            "start": int(feat.location.start) + 1,  # 1-based
            "end": int(feat.location.end),
            "locus": feat.qualifiers.get("locus_tag", ["NA"])[0],
            "product": feat.qualifiers.get("product",
                                           ["hypothetical protein"])[0],
            "translation": feat.qualifiers["translation"][0],
        })
    return len(rec.seq), cds_list


def overlap(a_start, a_end, b_start, b_end):
    return a_start <= b_end and b_start <= a_end


def span_ok(span):
    return span > 0 and (MAX_SPAN == 0 or span <= MAX_SPAN)


def find_arcs(elements, plen):
    """Yield (family, left_type, right_type, span, arc, lo, hi).
    arc='linear': cargo region [lo, hi]; arc='wrap': >= lo or <= hi.
    Relaxed: all same-family pairs, both arcs. Strict: adjacent pairs only,
    plus an explicit circular (last, first) adjacency check."""
    n = len(elements)
    for i in range(n):
        js = (i + 1,) if STRICT else range(i + 1, n)
        for j in js:
            if j >= n:
                break
            a, b = elements[i], elements[j]
            if a["family"] != b["family"]:
                continue
            lin = b["start"] - a["end"] - 1
            if span_ok(lin):
                yield (a["family"], a["type"], b["type"], lin,
                       "linear", a["end"] + 1, b["start"] - 1)
            if not STRICT:
                wrap = (plen - b["end"]) + (a["start"] - 1)
                if span_ok(wrap):
                    yield (a["family"], b["type"], a["type"], wrap,
                           "wrap", b["end"] + 1, a["start"] - 1)
    if STRICT and n >= 2:
        # circular adjacency: last element wraps to first across the origin
        b, a = elements[-1], elements[0]
        if a["family"] == b["family"]:
            wrap = (plen - b["end"]) + (a["start"] - 1)
            if span_ok(wrap):
                yield (a["family"], b["type"], a["type"], wrap,
                       "wrap", b["end"] + 1, a["start"] - 1)


def main():
    md5_to_acc = {}
    with open(MAP_CSV) as f:
        for row in csv.DictReader(f):
            md5_to_acc[row["md5"]] = row["representative_acc"]

    tsv_files = sorted(glob.glob(os.path.join(TSV_DIR, "*.fna.tsv")))
    if LIMIT:
        tsv_files = tsv_files[:LIMIT]

    cargo_records, bg_records, table_rows, span_rows = [], [], [], []
    n_pl = n_no_gbff = n_no_is = n_with_pair = n_with_cargo = 0
    n_span_overcap = 0

    for k, tsv_path in enumerate(tsv_files):
        md5 = os.path.basename(tsv_path)[:-len(".fna.tsv")]
        acc = md5_to_acc.get(md5)
        if not acc:
            continue
        gbff = os.path.join(GBFF_DIR, acc + ".gbff")
        if not os.path.exists(gbff):
            n_no_gbff += 1
            continue
        elements = parse_isescan_tsv(tsv_path)
        if not elements:
            n_no_is += 1
            continue
        n_pl += 1
        plen, cds_list = load_cds(gbff)

        def overlaps_is(c):
            return any(overlap(c["start"], c["end"], e["start"], e["end"])
                       for e in elements)

        # span distribution (all same-family pairs, linear arc)
        n_pair = 0
        for i in range(len(elements)):
            for j in range(i + 1, len(elements)):
                a, b = elements[i], elements[j]
                if a["family"] != b["family"]:
                    continue
                n_pair += 1
                span = b["start"] - a["end"] - 1
                if 0 < span <= SPAN_DIST_MAX:
                    span_rows.append([acc, md5, a["family"], span,
                                      1 if j == i + 1 else 0])
                elif span > SPAN_DIST_MAX:
                    n_span_overcap += 1
        if n_pair:
            n_with_pair += 1

        cargo_flags = [False] * len(cds_list)
        for fam, lt, rt, span, arc, lo, hi in find_arcs(elements, plen):
            for idx, c in enumerate(cds_list):
                if cargo_flags[idx] or overlaps_is(c):
                    continue
                inside = (lo <= c["start"] and c["end"] <= hi) if arc == "linear" \
                    else (c["start"] >= lo or c["end"] <= hi)
                if inside:
                    cargo_flags[idx] = True
                    table_rows.append([acc, md5, c["locus"], c["product"],
                                       fam, lt, rt, span, arc,
                                       c["start"], c["end"]])

        if any(cargo_flags):
            n_with_cargo += 1
        for idx, c in enumerate(cds_list):
            if overlaps_is(c):
                continue
            line = f">{acc}|{c['locus']}|{c['product'][:60]}\n{c['translation']}\n"
            (cargo_records if cargo_flags[idx] else bg_records).append(line)

        if (k + 1) % 2000 == 0:
            print(f"{k + 1}/{len(tsv_files)} plasmids processed", flush=True)

    with open(OUT_CARGO, "w") as f:
        f.writelines(cargo_records)
    with open(OUT_BG, "w") as f:
        f.writelines(bg_records)
    with open(OUT_TABLE, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["accession", "md5", "locus_tag", "product", "is_family",
                    "left_flank_type", "right_flank_type", "span_bp", "arc",
                    "cds_start", "cds_end"])
        w.writerows(table_rows)
    with open(OUT_SPAN, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["accession", "md5", "is_family", "span_bp", "adjacent"])
        w.writerows(span_rows)

    summary = "\n".join([
        "# Step 5: cargo extraction summary",
        f"mode: {'STRICT (adjacent-only)' if STRICT else 'relaxed (non-adjacent allowed)'}",
        f"max_span_bp: {MAX_SPAN if MAX_SPAN else 'no limit'}",
        f"plasmids_with_IS_processed: {n_pl}  (no IS: {n_no_is} skipped; "
        f"gbff missing: {n_no_gbff})",
        f"plasmids_with_samefamily_pair: {n_with_pair}",
        f"plasmids_with_cargo: {n_with_cargo}",
        f"cargo_proteins: {len(cargo_records)}",
        f"background_proteins: {len(bg_records)}",
        f"samefamily_pairs_in_span_distribution: {len(span_rows)} "
        f"(>{SPAN_DIST_MAX} bp not listed: {n_span_overcap})",
        "outputs: " + ", ".join(os.path.basename(path) for path in
                                (OUT_CARGO, OUT_BG, OUT_TABLE)),
    ])
    print(summary)
    with open(OUT_SUM, "w") as f:
        f.write(summary + "\n")

    # span distribution figure (log-x histogram, adjacent vs non-adjacent)
    adj = sorted(r[3] for r in span_rows if r[4] == 1)
    nonadj = sorted(r[3] for r in span_rows if r[4] == 0)
    if adj or nonadj:
        import numpy as np
        bins = np.logspace(0, np.log10(SPAN_DIST_MAX), 60)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist([adj, nonadj], bins=bins, stacked=True,
                label=["adjacent pairs", "non-adjacent pairs"],
                color=["#4C72B0", "#C44E52"])
        if MAX_SPAN:
            ax.axvline(MAX_SPAN, color="black", ls="--", lw=1,
                       label=f"MAX_SPAN={MAX_SPAN}")
        ax.set_xscale("log")
        ax.set_xlabel("same-family pair span (bp)")
        ax.set_ylabel("number of pairs")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT_SPANPDF)
        print(f"figure -> {OUT_SPANPDF}")


if __name__ == "__main__":
    main()
