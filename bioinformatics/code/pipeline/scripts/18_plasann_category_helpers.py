#!/usr/bin/env python3
"""Summarize raw PlasAnn categories against the frozen cargo/background groups."""
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PROCESSED_DIR = Path(os.environ.get("NEE_RESULTS_DIR", ROOT / "work" / "processed"))
RAW_DIR = Path(os.environ.get("NEE_DATA_DIR", ROOT / "work" / "raw"))
MAP_CSV = (PROCESSED_DIR / "03_plasmid_sequence_hashes.csv"
           if (PROCESSED_DIR / "03_plasmid_sequence_hashes.csv").is_file()
           else ROOT / "data" / "03_plasmid_sequence_hashes.csv")
CARGO_FASTA = PROCESSED_DIR / "07_cargo_relaxed.faa"
BACKGROUND_FASTA = PROCESSED_DIR / "07_background_relaxed.faa"
CARGO_TABLE = PROCESSED_DIR / "07_cargo_membership_relaxed.csv.gz"
GBFF_DIR = RAW_DIR / "gbff"
DEFAULT_INPUT = ROOT / "results" / "19_plasann_coverage_annotation"
DEFAULT_OUTPUT = ROOT / "results" / "18_plasann_category_summary"
GROUPS = ("all_CDS", "translated_CDS", "cargo", "background", "outside_original_analysis")
FLAGS = ("conjugation", "plasmid_maintenance_replication_regulation", "hypothetical",
         "unclassified", "exclusion_union", "original_hypothetical")
PROVENANCE_KEYS = ("reference_sha256", "plasann_version", "plasann_commit", "config_sha256", "config")


def norm(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).casefold())


def product_for_header(product):
    return product or "hypothetical protein"


def header_key(accession, locus, product):
    return str(accession), (str(locus) or "NA"), product_for_header(str(product))[:60]


def is_unclassified(category):
    return norm(category) in {"", "na", "unknown", "unclassified", "uncharacterized", "openreadingframe"}


def is_maintenance_category(category):
    return norm(category) in {
        norm("Plasmid backbone and regulation"),
        norm("Plasmid Maintenance, Replication and Regulation"),
        norm("4 Plasmid backbone and regulation"),
    }


def classify(row):
    category = norm(row["category"])
    new_text = f"{row['gene_name']} {row['product']}".casefold()
    original_text = row["original_product"].casefold()
    conjugation = category == norm("Conjugation")
    maintenance = is_maintenance_category(row["category"])
    hypothetical = "hypothetical" in new_text
    unclassified = is_unclassified(row["category"])
    return {
        "conjugation": conjugation,
        "plasmid_maintenance_replication_regulation": maintenance,
        "hypothetical": hypothetical,
        "unclassified": unclassified,
        "exclusion_union": conjugation or maintenance or hypothetical or unclassified,
        "original_hypothetical": "hypothetical" in original_text,
    }


def groups_for(row, coord, cargo_coords, background_coords):
    groups = {"all_CDS"}
    if row["eligible"]:
        groups.add("translated_CDS")
    if coord in cargo_coords:
        groups.add("cargo")
    elif coord in background_coords:
        groups.add("background")
    else:
        groups.add("outside_original_analysis")
    return groups


def resolve_background(header_count, candidates, accession, key):
    if header_count > len(candidates):
        raise ValueError(f"unmatched background FASTA header: {accession} {key}")
    if header_count and header_count != len(candidates):
        raise ValueError(f"ambiguous repeated background header key: {accession} {key}")
    return candidates if header_count else []


def read_expected():
    expected = {}
    with MAP_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            accession, md5 = row["representative_acc"], row["md5"]
            if accession in expected:
                raise ValueError(f"duplicate representative accession in sequence map: {accession}")
            expected[accession] = md5
    return expected


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_runs(input_dir, expected, allow_partial):
    plasmids = input_dir / "plasmids"
    if not plasmids.is_dir():
        raise ValueError(f"missing per-plasmid output directory: {plasmids}")
    runs, statuses = {}, {}
    for accession, md5 in expected.items():
        folder = plasmids / accession
        manifest_path, table_path = folder / "success.json", folder / "cds_annotations.csv"
        try:
            if (folder / "failure.json").exists():
                raise ValueError("failure.json is present alongside output")
            if not manifest_path.is_file() or not table_path.is_file():
                raise ValueError("missing success.json or cds_annotations.csv")
            manifest = json.loads(manifest_path.read_text())
            if str(manifest.get("status", "")).casefold() != "success":
                raise ValueError("manifest status is not success")
            if manifest.get("accession") != accession or manifest.get("md5") != md5:
                raise ValueError("manifest accession/md5 does not match sequence map")
            for key in (*PROVENANCE_KEYS, "input_sha256", "cds_annotations_sha256", "counts"):
                if key not in manifest or manifest[key] in (None, ""):
                    raise ValueError(f"manifest lacks {key}")
            source_path = GBFF_DIR / f"{accession}.gbff"
            if not source_path.is_file() or sha256_file(source_path) != manifest["input_sha256"]:
                raise ValueError("source GenBank checksum does not match success.json")
            if sha256_file(table_path) != manifest["cds_annotations_sha256"]:
                raise ValueError("cds_annotations.csv checksum does not match success.json")
            counts = manifest["counts"]
            if not isinstance(counts, dict) or "cds_total" not in counts or "protein_eligible" not in counts:
                raise ValueError("manifest counts lack cds_total/protein_eligible")
            signature = json.dumps({k: manifest[k] for k in PROVENANCE_KEYS}, sort_keys=True, separators=(",", ":"))
            runs[accession] = (md5, table_path, manifest, signature)
            statuses[accession] = "success"
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            statuses[accession] = str(exc)
    for folder in plasmids.iterdir():
        if folder.is_dir() and folder.name not in expected:
            statuses[f"EXTRA:{folder.name}"] = "unexpected plasmid directory"
    issues = [f"{acc}: {status}" for acc, status in statuses.items() if status != "success"]
    if issues and not allow_partial:
        raise ValueError("step11 output is incomplete; rerun or pass --allow-partial. " + "; ".join(issues[:8]))
    if not runs:
        raise ValueError("no successful step11 plasmid outputs found; " + "; ".join(issues[:8]))
    if len({run[3] for run in runs.values()}) != 1:
        raise ValueError("step11 manifests mix database/version/script configuration provenance")
    return runs, statuses


def read_headers(path, selected):
    counts = defaultdict(Counter)
    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            parts = line[1:].rstrip("\r\n").split("|", 2)
            if parts[0] not in selected:
                continue
            if len(parts) != 3:
                raise ValueError(f"malformed FASTA header in {path}: {line[:120].rstrip()}")
            accession, locus, product = parts
            counts[accession][(locus or "NA", product[:60])] += 1
    return counts


def read_cargo_table(selected):
    coords, header_counts = defaultdict(set), defaultdict(Counter)
    opener = gzip.open if CARGO_TABLE.suffix == ".gz" else open
    with opener(CARGO_TABLE, "rt", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            accession = row["accession"]
            if accession not in selected:
                continue
            locus = row["locus_tag"] or "NA"
            start, end = int(row["cds_start"]), int(row["cds_end"])
            coord = (locus, start, end)
            if coord in coords[accession]:
                raise ValueError(f"duplicate cargo coordinates: {accession} {coord}")
            coords[accession].add(coord)
            header_counts[accession][(locus, row["product"][:60])] += 1
    return coords, header_counts


def columns(fieldnames):
    actual = {norm(name): name for name in (fieldnames or [])}
    required = ("accession", "md5", "cds_index", "locus_tag", "original_product", "start", "end",
                "protein_eligible", "gene_name", "product", "category")
    out = {}
    for name in required:
        if norm(name) not in actual:
            raise ValueError(f"cds_annotations.csv lacks required column {name}")
        out[name] = actual[norm(name)]
    return out


def read_annotations(accession, md5, path, manifest):
    rows, coords, by_header, cds_indices = [], {}, defaultdict(list), set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        cols = columns(reader.fieldnames)
        for raw in reader:
            row = {name: raw[col] or "" for name, col in cols.items()}
            if row["accession"] != accession or row["md5"] != md5:
                raise ValueError(f"annotation identity mismatch: {path}")
            row["start"], row["end"] = int(row["start"]), int(row["end"])
            cds_index = int(row["cds_index"])
            if cds_index in cds_indices:
                raise ValueError(f"duplicate CDS index: {accession} {cds_index}")
            cds_indices.add(cds_index)
            row["eligible"] = str(row["protein_eligible"]).casefold() in {"true", "1", "yes", "y"}
            row["locus_tag"] = row["locus_tag"] or "NA"
            row["key"] = (row["locus_tag"], product_for_header(row["original_product"])[:60])
            coord = (row["locus_tag"], row["start"], row["end"])
            if coord in coords:
                raise ValueError(f"duplicate annotation coordinates: {accession} {coord}")
            coords[coord] = row
            by_header[row["key"]].append(row)
            rows.append(row)
    counts = manifest["counts"]
    if len(rows) != int(counts["cds_total"]):
        raise ValueError(f"annotation CSV row count disagrees with manifest: {accession}")
    if sum(row["eligible"] for row in rows) != int(counts["protein_eligible"]):
        raise ValueError(f"translated CDS row count disagrees with manifest: {accession}")
    return rows, coords, by_header


def update_counts(category_counts, flag_counts, group_counts, accession, row, groups):
    flags = classify(row)
    for group in groups:
        total, translated, plasmids = group_counts[group]
        group_counts[group] = (total + 1, translated + int(row["eligible"]), plasmids | {accession})
        category = row["category"]
        cds, eligible, plasmids = category_counts[(group, category)]
        category_counts[(group, category)] = (cds + 1, eligible + int(row["eligible"]), plasmids | {accession})
        for flag, value in flags.items():
            if value:
                cds, eligible, plasmids = flag_counts[(group, flag)]
                flag_counts[(group, flag)] = (cds + 1, eligible + int(row["eligible"]), plasmids | {accession})


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields,
                                delimiter="\t" if path.suffix == ".tsv" else ",")
        writer.writeheader()
        writer.writerows(rows)


def self_test():
    assert is_unclassified("Uncharacterized") and is_unclassified("")
    assert is_unclassified("open reading frame")
    assert is_maintenance_category("4 PLasmid backbone and regulation")
    assert is_maintenance_category("Plasmid Maintenance, Replication and Regulation")
    assert not is_maintenance_category("other plasmid backbone-adjacent description")
    assert classify({"category": "Conjugation", "gene_name": "", "product": "", "original_product": ""})["conjugation"]
    ta = classify({"category": "Toxin-Antitoxin System", "gene_name": "toxin", "product": "antitoxin", "original_product": "known protein"})
    assert not ta["exclusion_union"]
    ta_hypothetical = classify({"category": "Toxin-Antitoxin System", "gene_name": "", "product": "hypothetical antitoxin", "original_product": "known protein"})
    assert ta_hypothetical["hypothetical"] and ta_hypothetical["exclusion_union"]
    overlap = classify({"category": "Uncharacterized", "gene_name": "", "product": "hypothetical protein", "original_product": "known enzyme"})
    assert overlap["hypothetical"] and overlap["unclassified"] and overlap["exclusion_union"]
    assert not overlap["original_hypothetical"]
    assert header_key("acc", "", "") == ("acc", "NA", "hypothetical protein")
    assert groups_for({"eligible": True}, ("L1", 1, 2), {("L1", 1, 2)}, set()) == {"all_CDS", "translated_CDS", "cargo"}
    assert groups_for({"eligible": False}, ("L2", 3, 4), set(), set()) == {"all_CDS", "outside_original_analysis"}
    assert len(resolve_background(2, [1, 2], "acc", "key")) == 2
    try:
        resolve_background(1, [1, 2], "acc", "key")
    except ValueError:
        pass
    else:
        raise AssertionError("ambiguous background key was accepted")
    print("self-test passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--allow-partial", action="store_true",
                        help="summarize successful accessions only; denominators are observed-subset counts")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return

    expected = read_expected()
    runs, statuses = read_runs(args.input_dir, expected, args.allow_partial)
    selected = set(runs)
    cargo_headers = read_headers(CARGO_FASTA, selected)
    background_headers = read_headers(BACKGROUND_FASTA, selected)
    cargo_coords, cargo_table_headers = read_cargo_table(selected)
    category_counts = defaultdict(lambda: (0, 0, set()))
    flag_counts = defaultdict(lambda: (0, 0, set()))
    group_counts = defaultdict(lambda: (0, 0, set()))

    for accession, (md5, path, manifest, _signature) in runs.items():
        rows, coords, by_header = read_annotations(accession, md5, path, manifest)
        if cargo_headers.get(accession, Counter()) != cargo_table_headers.get(accession, Counter()):
            raise ValueError(f"cargo FASTA headers do not match Cargo_table.csv: {accession}")
        matched_cargo_headers = Counter()
        for coord in cargo_coords.get(accession, set()):
            row = coords.get(coord)
            if row is None or not row["eligible"]:
                raise ValueError(f"cargo table coordinate has no translated CDS: {accession} {coord}")
            matched_cargo_headers[row["key"]] += 1
        if matched_cargo_headers != cargo_table_headers.get(accession, Counter()):
            raise ValueError(f"cargo table products do not match annotations: {accession}")

        background_coords = set()
        for key, count in background_headers.get(accession, Counter()).items():
            candidates = [row for row in by_header.get(key, []) if row["eligible"] and
                          (row["locus_tag"], row["start"], row["end"]) not in cargo_coords.get(accession, set())]
            resolved = resolve_background(count, candidates, accession, key)
            background_coords.update((row["locus_tag"], row["start"], row["end"]) for row in resolved)

        for row in rows:
            coord = (row["locus_tag"], row["start"], row["end"])
            groups = groups_for(row, coord, cargo_coords.get(accession, set()), background_coords)
            update_counts(category_counts, flag_counts, group_counts, accession, row, groups)

    all_cds = group_counts["all_CDS"][0]
    if all_cds != sum(group_counts[group][0] for group in
                      ("cargo", "background", "outside_original_analysis")):
        raise ValueError("cargo, background, and outside-original groups do not partition all_CDS")
    if group_counts["translated_CDS"][0] < group_counts["cargo"][0] + group_counts["background"][0]:
        raise ValueError("cargo/background membership exceeds translated_CDS")

    partial = len(runs) != len(expected) or any(value != "success" for value in statuses.values())
    status = "PARTIAL" if partial else "COMPLETE"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cat_rows = []
    categories = {category for _, category in category_counts}
    for group in GROUPS:
        for category in sorted(categories, key=lambda value: (value.casefold(), value)):
            cds, eligible, plasmids = category_counts.get((group, category), (0, 0, set()))
            total, total_eligible, _ = group_counts[group]
            cat_rows.append({"run_status": status, "group": group, "raw_category": category,
                             "cds_count": cds, "translated_cds_count": eligible,
                             "group_cds_total": total, "group_translated_cds_total": total_eligible,
                             "cds_fraction": f"{cds / total:.8f}" if total else "",
                             "translated_cds_fraction": f"{eligible / total_eligible:.8f}" if total_eligible else "",
                             "unique_plasmids": len(plasmids)})
    write_tsv(args.output_dir / "category_summary.csv",
              ["run_status", "group", "raw_category", "cds_count", "translated_cds_count",
               "group_cds_total", "group_translated_cds_total", "cds_fraction",
               "translated_cds_fraction", "unique_plasmids"], cat_rows)

    flag_rows = []
    for group in GROUPS:
        for flag in FLAGS:
            cds, eligible, plasmids = flag_counts.get((group, flag), (0, 0, set()))
            total, total_eligible, _ = group_counts[group]
            flag_rows.append({"run_status": status, "group": group, "flag": flag,
                              "is_exclusion": flag != "original_hypothetical", "cds_count": cds,
                              "translated_cds_count": eligible, "group_cds_total": total,
                              "group_translated_cds_total": total_eligible,
                              "cds_fraction": f"{cds / total:.8f}" if total else "",
                              "translated_cds_fraction": f"{eligible / total_eligible:.8f}" if total_eligible else "",
                              "unique_plasmids": len(plasmids)})
    write_tsv(args.output_dir / "flag_summary.csv",
              ["run_status", "group", "flag", "is_exclusion", "cds_count", "translated_cds_count",
               "group_cds_total", "group_translated_cds_total", "cds_fraction",
               "translated_cds_fraction", "unique_plasmids"], flag_rows)

    write_tsv(args.output_dir / "completeness.tsv", ["accession", "expected_md5", "status"],
              [{"accession": accession, "expected_md5": expected.get(accession, ""), "status": value}
               for accession, value in sorted(statuses.items())])
    signature = next(iter(runs.values()))[3]
    provenance = json.loads(signature)
    frozen_inputs = {path.name: sha256_file(path) for path in
                     (MAP_CSV, CARGO_FASTA, BACKGROUND_FASTA, CARGO_TABLE)}
    lines = [f"status: {status}", f"expected_representative_plasmids: {len(expected)}",
             f"successful_plasmids_summarized: {len(runs)}",
             f"processed_fraction: {len(runs) / len(expected):.6f}",
             "denominator_scope: observed successful accessions only" if partial else
             "denominator_scope: complete representative plasmid set",
             f"provenance: {json.dumps(provenance, sort_keys=True)}",
             f"frozen_step5_sha256: {json.dumps(frozen_inputs, sort_keys=True)}",
             "membership: frozen step-5 FASTA headers; cargo coordinates validated against Cargo_table.csv",
             "outside_original_analysis: CDS not present in frozen cargo/background FASTAs",
             "exclusion_union: conjugation OR plasmid maintenance/replication/regulation OR PlasAnn hypothetical OR unclassified",
             "toxin-antitoxin categories are retained in the raw inventory and are not part of exclusion_union",
             "original_hypothetical is descriptive only and uses explicit original GenBank product text"]
    for group in GROUPS:
        total, eligible, plasmids = group_counts[group]
        lines.append(f"{group}: cds={total}; translated_cds={eligible}; unique_plasmids={len(plasmids)}")
    (args.output_dir / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:5]))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        sys.exit(f"ERROR: {exc}")
