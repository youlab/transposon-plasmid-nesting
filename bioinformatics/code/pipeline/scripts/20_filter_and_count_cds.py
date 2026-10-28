#!/usr/bin/env python3
"""Count original functional hits before and after the two-class PlasAnn filter.

Default input: data/19_cds_annotation_membership.csv.gz, streamed without
extracting it. Use --from-annotation to rebuild that membership from existing
per-accession annotation outputs and work/processed inputs.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from scipy.special import logsumexp
from scipy.stats import fisher_exact, hypergeom

ROOT = Path(__file__).resolve().parents[3]
RESULTS = Path(os.environ.get("NEE_RESULTS_DIR", ROOT / "work" / "processed"))
DEFAULT_INPUT = ROOT / "results" / "19_plasann_coverage_annotation"
DEFAULT_OUTPUT = ROOT / "results"
DEFAULT_MEMBERSHIP = ROOT / "data" / "19_cds_annotation_membership.csv.gz"
DEFAULT_CONFIG = ROOT / "code" / "pipeline" / "config" / "analysis.json"
THRESHOLDS = (90.0, 60.0)
DATABASES = ("card", "vfdb", "bacmet", "tadb")
LABELS = {
    "card": "Antibiotic resistance (CARD)",
    "vfdb": "Virulence (VFDB)",
    "bacmet": "Metal/biocide resistance (BacMet)",
    "tadb": "Toxin-antitoxin (TADB)",
}
STEP12_PATH = Path(__file__).with_name("18_plasann_category_helpers.py")
_SPEC = importlib.util.spec_from_file_location("step12_helpers", STEP12_PATH)
STEP12 = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(STEP12)

MODES = {
    "relaxed": {
        "prefix": "08", "cargo_table": "07_cargo_membership_relaxed.csv.gz",
        "table": "11_enrichment_relaxed.csv",
    },
    "strict": {
        "prefix": "13", "cargo_table": "07_cargo_membership_strict.csv.gz",
        "table": "13_enrichment_strict.csv",
    },
}


def md5_text(sequence):
    return hashlib.md5(sequence.encode()).hexdigest()


def classify(category):
    """Return (backbone-filter membership, reason); other labels remain retained."""
    conjugation = STEP12.norm(category) == STEP12.norm("Conjugation")
    backbone = STEP12.is_maintenance_category(category)
    reasons = []
    if conjugation:
        reasons.append("conjugation")
    if backbone:
        reasons.append("plasmid_maintenance_replication_regulation")
    return bool(reasons), ";".join(reasons)


def hit_status(protein_md5, known_hashes, hits):
    if protein_md5 not in known_hashes:
        raise ValueError(f"protein hash is missing from the original unique map: {protein_md5}")
    return protein_md5 in hits


def read_fasta_members(path, selected, all_accessions):
    """Adapt step12.read_headers, adding per-key sequence hashes for CDS audit."""
    headers = defaultdict(Counter)
    by_key_hashes = defaultdict(Counter)
    hash_counts = Counter()
    digest = hashlib.sha256()
    header, sequence = None, []

    def save_record():
        if header is None:
            return
        parts = header.split("|", 2)
        if len(parts) != 3:
            raise ValueError(f"malformed FASTA header in {path}: {header[:120]}")
        accession, locus, product = parts
        if accession not in all_accessions:
            raise ValueError(f"unknown representative accession in {path}: {accession}")
        if accession not in selected:
            return
        key = (locus or "NA", product[:60])
        protein_md5 = md5_text("".join(sequence))
        headers[accession][key] += 1
        by_key_hashes[(accession, key)][protein_md5] += 1
        hash_counts[protein_md5] += 1

    with path.open("rb") as handle:
        for raw in handle:
            digest.update(raw)
            line = raw.decode("utf-8").rstrip("\r\n")
            if line.startswith(">"):
                save_record()
                header, sequence = line[1:], []
            elif header is not None and (header.split("|", 1)[0] in selected):
                sequence.append(line.strip())
        save_record()
    return {
        "headers": headers,
        "by_key_hashes": by_key_hashes,
        "hash_counts": hash_counts,
        "sha256": digest.hexdigest(),
    }


def read_cargo_table(path, selected, accession_md5):
    coords, header_counts = defaultdict(set), defaultdict(Counter)
    required = {"accession", "locus_tag", "product", "cds_start", "cds_end"}
    digest = STEP12.sha256_file(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return _read_cargo_table_rows(reader, required, selected, accession_md5, digest)


def _read_cargo_table_rows(reader, required, selected, accession_md5, digest):
    coords, header_counts = defaultdict(set), defaultdict(Counter)
    actual = {STEP12.norm(name): name for name in (reader.fieldnames or [])}
    required = {STEP12.norm(name) for name in required}
    if not required.issubset(actual):
        raise ValueError(f"cargo table lacks columns: {sorted(required - set(actual))}")
    accession_col = actual[STEP12.norm("accession")]
    locus_col = actual[STEP12.norm("locus_tag")]
    product_col = actual[STEP12.norm("product")]
    start_col, end_col = actual[STEP12.norm("cds_start")], actual[STEP12.norm("cds_end")]
    md5_column = actual.get(STEP12.norm("md5"))
    for row in reader:
        accession = row[accession_col]
        if accession not in selected:
            continue
        if md5_column and row[md5_column] != accession_md5[accession]:
            raise ValueError(f"cargo table accession/md5 mismatch: {accession}")
        locus = row[locus_col] or "NA"
        start, end = int(row[start_col]), int(row[end_col])
        coord = (locus, start, end)
        if coord in coords[accession]:
            raise ValueError(f"duplicate cargo coordinates: {accession} {coord}")
        coords[accession].add(coord)
        header_counts[accession][(locus, row[product_col][:60])] += 1
    return coords, header_counts, digest


def read_annotations(accession, plasmid_md5, path, manifest):
    rows, coords, by_header, indices = [], {}, defaultdict(list), set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = STEP12.columns(reader.fieldnames)
        normalized = {STEP12.norm(name): name for name in (reader.fieldnames or [])}
        translation_col = normalized.get("originaltranslation")
        if translation_col is None:
            raise ValueError(f"cds_annotations.csv lacks original_translation: {path}")
        for raw in reader:
            row = {name: raw[column] or "" for name, column in columns.items()}
            row["original_translation"] = raw[translation_col] or ""
            if row["accession"] != accession or row["md5"] != plasmid_md5:
                raise ValueError(f"annotation identity mismatch: {path}")
            row["start"], row["end"] = int(row["start"]), int(row["end"])
            cds_index = int(row["cds_index"])
            if cds_index in indices:
                raise ValueError(f"duplicate CDS index: {accession} {cds_index}")
            indices.add(cds_index)
            row["eligible"] = row["protein_eligible"].casefold() in {"true", "1", "yes", "y"}
            row["locus_tag"] = row["locus_tag"] or "NA"
            row["protein_md5"] = md5_text(row["original_translation"])
            row["key"] = (row["locus_tag"], STEP12.product_for_header(row["original_product"])[:60])
            coord = (row["locus_tag"], row["start"], row["end"])
            if coord in coords:
                raise ValueError(f"duplicate annotation coordinates: {accession} {coord}")
            coords[coord] = row
            by_header[row["key"]].append(row)
            rows.append(row)
    counts = manifest["counts"]
    if len(rows) != int(counts["cds_total"]):
        raise ValueError(f"annotation row count disagrees with manifest: {accession}")
    if sum(row["eligible"] for row in rows) != int(counts["protein_eligible"]):
        raise ValueError(f"translated CDS count disagrees with manifest: {accession}")
    return rows, coords, by_header


def read_unique_map(path):
    mapping = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            protein_md5 = row["md5"]
            if protein_md5 in mapping:
                raise ValueError(f"duplicate protein hash in {path}: {protein_md5}")
            copies = int(row["n_copies"])
            if copies < 1:
                raise ValueError(f"nonpositive n_copies in {path}: {protein_md5}")
            mapping[protein_md5] = copies
    if not mapping:
        raise ValueError(f"empty unique-protein map: {path}")
    return mapping


def verify_unique_fasta(path, known_hashes):
    seen, header, sequence = set(), None, []

    def save_record():
        if header is None:
            return
        if header in seen or header not in known_hashes:
            raise ValueError(f"unexpected or duplicate unique FASTA hash in {path}: {header}")
        if md5_text("".join(sequence)) != header:
            raise ValueError(f"unique FASTA sequence does not match its protein hash: {path} {header}")
        seen.add(header)

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                save_record()
                header, sequence = line[1:], []
            elif header is not None:
                sequence.append(line.strip())
        save_record()
    if seen != set(known_hashes):
        raise ValueError(f"unique FASTA/map hash sets differ: {path}")


def read_hits(path, known_hashes):
    hits, seen_queries = set(), set()
    with path.open(errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 7:
                raise ValueError(f"malformed DIAMOND row: {path}:{line_no}")
            protein_md5 = fields[0]
            if protein_md5 not in known_hashes:
                raise ValueError(f"DIAMOND query absent from original unique map: {path}:{line_no} {protein_md5}")
            if protein_md5 in seen_queries:
                raise ValueError(f"duplicate DIAMOND query row (expected -k 1): {path}:{line_no} {protein_md5}")
            seen_queries.add(protein_md5)
            try:
                pident, qcov = float(fields[2]), float(fields[4])
            except ValueError as exc:
                raise ValueError(f"invalid DIAMOND identity/coverage: {path}:{line_no}") from exc
            if not math.isfinite(pident) or not math.isfinite(qcov):
                raise ValueError(f"nonfinite DIAMOND identity/coverage: {path}:{line_no}")
            if pident >= THRESHOLDS[0] and qcov >= THRESHOLDS[1]:
                hits.add(protein_md5)
    return hits


def legacy_inputs(mode, group):
    prefix = MODES[mode]["prefix"]
    if mode == "relaxed":
        hit_prefix = "09"
        suffix = "relaxed"
    else:
        hit_prefix = "13"
        suffix = "strict"
    def locate(name):
        for folder in (RESULTS, ROOT / "data"):
            candidate = folder / name
            if candidate.is_file():
                return candidate
        return RESULTS / name
    return {
        "map": locate(f"{prefix}_{group}_{suffix}_unique_map.csv"),
        "fasta": locate(f"{prefix}_{group}_{suffix}_unique.faa"),
        "hits": {db: locate(f"{hit_prefix}_{group}_{suffix}_vs_{db}.tsv") for db in DATABASES},
    }


def load_legacy(mode, complete):
    maps, hit_sets, totals, files = {}, {}, {}, {}
    for group in ("cargo", "background"):
        paths = legacy_inputs(mode, group)
        counts = read_unique_map(paths["map"])
        if complete:
            verify_unique_fasta(paths["fasta"], counts)
        maps[group] = counts
        totals[group] = sum(counts.values())
        hit_sets[group] = {}
        for db in DATABASES:
            path = paths["hits"][db]
            hit_sets[group][db] = read_hits(path, counts)
            files[str(path)] = STEP12.sha256_file(path)
        files[str(paths["map"])] = STEP12.sha256_file(paths["map"])
        if complete:
            files[str(paths["fasta"])] = STEP12.sha256_file(paths["fasta"])
    table_name = MODES[mode]["table"]
    table_path = next((folder / table_name for folder in (RESULTS, ROOT / "data")
                       if (folder / table_name).is_file()), RESULTS / table_name)
    table_rows = {}
    with table_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            db = row["category"].casefold()
            if db not in DATABASES or db in table_rows:
                raise ValueError(f"unexpected or duplicate category in {table_path}: {db}")
            table_rows[db] = {
                "cargo_hits": int(row["cargo_hits"]),
                "background_hits": int(row["bg_hits"]),
                "cargo_total": int(row["cargo_total"]),
                "background_total": int(row["bg_total"]),
            }
    if set(table_rows) != set(DATABASES):
        raise ValueError(f"original enrichment table does not contain all four databases: {table_path}")
    if complete:
        for group in ("cargo", "background"):
            column = f"{group}_total" if group == "cargo" else "background_total"
            if any(values[column] != totals[group] for values in table_rows.values()):
                raise ValueError(f"original map totals disagree with {table_path}: {group}")
        for db in DATABASES:
            for group in ("cargo", "background"):
                count = sum(maps[group][protein_md5] for protein_md5 in hit_sets[group][db])
                expected = table_rows[db][f"{group}_hits" if group == "cargo" else "background_hits"]
                if count != expected:
                    raise ValueError(f"original {mode} {group}/{db} weighted hit count is {count}, table says {expected}")
    files[str(table_path)] = STEP12.sha256_file(table_path)
    return maps, hit_sets, totals, table_rows, files


def resolve_mode_accession(accession, rows, coords, by_header, mode_members, cargo_table, accession_md5):
    cargo = mode_members["cargo"]
    background = mode_members["background"]
    cargo_headers = cargo["headers"].get(accession, Counter())
    if cargo_headers != cargo_table[1].get(accession, Counter()):
        raise ValueError(f"cargo FASTA headers do not match cargo table: {accession}")
    cargo_coords = cargo_table[0].get(accession, set())
    matched_headers, matched_hashes = Counter(), defaultdict(Counter)
    cargo_rows = {}
    for coord in cargo_coords:
        row = coords.get(coord)
        if row is None or not row["eligible"]:
            raise ValueError(f"cargo coordinate has no translated CDS: {accession} {coord}")
        matched_headers[row["key"]] += 1
        matched_hashes[row["key"]][row["protein_md5"]] += 1
        cargo_rows[coord] = row
    if matched_headers != cargo_table[1].get(accession, Counter()):
        raise ValueError(f"cargo table products do not match annotations: {accession}")
    for key, expected_hashes in matched_hashes.items():
        if cargo["by_key_hashes"].get((accession, key), Counter()) != expected_hashes:
            raise ValueError(f"cargo protein sequence hashes do not match annotations: {accession} {key}")
    for key, count in cargo_headers.items():
        if key not in matched_hashes:
            raise ValueError(f"cargo FASTA header has no cargo-table CDS: {accession} {key}")

    background_coords = set()
    for key, count in background["headers"].get(accession, Counter()).items():
        candidates = [row for row in by_header.get(key, [])
                      if row["eligible"] and (row["locus_tag"], row["start"], row["end"]) not in cargo_coords]
        resolved = STEP12.resolve_background(count, candidates, accession, key)
        expected_hashes = Counter(row["protein_md5"] for row in resolved)
        observed_hashes = background["by_key_hashes"].get((accession, key), Counter())
        if observed_hashes != expected_hashes:
            raise ValueError(f"background protein sequence hashes do not match annotations: {accession} {key}")
        background_coords.update((row["locus_tag"], row["start"], row["end"]) for row in resolved)
    return cargo_coords, background_coords


def _fisher_log_p(table):
    a, b, c, d = (int(value) for row in table for value in row)
    total, successes, draws = a + b + c + d, a + c, a + b
    low, high = max(0, draws - (total - successes)), min(draws, successes)
    support = range(low, high + 1)
    log_probs = hypergeom.logpmf(support, total, successes, draws)
    observed = float(hypergeom.logpmf(a, total, successes, draws))
    return min(0.0, float(logsumexp(log_probs[log_probs <= observed + 1e-12])))


def _bh_adjust_log_probabilities(log10_ps):
    order = sorted(range(len(log10_ps)), key=log10_ps.__getitem__)
    q_values, running, count = [0.0] * len(order), 0.0, len(order)
    for rank in range(count, 0, -1):
        index = order[rank - 1]
        running = min(running, log10_ps[index] + math.log(count / rank, 10.0))
        q_values[index] = min(0.0, running)
    return q_values


def _probability_text(log10_value):
    exponent = math.floor(log10_value)
    return f"{10 ** (log10_value - exponent):.6f}e{exponent:+d}"


def _csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _bool(value, field, line_number):
    text = str(value).strip().casefold()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n"}:
        return False
    raise ValueError(f"invalid boolean in {field} at membership row {line_number}: {value!r}")


def _load_config(path):
    config = json.loads(path.read_text(encoding="utf-8"))
    plasann = config.get("plasann", {})
    expected = {"identity_min_percent": 80, "query_coverage_min_percent": 90,
                "subject_coverage_min_percent": 90, "fallback_enabled": False,
                "evalue": 1e-5, "max_targets": 10}
    if any(plasann.get(key) != value for key, value in expected.items()):
        raise ValueError(f"analysis config differs from the frozen PlasAnn 80/90/90 rule: {plasann}")
    func = config.get("functional_hits", {})
    if func.get("identity_min_percent") != 90 or func.get("query_coverage_min_percent") != 60:
        raise ValueError(f"analysis config differs from frozen functional-hit 90/60 rule: {func}")
    categories = [str(value).casefold() for value in config.get("backbone_categories", [])]
    if not categories or any(value not in DATABASES for value in categories):
        raise ValueError("config backbone_categories must select one or more original databases")
    if config.get("remove_ta_positive_cds") is not False:
        raise ValueError("TA-positive CDSs must remain in the frozen denominators")
    return config, categories


def count_membership(membership_path, output_dir, config_path):
    """Stream the bundled per-CDS membership file into the two count tables."""
    config, backbone_categories = _load_config(config_path)
    opener = gzip.open if str(membership_path).endswith(".gz") else open
    totals = Counter()
    hits = Counter()
    carriers = defaultdict(set)
    filtered_totals = Counter()
    filtered_hits = Counter()
    filtered_carriers = defaultdict(set)
    states = defaultdict(lambda: {name: 0 for name in ("raw", "excluded", "retained")})
    state_carriers = defaultdict(lambda: {name: set() for name in ("raw", "excluded", "retained")})
    required = {"mode", "group", "accession", "category", "plasmid_backbone", "retained",
                *(f"{db}_hit" for db in DATABASES)}
    with opener(membership_path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        actual = {str(name).strip().casefold(): name for name in (reader.fieldnames or [])}
        missing = required - set(actual)
        if missing:
            raise ValueError(f"membership input lacks required columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, 2):
            mode = str(row[actual["mode"]]).strip().casefold()
            group = str(row[actual["group"]]).strip().casefold()
            accession = str(row[actual["accession"]]).strip()
            if mode not in {"relaxed", "strict"} or group not in {"cargo", "background"}:
                raise ValueError(f"unexpected mode/group at membership row {line_number}: {mode}/{group}")
            if not accession:
                raise ValueError(f"empty accession at membership row {line_number}")
            category = row[actual["category"]]
            expected_backbone, _reason = classify(category)
            stored_backbone = _bool(row[actual["plasmid_backbone"]], "plasmid_backbone", line_number)
            retained = _bool(row[actual["retained"]], "retained", line_number)
            if stored_backbone != expected_backbone or retained == expected_backbone:
                raise ValueError(f"stored backbone/retained flags disagree with category at membership row {line_number}: {category!r}")
            key = (mode, group)
            totals[key] += 1
            states[key]["raw"] += 1
            state_carriers[key]["raw"].add(accession)
            state = "excluded" if expected_backbone else "retained"
            states[key][state] += 1
            state_carriers[key][state].add(accession)
            for db in DATABASES:
                is_hit = _bool(row[actual[f"{db}_hit"]], f"{db}_hit", line_number)
                count_key = (mode, group, db)
                if is_hit:
                    hits[count_key] += 1
                    carriers[count_key].add(accession)
                if not expected_backbone and db in backbone_categories:
                    filtered_totals[count_key] += 1
                    if is_hit:
                        filtered_hits[count_key] += 1
                        filtered_carriers[count_key].add(accession)

    expected_pairs = {(mode, group) for mode in ("relaxed", "strict")
                      for group in ("cargo", "background")}
    if set(totals) != expected_pairs:
        raise ValueError(f"membership input does not contain all mode/group pairs: {sorted(totals)}")
    before_rows = []
    after_rows = []
    for mode, group in sorted(expected_pairs):
        for db in DATABASES:
            count_key = (mode, group, db)
            n = totals[(mode, group)]
            h = hits[count_key]
            before_rows.append({"mode": mode, "group": group, "database": db,
                                "hit_cds": h, "cds_total": n,
                                "hit_fraction": format(h / n, ".17g") if n else "",
                                "hit_carriers": len(carriers[count_key])})
            if db in backbone_categories:
                retained_n = filtered_totals[count_key]
                retained_h = filtered_hits[count_key]
                after_rows.append({"mode": mode, "group": group, "database": db,
                                   "hit_cds": retained_h, "cds_total": retained_n,
                                   "hit_fraction": format(retained_h / retained_n, ".17g") if retained_n else "",
                                   "hit_carriers": len(filtered_carriers[count_key])})

    fields = ["mode", "group", "database", "hit_cds", "cds_total", "hit_fraction", "hit_carriers"]
    output_dir.mkdir(parents=True, exist_ok=True)
    _csv(output_dir / "20_cds_counts_before_filtering.csv", fields, before_rows)
    _csv(output_dir / "20_cds_counts_after_filtering.csv", fields, after_rows)
    mode_rows = []
    for mode, group in sorted(expected_pairs):
        item = states[(mode, group)]
        if item["raw"] != item["excluded"] + item["retained"]:
            raise ValueError(f"raw/excluded/retained counts do not partition {mode}/{group}")
        mode_rows.append({
            "mode": mode, "group": group, "raw_cds_count": item["raw"],
            "excluded_cds_count": item["excluded"], "retained_cds_count": item["retained"],
            "raw_carriers": len(state_carriers[(mode, group)]["raw"]),
            "excluded_carriers": len(state_carriers[(mode, group)]["excluded"]),
            "retained_carriers": len(state_carriers[(mode, group)]["retained"]),
        })
    _csv(output_dir / "20_backbone_filter_counts.csv",
         ["mode", "group", "raw_cds_count", "excluded_cds_count", "retained_cds_count",
          "raw_carriers", "excluded_carriers", "retained_carriers"], mode_rows)
    print(f"Streamed {membership_path}; wrote {len(before_rows)} raw and {len(after_rows)} filtered count rows")


def self_test():
    assert classify("Conjugation") == (True, "conjugation")
    assert classify("PLasmid Maintenance, Replication and Regulation")[0]
    assert classify("Uncharacterized") == (False, "")
    assert classify("hypothetical protein") == (False, "")
    table_reader = csv.DictReader([
        "accession,md5,locus_tag,product,cds_start,cds_end",
        "acc,plasmid-hash,L1,enzyme,1,9",
    ])
    table_coords, table_headers, _ = _read_cargo_table_rows(
        table_reader, {"accession", "locus_tag", "product", "cds_start", "cds_end"},
        {"acc"}, {"acc": "plasmid-hash"}, "test")
    assert table_coords["acc"] == {("L1", 1, 9)} and table_headers["acc"][("L1", "enzyme")] == 1
    assert [hit_status("known", {"known"}, {"known"}), hit_status("known", {"known"}, set())] == [True, False]
    try:
        hit_status("unknown", {"known"}, set())
    except ValueError:
        pass
    else:
        raise AssertionError("unknown protein hash was treated as a no-hit")
    instances = [("same", "Conjugation"), ("same", "Uncharacterized")]
    original_hits = sum(hit_status(protein_md5, {"same"}, {"same"}) for protein_md5, _ in instances)
    retained_hits = sum(not classify(category)[0] and hit_status(protein_md5, {"same"}, {"same"})
                        for protein_md5, category in instances)
    assert original_hits == 2 and retained_hits == 1
    p = math.exp(_fisher_log_p([[1, 9], [11, 3]]))
    assert math.isclose(p, 0.0027594561852200836, rel_tol=1e-10)
    print("self-test passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--membership-file", type=Path,
                        help="per-CDS membership CSV/CSV.GZ; defaults to the configured bundle")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--from-annotation", action="store_true",
                        help="join per-accession annotations and legacy hit inputs instead of streaming the bundle")
    parser.add_argument("--limit", type=int,
                        help="pilot only: process the first N accessions and mark output PARTIAL")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.from_annotation:
        config, _ = _load_config(args.config)
        membership = args.membership_file
        if membership is None:
            configured = config.get("membership_file", "data/19_cds_annotation_membership.csv.gz")
            membership = Path(configured)
            if not membership.is_absolute():
                membership = ROOT / membership
        count_membership(membership, args.output_dir, args.config)
        return
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "run_summary.txt").write_text(
        "status: RUNNING\nanalysis: original four-database hit recount after the two-class filter\n",
        encoding="utf-8")
    (args.output_dir / "provenance.json").write_text('{"status":"RUNNING"}\n', encoding="utf-8")

    expected = STEP12.read_expected()
    selected_expected = expected
    if args.limit is not None:
        if args.limit > len(expected):
            parser.error(f"--limit exceeds representative count ({len(expected)})")
        selected_expected = dict(sorted(expected.items())[:args.limit])
    runs, statuses = STEP12.read_runs(args.input_dir, selected_expected, args.limit is not None)
    if set(runs) != set(selected_expected):
        raise ValueError("one or more selected step11 annotations are missing; pilot/full run requires every selected accession")
    selected = set(runs)
    accession_md5 = {accession: item[0] for accession, item in runs.items()}
    signature = json.loads(next(iter(runs.values()))[3])
    complete = args.limit is None

    per_cds_path = args.output_dir / "per_cds_membership.csv.gz"
    temp_path = per_cds_path.with_name(f".{per_cds_path.name}.tmp.{os.getpid()}")
    per_cds_fields = ["mode", "group", "accession", "locus_tag", "start", "end", "protein_md5",
                      "category", "plasmid_backbone", "exclusion_reason", "retained",
                      "card_hit", "vfdb_hit", "bacmet_hit", "tadb_hit"]
    group_counts = defaultdict(lambda: {"raw": 0, "excluded": 0, "retained": 0,
                                         "raw_carriers": set(), "excluded_carriers": set(),
                                         "retained_carriers": set(), "proteins": set()})
    hit_counts = Counter()
    raw_hit_counts = Counter()
    excluded_hit_counts = Counter()
    hit_carriers = defaultdict(set)
    category_counts = defaultdict(lambda: {state: 0 for state in ("raw", "excluded", "retained")})
    category_carriers = defaultdict(lambda: {state: set() for state in ("raw", "excluded", "retained")})
    cross_mode_hits = {}
    mode_files, old_totals_by_mode, original_tables, union_digests, union_counts = {}, {}, {}, {}, {}
    manifest_digest = hashlib.sha256()
    for accession in sorted(runs):
        manifest_digest.update(
            f"{accession}\t{runs[accession][2]['input_sha256']}\t{runs[accession][2]['cds_annotations_sha256']}\n".encode())

    try:
        with gzip.open(temp_path, "wt", newline="", encoding="utf-8") as compressed:
            writer = csv.DictWriter(compressed, fieldnames=per_cds_fields, lineterminator="\n")
            writer.writeheader()
            for mode in ("relaxed", "strict"):
                started_files = {}
                cargo_path = RESULTS / f"07_cargo_{mode}.faa"
                background_path = RESULTS / f"07_background_{mode}.faa"
                cargo = read_fasta_members(cargo_path, selected, set(expected))
                background = read_fasta_members(background_path, selected, set(expected))
                started_files[str(cargo_path)] = cargo["sha256"]
                started_files[str(background_path)] = background["sha256"]
                cargo_table_path = RESULTS / MODES[mode]["cargo_table"]
                if not cargo_table_path.is_file() and (ROOT / "data" / cargo_table_path.name).is_file():
                    cargo_table_path = ROOT / "data" / cargo_table_path.name
                cargo_table = read_cargo_table(cargo_table_path, selected, accession_md5)
                started_files[str(cargo_table_path)] = cargo_table[2]
                members = {"cargo": cargo, "background": background}
                maps, hit_sets, old_totals, original_table, legacy_files = load_legacy(mode, complete)
                old_totals_by_mode[mode] = old_totals
                original_tables[mode] = original_table
                started_files.update(legacy_files)
                mode_files[mode] = started_files

                if complete:
                    for group in ("cargo", "background"):
                        observed = members[group]["hash_counts"]
                        if observed != Counter(maps[group]):
                            raise ValueError(f"{mode} {group} frozen FASTA protein multiplicities differ from original unique map")

                union_digest = hashlib.sha256()
                mode_union_count = 0
                for accession in sorted(runs):
                    plasmid_md5, annotation_path, manifest, _ = runs[accession]
                    rows, coords, by_header = read_annotations(accession, plasmid_md5, annotation_path, manifest)
                    cargo_coords, background_coords = resolve_mode_accession(
                        accession, rows, coords, by_header, members, cargo_table, accession_md5)
                    observed_count = len(cargo_coords) + len(background_coords)
                    expected_count = sum(members[group]["headers"].get(accession, Counter()).total()
                                         for group in ("cargo", "background"))
                    if observed_count != expected_count:
                        raise ValueError(f"{mode} membership rows do not account for every FASTA instance: {accession}")
                    group_coords = {"cargo": cargo_coords, "background": background_coords}
                    for row in rows:
                        coord = (row["locus_tag"], row["start"], row["end"])
                        group = next((name for name, coords_set in group_coords.items() if coord in coords_set), None)
                        if group is None:
                            continue
                        protein_md5 = row["protein_md5"]
                        if protein_md5 not in maps[group]:
                            raise ValueError(f"frozen protein hash is absent from original unique map: {mode} {group} {accession} {coord} {protein_md5}")
                        union_digest.update(
                            f"{accession}\t{row['locus_tag']}\t{row['start']}\t{row['end']}\t{protein_md5}\n".encode())
                        mode_union_count += 1
                        excluded, reason = classify(row["category"])
                        retained = not excluded
                        state = "retained" if retained else "excluded"
                        stats = group_counts[(mode, group)]
                        stats["raw"] += 1
                        stats["raw_carriers"].add(accession)
                        stats["proteins"].add(protein_md5)
                        stats[state] += 1
                        stats[f"{state}_carriers"].add(accession)
                        category_key = (mode, group, row["category"])
                        category_counts[category_key]["raw"] += 1
                        category_counts[category_key][state] += 1
                        category_carriers[category_key]["raw"].add(accession)
                        category_carriers[category_key][state].add(accession)

                        status_tuple = tuple(hit_status(protein_md5, maps[group], hit_sets[group][db])
                                             for db in DATABASES)
                        previous = cross_mode_hits.setdefault(protein_md5, status_tuple)
                        if previous != status_tuple:
                            raise ValueError(f"same protein hash has inconsistent original database hit status across sets: {protein_md5}")
                        for db, is_hit in zip(DATABASES, status_tuple):
                            if is_hit:
                                raw_hit_counts[(mode, group, db)] += 1
                                if not retained:
                                    excluded_hit_counts[(mode, group, db)] += 1
                            if is_hit and retained:
                                hit_counts[(mode, group, db)] += 1
                                hit_carriers[(mode, group, db)].add(accession)
                        writer.writerow({
                            "mode": mode, "group": group, "accession": accession,
                            "locus_tag": row["locus_tag"], "start": row["start"], "end": row["end"],
                            "protein_md5": protein_md5, "category": row["category"],
                            "plasmid_backbone": str(excluded).lower(), "exclusion_reason": reason,
                            "retained": str(retained).lower(),
                            **{f"{db}_hit": str(is_hit).lower() for db, is_hit in zip(DATABASES, status_tuple)},
                        })
                union_digests[mode] = union_digest.hexdigest()
                union_counts[mode] = mode_union_count
                print(f"{mode}: joined {mode_union_count:,} CDS instances across {len(runs):,} accessions", flush=True)
        os.replace(temp_path, per_cds_path)
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise

    if complete:
        if union_counts != {"relaxed": 3_323_101, "strict": 3_323_101}:
            raise ValueError(f"cargo/background CDS union does not match frozen all-mode total: {union_counts}")
        if union_digests["relaxed"] != union_digests["strict"]:
            raise ValueError("relaxed and strict frozen groups do not contain the same CDS instances")
        for mode in ("relaxed", "strict"):
            for group in ("cargo", "background"):
                if group_counts[(mode, group)]["raw"] != old_totals_by_mode[mode][group]:
                    raise ValueError(f"{mode} {group} CDS total differs from original weighted analysis: {group_counts[(mode, group)]['raw']} vs {old_totals_by_mode[mode][group]}")
                for db in DATABASES:
                    expected_hits = original_tables[mode][db]["cargo_hits" if group == "cargo" else "background_hits"]
                    observed_hits = raw_hit_counts[(mode, group, db)]
                    if observed_hits != expected_hits:
                        raise ValueError(f"{mode} {group}/{db} per-CDS hit count is {observed_hits}, original table says {expected_hits}")
        for state in ("excluded", "retained"):
            relaxed_count = sum(group_counts[("relaxed", group)][state] for group in ("cargo", "background"))
            strict_count = sum(group_counts[("strict", group)][state] for group in ("cargo", "background"))
            if relaxed_count != strict_count:
                raise ValueError(f"relaxed/strict {state} CDS totals differ despite matching frozen union")

    mode_rows = []
    for mode in ("relaxed", "strict"):
        for group in ("cargo", "background"):
            item = group_counts[(mode, group)]
            if item["raw"] != item["excluded"] + item["retained"]:
                raise ValueError(f"raw/excluded/retained counts do not partition {mode}/{group}")
            mode_rows.append({
                "mode": mode, "group": group, "raw_cds_count": item["raw"],
                "excluded_cds_count": item["excluded"], "retained_cds_count": item["retained"],
                "raw_unique_proteins": len(item["proteins"]),
                "raw_carriers": len(item["raw_carriers"]),
                "excluded_carriers": len(item["excluded_carriers"]),
                "retained_carriers": len(item["retained_carriers"]),
            })
    _csv(args.output_dir / "mode_counts.csv",
         ["mode", "group", "raw_cds_count", "excluded_cds_count", "retained_cds_count",
          "raw_unique_proteins", "raw_carriers", "excluded_carriers", "retained_carriers"], mode_rows)

    hit_rows, stats_rows = [], []
    for mode in ("relaxed", "strict"):
        log10_ps = []
        mode_stats = []
        for db in DATABASES:
            cargo_total = group_counts[(mode, "cargo")]["retained"]
            background_total = group_counts[(mode, "background")]["retained"]
            if cargo_total == 0 or background_total == 0:
                raise ValueError(f"cannot calculate retained hit fractions with an empty group: {mode} cargo={cargo_total}, background={background_total}")
            cargo_hit = hit_counts.get((mode, "cargo", db), 0)
            background_hit = hit_counts.get((mode, "background", db), 0)
            cargo_raw = group_counts[(mode, "cargo")]["raw"]
            background_raw = group_counts[(mode, "background")]["raw"]
            cargo_raw_hit = raw_hit_counts[(mode, "cargo", db)]
            background_raw_hit = raw_hit_counts[(mode, "background", db)]
            cargo_excluded_hit = excluded_hit_counts[(mode, "cargo", db)]
            background_excluded_hit = excluded_hit_counts[(mode, "background", db)]
            hit_rows.extend((
                {"mode": mode, "group": "cargo", "database": db, "retained_cds_total": cargo_total,
                 "retained_hit_cds": cargo_hit, "retained_hit_fraction": f"{cargo_hit / cargo_total:.12g}",
                 "retained_hit_carriers": len(hit_carriers[(mode, "cargo", db)]),
                 "raw_cds_total": cargo_raw, "raw_hit_cds": cargo_raw_hit,
                 "excluded_hit_cds": cargo_excluded_hit},
                {"mode": mode, "group": "background", "database": db, "retained_cds_total": background_total,
                 "retained_hit_cds": background_hit, "retained_hit_fraction": f"{background_hit / background_total:.12g}",
                 "retained_hit_carriers": len(hit_carriers[(mode, "background", db)]),
                 "raw_cds_total": background_raw, "raw_hit_cds": background_raw_hit,
                 "excluded_hit_cds": background_excluded_hit},
            ))
            table = [[cargo_hit, cargo_total - cargo_hit],
                     [background_hit, background_total - background_hit]]
            odds = float(fisher_exact(table, alternative="two-sided").statistic)
            log10_p = _fisher_log_p(table) / math.log(10.0)
            log10_ps.append(log10_p)
            mode_stats.append({
                "mode": mode, "database": db, "label": LABELS[db],
                "cargo_hit_cds": cargo_hit, "cargo_other_retained_cds": cargo_total - cargo_hit,
                "cargo_retained_cds_total": cargo_total, "cargo_hit_fraction": f"{cargo_hit / cargo_total:.12g}",
                "background_hit_cds": background_hit,
                "background_other_retained_cds": background_total - background_hit,
                "background_retained_cds_total": background_total,
                "background_hit_fraction": f"{background_hit / background_total:.12g}",
                "odds_ratio_cargo_vs_background": f"{odds:.12g}",
                "p_value_two_sided": _probability_text(log10_p), "log10_p": f"{log10_p:.12g}",
            })
        for row, log10_q in zip(mode_stats, _bh_adjust_log_probabilities(log10_ps)):
            row["q_value_BH"] = _probability_text(log10_q)
            row["log10_q"] = f"{log10_q:.12g}"
        stats_rows.extend(mode_stats)
    _csv(args.output_dir / "retained_hits_by_database.csv",
         ["mode", "group", "database", "retained_cds_total", "retained_hit_cds",
          "retained_hit_fraction", "retained_hit_carriers", "raw_cds_total", "raw_hit_cds",
          "excluded_hit_cds"], hit_rows)
    _csv(args.output_dir / "fisher_by_mode.csv",
         ["mode", "database", "label", "cargo_hit_cds", "cargo_other_retained_cds",
          "cargo_retained_cds_total", "cargo_hit_fraction", "background_hit_cds",
          "background_other_retained_cds", "background_retained_cds_total", "background_hit_fraction",
          "odds_ratio_cargo_vs_background", "p_value_two_sided", "log10_p", "q_value_BH", "log10_q"],
         stats_rows)

    category_rows = []
    for (mode, group, category), values in sorted(category_counts.items()):
        carrier_values = category_carriers[(mode, group, category)]
        category_rows.append({
            "mode": mode, "group": group, "category": category,
            "raw_cds_count": values["raw"], "excluded_cds_count": values["excluded"],
            "retained_cds_count": values["retained"],
            "raw_carriers": len(carrier_values["raw"]),
            "excluded_carriers": len(carrier_values["excluded"]),
            "retained_carriers": len(carrier_values["retained"]),
        })
    _csv(args.output_dir / "plasann_category_counts.csv",
         ["mode", "group", "category", "raw_cds_count", "excluded_cds_count", "retained_cds_count",
          "raw_carriers", "excluded_carriers", "retained_carriers"], category_rows)

    count_membership(per_cds_path, args.output_dir, args.config)

    if complete:
        # Both step-5 modes partition the same CDS universe; identical hashes must
        # also produce identical old DIAMOND hit status wherever they recur.
        for db in DATABASES:
            relaxed_raw = sum(raw_hit_counts[("relaxed", group, db)] for group in ("cargo", "background"))
            strict_raw = sum(raw_hit_counts[("strict", group, db)] for group in ("cargo", "background"))
            relaxed_retained = sum(hit_counts[("relaxed", group, db)] for group in ("cargo", "background"))
            strict_retained = sum(hit_counts[("strict", group, db)] for group in ("cargo", "background"))
            relaxed_excluded = sum(excluded_hit_counts[("relaxed", group, db)] for group in ("cargo", "background"))
            strict_excluded = sum(excluded_hit_counts[("strict", group, db)] for group in ("cargo", "background"))
            if (relaxed_raw, relaxed_excluded, relaxed_retained) != (strict_raw, strict_excluded, strict_retained):
                raise ValueError(f"raw/excluded/retained {db} hits differ across frozen modes")

    status = "PARTIAL" if not complete else "COMPLETE"
    summary = [
        f"status: {status}",
        "analysis: original CARD/VFDB/BacMet/TADB hit status (pident>=90, qcov>=60) after excluding only Conjugation and step12.is_maintenance_category matches",
        "hypothetical and unclassified PlasAnn records are retained unless they match one of the two explicit excluded classes",
        f"expected_representative_plasmids: {len(expected)}",
        f"processed_representative_plasmids: {len(runs)}",
        f"step11_provenance: {json.dumps(signature, sort_keys=True)}",
        "old_step7_step8_weighted_count_validation: SKIPPED for --limit pilot" if not complete else
        "old_step7_step8_weighted_count_validation: PASSED for all four databases in both modes",
        f"relaxed_strict_CDS_union: {union_counts}; matching_sha256={union_digests['relaxed'] == union_digests['strict']}",
        "protein_md5: md5(exact original_translation UTF-8 bytes); row md5 is the plasmid sequence checksum",
        "zero-hit status is assigned only for hashes present in the matching old unique-protein map",
        f"per_CDS_audit: {per_cds_path}",
    ]
    for row in mode_rows:
        summary.append(
            f"{row['mode']}/{row['group']}: raw={row['raw_cds_count']}; excluded={row['excluded_cds_count']}; retained={row['retained_cds_count']}; "
            f"raw_carriers={row['raw_carriers']}; excluded_carriers={row['excluded_carriers']}; retained_carriers={row['retained_carriers']}"
        )
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "script": str(Path(__file__).resolve()),
        "script_sha256": STEP12.sha256_file(Path(__file__)),
        "step12_helper_script": str(STEP12_PATH.resolve()),
        "step12_helper_script_sha256": STEP12.sha256_file(STEP12_PATH),
        "input_step11_dir": str(args.input_dir.resolve()),
        "validated_step11_manifests": len(runs),
        "validated_step11_manifest_inputs_sha256": manifest_digest.hexdigest(),
        "plasann_provenance": signature,
        "protein_md5_definition": "md5(exact original_translation UTF-8 bytes)",
        "identity_threshold_percent": THRESHOLDS[0],
        "query_coverage_threshold_percent": THRESHOLDS[1],
        "mode_inputs": mode_files,
        "union_cds_counts": union_counts,
        "union_membership_sha256": union_digests,
        "limit": args.limit,
    }
    (args.output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (args.output_dir / "run_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary), flush=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, csv.Error, json.JSONDecodeError) as exc:
        sys.exit(f"ERROR: {exc}")
