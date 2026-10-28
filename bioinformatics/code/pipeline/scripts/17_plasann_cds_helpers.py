#!/usr/bin/env python3
"""Run the pinned PlasAnn CDS classifier on the frozen GBFF set.

Inputs: work/raw/gbff and work/raw/PlasAnn_database, plus the sequence map in
work/processed (or the bundled data fallback). Outputs: per-accession CDS
annotations under results/17_plasann_cds_annotations. Use --self-test or
--prepare-db before a full run. This is the original identity-tier core;
coverage-qualified results are produced by 19_annotate_plasann_cds.py.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import importlib.metadata
import json
import os
import shutil
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PROCESSED_DIR = Path(os.environ.get("NEE_RESULTS_DIR", ROOT / "work" / "processed"))
RAW_DIR = Path(os.environ.get("NEE_DATA_DIR", ROOT / "work" / "raw"))
MAP_CSV = (PROCESSED_DIR / "03_plasmid_sequence_hashes.csv"
           if (PROCESSED_DIR / "03_plasmid_sequence_hashes.csv").is_file()
           else ROOT / "data" / "03_plasmid_sequence_hashes.csv")
GBFF_DIR = RAW_DIR / "gbff"
DATABASE_DIR = RAW_DIR / "PlasAnn_database"
DATABASE_CSV = DATABASE_DIR / "Database.csv"
BLAST_DIR = DATABASE_DIR / "database_blast"
BLAST_PREFIX = BLAST_DIR / "translations_db"
CACHE_MANIFEST = BLAST_DIR / "plasann_cache_manifest.json"
DEFAULT_OUTPUT = ROOT / "results" / "17_plasann_cds_annotations"
PLASANN_VERSION = "1.1.6"
PLASANN_COMMIT = "35199e00d90462b9cf84e9f924ef098edd5da122"
RESULT_COLUMNS = [
    "accession", "md5", "cds_index", "feature_index", "locus_tag",
    "original_product", "original_translation", "location", "start", "end",
    "strand", "codon_start", "transl_table", "protein_eligible",
    "Gene Name", "Product", "Category", "plasann_translation",
]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_sha256(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with open(tmp, "w") as f:
        json.dump(value, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


def software_modules():
    env_bin = str(Path(sys.executable).resolve().parent)
    os.environ["PATH"] = env_bin + os.pathsep + os.environ.get("PATH", "")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(name, "1")
    import pandas as pd
    from Bio import SeqIO
    from plasann import blast_search, database_download
    from plasann.essential_annotation import generate_orf_annotation

    installed_version = importlib.metadata.version("plasann")
    if installed_version != PLASANN_VERSION:
        raise RuntimeError(
            f"Expected PlasAnn {PLASANN_VERSION}, found {installed_version}"
        )
    return pd, SeqIO, blast_search, database_download, generate_orf_annotation


def reference_sha256():
    if not DATABASE_CSV.is_file():
        raise FileNotFoundError(f"Local PlasAnn reference CSV not found: {DATABASE_CSV}")
    digest = sha256_file(DATABASE_CSV)
    source_manifest = DATABASE_DIR / "database_manifest.json"
    if source_manifest.is_file():
        recorded = json.loads(source_manifest.read_text()).get("sha256")
        if recorded and recorded != digest:
            raise RuntimeError("Database.csv does not match database_manifest.json")
    return digest


def blast_components(blast_dir=BLAST_DIR):
    paths = [blast_dir / "translations.fasta"]
    paths.extend(sorted(blast_dir.glob("translations_db.*")))
    return {p.name: sha256_file(p) for p in paths if p.is_file()}


def prepare_database():
    """Build/check the official local BLAST cache, never downloading a database."""
    _, _, _, database_download, _ = software_modules()
    source_hash = reference_sha256()
    BLAST_DIR.mkdir(parents=True, exist_ok=True)
    expected = {
        "database_sha256": source_hash,
        "plasann_version": PLASANN_VERSION,
        "plasann_commit": PLASANN_COMMIT,
        "blast_prefix": str(BLAST_PREFIX),
    }
    required = [BLAST_DIR / f"translations_db.{ext}" for ext in ("pin", "phr", "psq")]
    try:
        cached = json.loads(CACHE_MANIFEST.read_text())
        if all(cached.get(k) == v for k, v in expected.items()) and all(p.is_file() for p in required):
            if cached.get("components") == blast_components():
                print(f"PlasAnn BLAST cache is current: {BLAST_PREFIX}")
                return
    except (OSError, json.JSONDecodeError):
        pass

    if not shutil.which("makeblastdb"):
        raise RuntimeError("makeblastdb is not on PATH; activate conda_envs/plasann or invoke its bin/python")
    # Build in a sibling directory so a failed or stale-cache rebuild cannot
    # remove the previously usable index. The official helper writes ./database_blast.
    stage = Path(tempfile.mkdtemp(prefix=".plasann_prepare_", dir=DATABASE_DIR))
    previous = Path.cwd()
    try:
        database_download.DATABASE_DIR = DATABASE_DIR
        os.chdir(stage)
        database_download.prepare_blast_database()
        staged_blast_dir = stage / "database_blast"
        if not all((staged_blast_dir / f"translations_db.{ext}").is_file() for ext in ("pin", "phr", "psq")):
            raise RuntimeError("Official PlasAnn builder did not create the protein BLAST index")
        manifest = dict(expected, components=blast_components(staged_blast_dir))
        write_json(staged_blast_dir / CACHE_MANIFEST.name, manifest)
        backup = DATABASE_DIR / f"database_blast.old.{os.getpid()}"
        if BLAST_DIR.exists():
            os.replace(BLAST_DIR, backup)
        try:
            os.replace(staged_blast_dir, BLAST_DIR)
        except Exception:
            if backup.exists() and not BLAST_DIR.exists():
                os.replace(backup, BLAST_DIR)
            raise
        if backup.exists():
            shutil.rmtree(backup)
    finally:
        os.chdir(previous)
        shutil.rmtree(stage, ignore_errors=True)
    print(f"Prepared checksum-bound PlasAnn BLAST cache: {BLAST_PREFIX}")


def require_database():
    source_hash = reference_sha256()
    cached = json.loads(CACHE_MANIFEST.read_text()) if CACHE_MANIFEST.is_file() else {}
    if cached.get("database_sha256") != source_hash or cached.get("plasann_commit") != PLASANN_COMMIT:
        raise RuntimeError("PlasAnn BLAST cache is absent or stale; run with --prepare-db first")
    if cached.get("components") != blast_components():
        raise RuntimeError("PlasAnn BLAST cache files changed; run with --prepare-db first")
    if not all((BLAST_DIR / f"translations_db.{ext}").is_file() for ext in ("pin", "phr", "psq")):
        raise RuntimeError("PlasAnn protein BLAST index is incomplete; run with --prepare-db first")
    return source_hash


def read_mapping():
    mapping = []
    seen_md5, seen_accessions = set(), set()
    with open(MAP_CSV, newline="") as f:
        for row in csv.DictReader(f):
            md5, accession = row["md5"], row["representative_acc"]
            if md5 in seen_md5 or accession in seen_accessions:
                raise ValueError(f"Non-unique sequence map entry: {md5}, {accession}")
            seen_md5.add(md5)
            seen_accessions.add(accession)
            mapping.append((md5, accession))
    return mapping


def select_shard(mapping, shard_index, shard_count):
    """Partition the frozen representative order, before any cache filtering."""
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("Require shard_count >= 1 and 0 <= shard_index < shard_count")
    return mapping[shard_index::shard_count]


def extract_cds(record, accession, md5):
    """Extract every CDS in feature order; Biopython handles strand/joins."""
    rows, sequences = [], []
    for feature_index, feature in enumerate(record.features):
        if feature.type != "CDS":
            continue
        sequence = str(feature.extract(record.seq)).upper()
        if not sequence:
            raise ValueError(f"Empty extracted CDS at feature {feature_index}")
        q = feature.qualifiers
        translation = q.get("translation", [""])[0]
        strand = feature.location.strand
        rows.append({
            "accession": accession,
            "md5": md5,
            "cds_index": len(rows),
            "feature_index": feature_index,
            "locus_tag": q.get("locus_tag", ["NA"])[0],
            "original_product": q.get("product", ["hypothetical protein"])[0],
            "original_translation": translation,
            "location": str(feature.location),
            "start": int(feature.location.start) + 1,
            "end": int(feature.location.end),
            "strand": "+" if strand == 1 else "-" if strand == -1 else "?",
            "codon_start": q.get("codon_start", [""])[0],
            "transl_table": q.get("transl_table", [""])[0],
            "protein_eligible": bool(str(translation).strip()),
        })
        sequences.append(sequence)
    return rows, sequences


def config_for():
    actual_sha = sha256_file(Path(__file__))
    compatibility_path = Path(__file__).with_name("11_PlasAnn_source_compatibility.json")
    compatibility = json.loads(compatibility_path.read_text()) if compatibility_path.exists() else {}
    return {
        "source": "PlasAnn CDS BLASTX core",
        "script_revision": "11_PlasAnn_annotation.py/1",
        "plasann_version": PLASANN_VERSION,
        "plasann_commit": PLASANN_COMMIT,
        # Only reviewed, exact-source cleanup fixes reuse an annotation cache.
        # Unlisted future edits keep their own SHA and invalidate cached results.
        "script_sha256": compatibility.get("compatible_sources", {}).get(actual_sha, actual_sha),
        "blast_evalue": "1e-5",
        "blast_max_target_seqs": 10,
        "identity_majority_thresholds_percent": [80, 60],
    }


def save_annotations(csv_path, rows):
    tmp = csv_path.with_name(csv_path.name + f".tmp.{os.getpid()}")
    with open(tmp, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=RESULT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, csv_path)


def load_csv_counts(csv_path):
    categories, eligible_categories = Counter(), Counter()
    rows = 0
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            rows += 1
            category = row.get("Category", "")
            categories[category] += 1
            if row.get("protein_eligible", "").lower() == "true":
                eligible_categories[category] += 1
    return rows, categories, eligible_categories


def cleanup_temp(temp_dir, log_event):
    """Temporary-file cleanup cannot invalidate committed CDS annotations."""
    try:
        shutil.rmtree(temp_dir)
    except OSError as exc:
        log_event("cleanup_warning", error=str(exc), blast_temp_dir=str(temp_dir))
        print(f"WARNING: annotation saved; temporary cleanup deferred: {temp_dir}: {exc}", flush=True)


def process_one(md5, accession, reference_hash, config, config_hash, output, blast_search, generate_orf_annotation, pd, SeqIO):
    start_time = time.time()
    out_dir = output / "plasmids" / accession
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "cds_annotations.csv"
    success_path = out_dir / "success.json"
    failure_path = out_dir / "failure.json"
    log_path = out_dir / "run.log"
    gbff_path = GBFF_DIR / f"{accession}.gbff"
    temp_dir = None
    def log_event(event, **values):
        with open(log_path, "a") as log:
            log.write(json.dumps(dict(time=time.time(), event=event, accession=accession, **values), sort_keys=True) + "\n")

    try:
        log_event("begin")
        input_hash = sha256_file(gbff_path)
        expected = {
            "accession": accession,
            "md5": md5,
            "input_sha256": input_hash,
            "reference_sha256": reference_hash,
            "plasann_version": PLASANN_VERSION,
            "plasann_commit": PLASANN_COMMIT,
            "config_sha256": config_hash,
        }
        if success_path.exists():
            old = json.loads(success_path.read_text())
            if all(old.get(k) == v for k, v in expected.items()) and csv_path.is_file() and old.get("cds_annotations_sha256") == sha256_file(csv_path):
                n_rows, category_counts, eligible_counts = load_csv_counts(csv_path)
                if n_rows == old.get("counts", {}).get("cds_total"):
                    failure_path.unlink(missing_ok=True)
                    log_event("cached", cds_total=n_rows)
                    return accession, old["counts"], True
            raise RuntimeError(f"Existing completed result has different identity/configuration; choose a new --output: {out_dir}")

        record = SeqIO.read(gbff_path, "genbank")
        observed_md5 = hashlib.md5(str(record.seq).upper().encode("ascii")).hexdigest()
        if observed_md5 != md5:
            raise ValueError(f"GBFF sequence MD5 mismatch: expected {md5}, observed {observed_md5}")
        base_rows, sequences = extract_cds(record, accession, md5)
        if sequences:
            cds_df = pd.DataFrame({
                "Sequence": sequences,
                "Start": [r["start"] for r in base_rows],
                "End": [r["end"] for r in base_rows],
                "Strand": [r["strand"] for r in base_rows],
            })
            temp_dir = Path(tempfile.mkdtemp(prefix="plasann_blast_", dir=out_dir))
            blast_search.perform_blast_multiprocessing(cds_df, str(BLAST_PREFIX), temp_dir, max_workers=1)
            files = list(temp_dir.glob("blast_result_*.csv"))
            if len(files) != len(sequences):
                raise RuntimeError(f"BLAST produced {len(files)} per-CDS result files for {len(sequences)} CDS features")
            blast_search.annotate_blast_results(temp_dir, str(DATABASE_CSV))
            required = {"identity", "subject_id", "Gene Name", "Product", "Category"}
            for i in range(len(sequences)):
                result_path = temp_dir / f"blast_result_{i}.csv"
                if not result_path.is_file():
                    raise RuntimeError(f"Missing expected BLAST result: {result_path.name}")
                if not required.issubset(pd.read_csv(result_path, nrows=0).columns):
                    raise RuntimeError(f"BLAST result was not enriched with required PlasAnn columns: {result_path.name}")
            annotated = generate_orf_annotation(cds_df, temp_dir)
            if len(annotated) != len(base_rows) or not {"Gene Name", "Product", "Category", "Translation"}.issubset(annotated.columns):
                raise RuntimeError("PlasAnn core returned incomplete CDS annotations")
        else:
            annotated = pd.DataFrame(columns=["Gene Name", "Product", "Category", "Translation"])

        rows = []
        for original, prediction in zip(base_rows, annotated.to_dict(orient="records")):
            row = dict(original)
            for key in ("Gene Name", "Product", "Category"):
                value = prediction.get(key, "")
                row[key] = "" if pd.isna(value) else value
            value = prediction.get("Translation", "")
            row["plasann_translation"] = "" if pd.isna(value) else value
            rows.append(row)
        save_annotations(csv_path, rows)
        n_rows, category_counts, eligible_counts = load_csv_counts(csv_path)
        counts = {
            "cds_total": n_rows,
            "protein_eligible": sum(bool(r["protein_eligible"]) for r in base_rows),
            "category_counts": dict(category_counts),
            "protein_eligible_category_counts": dict(eligible_counts),
        }
        manifest = dict(expected, status="success", config=config, counts=counts,
                        runner_source_sha256=sha256_file(Path(__file__)),
                        cds_annotations_sha256=sha256_file(csv_path),
                        timing_seconds=round(time.time() - start_time, 3))
        write_json(success_path, manifest)
        failure_path.unlink(missing_ok=True)
        if temp_dir is not None:
            cleanup_temp(temp_dir, log_event)
        log_event("success", cds_total=n_rows, protein_eligible=counts["protein_eligible"], seconds=manifest["timing_seconds"])
        print(f"{accession}: {n_rows} CDS, {counts['protein_eligible']} protein-eligible", flush=True)
        return accession, counts, False
    except Exception as exc:
        failure = {
            "status": "failed", "accession": accession, "md5": md5,
            "error": f"{type(exc).__name__}: {exc}",
            "timing_seconds": round(time.time() - start_time, 3),
            "blast_temp_dir": str(temp_dir) if temp_dir else None,
        }
        write_json(failure_path, failure)
        log_event("error", error=failure["error"], blast_temp_dir=failure["blast_temp_dir"])
        print(f"FAILED {accession}: {failure['error']}", flush=True)
        raise


def self_test():
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    _, SeqIO, _, _, _ = software_modules()
    del SeqIO
    record = SeqRecord(Seq("AAACCCGGGTTT"), id="synthetic")
    record.features = [
        SeqFeature(CompoundLocation([FeatureLocation(0, 3, strand=1), FeatureLocation(6, 9, strand=1)]), type="CDS", qualifiers={"locus_tag": ["plus"], "translation": ["KG"]}),
        SeqFeature(FeatureLocation(0, 3, strand=-1), type="CDS", qualifiers={"locus_tag": ["minus"]}),
    ]
    digest = hashlib.md5(str(record.seq).upper().encode()).hexdigest()
    rows, seqs = extract_cds(record, "synthetic", digest)
    assert seqs == ["AAAGGG", "TTT"]
    assert [(r["cds_index"], r["feature_index"], r["locus_tag"]) for r in rows] == [(0, 0, "plus"), (1, 1, "minus")]
    assert rows[0]["protein_eligible"] is True and rows[1]["protein_eligible"] is False
    assert all(r["md5"] == digest for r in rows)
    from unittest.mock import patch
    events = []
    with patch.object(shutil, "rmtree", side_effect=OSError(39, "Directory not empty")):
        cleanup_temp(Path("synthetic_temp"), lambda event, **values: events.append(event))
    assert events == ["cleanup_warning"]
    sample = list(range(25))
    shards = [select_shard(sample, i, 3) for i in range(3)]
    assert sorted(value for shard in shards for value in shard) == sample
    assert all(set(shards[i]).isdisjoint(shards[j]) for i in range(3) for j in range(i))
    print("Self-test passed: compound/strand extraction, row order, eligibility, and identity.")


def write_counts(output, completed):
    categories, eligible = Counter(), Counter()
    for counts in completed.values():
        categories.update(counts["category_counts"])
        eligible.update(counts["protein_eligible_category_counts"])
    path = output / "category_counts.csv"
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with open(tmp, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Category", "cds_count", "protein_eligible_count"])
        for category in sorted(categories):
            writer.writerow([category, categories[category], eligible[category]])
    os.replace(tmp, path)
    return categories, eligible


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=4, help="Concurrent plasmids (each official BLASTX uses up to 4 threads)")
    parser.add_argument("--shard-index", type=int, default=0, help="Zero-based shard index")
    parser.add_argument("--shard-count", type=int, default=1, help="Number of disjoint shards; resume uses the same shard count")
    parser.add_argument("--limit", type=int, help="Process the first N mapped representatives")
    parser.add_argument("--accession", action="append", help="Process this representative accession; repeatable")
    parser.add_argument("--prepare-db", action="store_true", help="Build or validate the local official BLAST database, then exit")
    parser.add_argument("--self-test", action="store_true", help="Run synthetic CDS extraction and identity checks, then exit")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.prepare_db:
        prepare_database()
        return 0
    if args.workers < 1 or (args.limit is not None and args.limit < 1):
        parser.error("--workers and --limit must be positive")
    if args.shard_count < 1 or not 0 <= args.shard_index < args.shard_count:
        parser.error("Require --shard-count >= 1 and 0 <= --shard-index < --shard-count")

    reference_hash = require_database()
    config = config_for()
    config_hash = json_sha256(config)
    full_mapping = read_mapping()
    mapping = select_shard(full_mapping, args.shard_index, args.shard_count)
    if args.accession:
        requested = set(args.accession)
        known = {accession for _, accession in full_mapping}
        unknown = requested - known
        if unknown:
            parser.error("Unknown representative accession(s): " + ", ".join(sorted(unknown)))
        mapping = [item for item in mapping if item[1] in requested]
    elif args.limit is not None:
        mapping = mapping[:args.limit]

    args.output.mkdir(parents=True, exist_ok=True)
    pd, SeqIO, blast_search, _, generate_orf_annotation = software_modules()
    completed, failures = {}, []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(process_one, md5, accession, reference_hash, config,
                        config_hash, args.output, blast_search,
                        generate_orf_annotation, pd, SeqIO): accession
            for md5, accession in mapping
        }
        for future in concurrent.futures.as_completed(futures):
            accession = futures[future]
            try:
                _, counts, _ = future.result()
                completed[accession] = counts
            except Exception as exc:
                failures.append({"accession": accession, "error": f"{type(exc).__name__}: {exc}"})

    summary_output = args.output
    if args.shard_count > 1:
        summary_output = args.output / "shards" / f"{args.shard_index}_of_{args.shard_count}"
        summary_output.mkdir(parents=True, exist_ok=True)
    categories, _ = write_counts(summary_output, completed)
    summary = {
        "runner_source_sha256": sha256_file(Path(__file__)),
        "status": "failed" if failures else "success",
        "reference_sha256": reference_hash,
        "plasann_version": PLASANN_VERSION,
        "plasann_commit": PLASANN_COMMIT,
        "config": config,
        "config_sha256": config_hash,
        "workers": args.workers,
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "plasmids_requested": len(mapping),
        "plasmids_completed": len(completed),
        "failures": failures,
        "cds_total": sum(c["cds_total"] for c in completed.values()),
        "protein_eligible_total": sum(c["protein_eligible"] for c in completed.values()),
        "category_counts": dict(categories),
    }
    write_json(summary_output / "run_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
