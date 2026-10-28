#!/usr/bin/env python3
"""Annotate frozen CDS features with PlasAnn BLASTX 80/90/90 coverage support.

Inputs: work/raw/gbff, work/raw/PlasAnn_database, and the sequence map in
work/processed (falling back to data/03_plasmid_sequence_hashes.csv).
Output: per-accession annotation/evidence files under results/19_plasann_coverage_annotation.
Run with the pinned PlasAnn environment; use --self-test for synthetic checks.
"""
import argparse
import concurrent.futures
import csv
import fcntl
import gzip
import hashlib
import importlib.util
import inspect
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CONFIG_PATH = ROOT / "code" / "pipeline" / "config" / "analysis.json"
HELPER_PATH = Path(__file__).with_name("17_plasann_cds_helpers.py")
_spec = importlib.util.spec_from_file_location("plasann_cds_helpers", HELPER_PATH)
HELPER = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(HELPER)

DEFAULT_OUTPUT = ROOT / "results" / "19_plasann_coverage_annotation"
EXPECTED_REFERENCE_SHA256 = "72b005cf392f6b22f969a1218a2f71c8ad4dfb399efd1327ce9d8c1a695ddd57"
MINIMUM_IDENTITY_PERCENT = 80.0
MINIMUM_QUERY_COVERAGE_PERCENT = 90.0
MINIMUM_SUBJECT_COVERAGE_PERCENT = 90.0
BLAST_FIELDS = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "evalue", "bitscore",
    "qstart", "qend", "sstart", "send", "qlen", "slen",
)
ANNOTATION_COLUMNS = [
    *HELPER.RESULT_COLUMNS,
    "coverage_hit_count", "best_identity_percent", "best_query_coverage_percent",
    "best_subject_coverage_percent", "best_bitscore",
]
EVIDENCE_COLUMNS = [
    "accession", "md5", "cds_index", "locus_tag", "qseqid", "sseqid",
    "subject_row_index", "pident_percent", "alignment_length_aa", "mismatch_count",
    "gap_open_count", "evalue", "bitscore",
    "qstart_nt", "qend_nt", "qlen_nt", "query_coverage_percent",
    "sstart_aa", "send_aa", "slen_aa", "subject_coverage_percent",
    "identity_pass", "query_coverage_pass", "subject_coverage_pass", "passes_thresholds",
    "Gene Name", "Product", "Category",
]


def sha256_file(path):
    return HELPER.sha256_file(Path(path))


def atomic_csv(path, fields, rows):
    path = Path(path)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with tmp.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def coverage_percent(start, end, length):
    if length <= 0:
        raise ValueError(f"sequence length must be positive, got {length}")
    return 100.0 * (abs(end - start) + 1) / length


def id_row_index(subject_id):
    match = re.search(r"(?:^|\|)(\d+)$", subject_id)
    if not match:
        raise ValueError(f"Expected numeric PlasAnn subject ID, got {subject_id!r}")
    return int(match.group(1))


def source_hashes(blast_search, database_download, generate_orf_annotation):
    paths = {
        "plasann_blast_search": Path(blast_search.__file__).resolve(),
        "plasann_database_download": Path(database_download.__file__).resolve(),
        "plasann_essential_annotation": Path(inspect.getsourcefile(generate_orf_annotation)).resolve(),
    }
    return {name: sha256_file(path) for name, path in paths.items()}


def blastx_version():
    path = shutil.which("blastx")
    if not path:
        raise RuntimeError("blastx is not on PATH; activate the pinned PlasAnn environment")
    result = subprocess.run([path, "-version"], check=True, capture_output=True, text=True)
    return result.stdout.splitlines()[0].strip()


def config_for(reference_hash, official_sources, blast_version):
    return {
        "source": "PlasAnn CDS BLASTX coverage core",
        "runner_source_sha256": sha256_file(Path(__file__)),
        "step11_helper_sha256": sha256_file(HELPER_PATH),
        "official_source_sha256": official_sources,
        "reference_sha256": reference_hash,
        "plasann_version": HELPER.PLASANN_VERSION,
        "plasann_commit": HELPER.PLASANN_COMMIT,
        "blastx_version": blast_version,
        "blast_evalue": "1e-5",
        "blast_max_target_seqs": 10,
        "minimum_identity_percent": MINIMUM_IDENTITY_PERCENT,
        "minimum_query_coverage_percent": MINIMUM_QUERY_COVERAGE_PERCENT,
        "minimum_subject_coverage_percent": MINIMUM_SUBJECT_COVERAGE_PERCENT,
        "query_coverage_formula": "100*(abs(qend-qstart)+1)/qlen_nt",
        "subject_coverage_formula": "100*(abs(send-sstart)+1)/slen_aa",
        "coverage_is_per_hsp": True,
        "fallback_enabled": False,
    }


def validate_analysis_config(path=CONFIG_PATH):
    config = json.loads(Path(path).read_text(encoding="utf-8"))
    plasann = config.get("plasann", {})
    required = {
        "identity_min_percent": MINIMUM_IDENTITY_PERCENT,
        "query_coverage_min_percent": MINIMUM_QUERY_COVERAGE_PERCENT,
        "subject_coverage_min_percent": MINIMUM_SUBJECT_COVERAGE_PERCENT,
        "fallback_enabled": False,
        "evalue": 1e-5,
        "max_targets": 10,
    }
    if any(plasann.get(key) != value for key, value in required.items()):
        raise ValueError(f"analysis config differs from the frozen PlasAnn 80/90/90 rule: {plasann}")


def load_database_metadata(path):
    required = {"Gene Name", "Product", "Category"}
    pandas_na = {"", "#n/a", "#n/a n/a", "#na", "-1.#ind", "-1.#qnan", "-nan",
                 "1.#ind", "1.#qnan", "<na>", "n/a", "na", "null", "nan"}
    metadata = {}
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"PlasAnn Database.csv lacks required columns: {sorted(required)}")
        for row_index, row in enumerate(reader):
            metadata[row_index] = {
                key: None if (row[key] or "").strip().casefold() in pandas_na else row[key]
                for key in ("Gene Name", "Product", "Category")
            }
    if not metadata:
        raise ValueError(f"PlasAnn Database.csv has no rows: {path}")
    return metadata


def parse_hsps(path, accession, md5, base_rows, sequences, metadata):
    qid_to_index = {f"cds{row['cds_index']}": row["cds_index"] for row in base_rows}
    by_cds, evidence = defaultdict(list), []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(BLAST_FIELDS):
                raise ValueError(f"Malformed blastx row at {path}:{line_number}: expected {len(BLAST_FIELDS)} fields, found {len(fields)}")
            qseqid, sseqid = fields[:2]
            if qseqid not in qid_to_index:
                raise ValueError(f"Unknown BLASTX query ID at {path}:{line_number}: {qseqid}")
            cds_index = qid_to_index[qseqid]
            pident = float(fields[2])
            alignment_length, mismatches, gap_opens = map(int, fields[3:6])
            evalue, bitscore = map(float, fields[6:8])
            qstart, qend, sstart, send, qlen, slen = map(int, fields[8:14])
            if not all(math.isfinite(value) for value in (pident, evalue, bitscore)):
                raise ValueError(f"Nonfinite BLASTX score at {path}:{line_number}")
            if not 0.0 <= pident <= 100.0 or not 0.0 <= evalue <= 1e-5 or bitscore < 0.0:
                raise ValueError(f"BLASTX identity/e-value/score is outside expected bounds at {path}:{line_number}")
            if alignment_length < 1 or mismatches < 0 or gap_opens < 0:
                raise ValueError(f"Invalid BLASTX alignment counts at {path}:{line_number}")
            if qlen != len(sequences[cds_index]):
                raise ValueError(f"BLASTX qlen is not nucleotide CDS length at {path}:{line_number}")
            if not (1 <= qstart <= qlen and 1 <= qend <= qlen and
                    1 <= sstart <= slen and 1 <= send <= slen):
                raise ValueError(f"BLASTX coordinates are outside sequence lengths at {path}:{line_number}")
            row_index = id_row_index(sseqid)
            if row_index not in metadata:
                raise ValueError(f"BLASTX subject row index is absent from Database.csv: {sseqid}")
            qcov = coverage_percent(qstart, qend, qlen)
            scov = coverage_percent(sstart, send, slen)
            pass_identity = pident >= MINIMUM_IDENTITY_PERCENT
            pass_query = qcov >= MINIMUM_QUERY_COVERAGE_PERCENT
            pass_subject = scov >= MINIMUM_SUBJECT_COVERAGE_PERCENT
            passed = pass_identity and pass_query and pass_subject
            base = base_rows[cds_index]
            hit = {
                "accession": accession, "md5": md5, "cds_index": cds_index,
                "locus_tag": base["locus_tag"], "qseqid": qseqid, "sseqid": sseqid,
                "subject_row_index": row_index, "pident_percent": pident,
                "alignment_length_aa": alignment_length, "mismatch_count": mismatches,
                "gap_open_count": gap_opens,
                "evalue": evalue, "bitscore": bitscore,
                "qstart_nt": qstart, "qend_nt": qend, "qlen_nt": qlen,
                "query_coverage_percent": qcov, "sstart_aa": sstart, "send_aa": send,
                "slen_aa": slen, "subject_coverage_percent": scov,
                "identity_pass": str(pass_identity).lower(),
                "query_coverage_pass": str(pass_query).lower(),
                "subject_coverage_pass": str(pass_subject).lower(),
                "passes_thresholds": str(passed).lower(),
                **metadata[row_index],
            }
            by_cds[cds_index].append(hit)
            evidence.append(hit)
    return by_cds, evidence


def classify_cds_rows(base_rows, by_cds):
    rows = []
    for base in base_rows:
        row = dict(base)
        passing = [hit for hit in by_cds.get(base["cds_index"], ()) if hit["passes_thresholds"] == "true"]
        if passing:
            # Counter.most_common preserves first-seen order for tied counts.
            defaults = {"Gene Name": "ORF", "Product": "Open reading frame", "Category": "Uncharacterized"}
            for key in ("Gene Name", "Product", "Category"):
                values = [hit[key] for hit in passing if hit[key] is not None]
                row[key] = Counter(values).most_common(1)[0][0] if values else defaults[key]
            best = max(passing, key=lambda hit: hit["bitscore"])
            row.update(
                coverage_hit_count=len(passing),
                best_identity_percent=best["pident_percent"],
                best_query_coverage_percent=best["query_coverage_percent"],
                best_subject_coverage_percent=best["subject_coverage_percent"],
                best_bitscore=best["bitscore"],
            )
        else:
            row.update({
                "Gene Name": "ORF", "Product": "Open reading frame",
                "Category": "Uncharacterized", "coverage_hit_count": 0,
                "best_identity_percent": "", "best_query_coverage_percent": "",
                "best_subject_coverage_percent": "", "best_bitscore": "",
            })
        # Coverage-only labels leave this PlasAnn-specific translation field empty;
        # the original GenBank translation remains preserved in original_translation.
        row["plasann_translation"] = ""
        rows.append(row)
    return rows


def write_evidence(path, rows):
    path = Path(path)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with tmp.open("wb") as raw:
            with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as compressed:
                import io
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=EVIDENCE_COLUMNS, lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(rows)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def verify_cached_result(success_path, csv_path, evidence_path, expected):
    old = json.loads(Path(success_path).read_text(encoding="utf-8"))
    if old.get("status") != "success" or any(old.get(key) != value for key, value in expected.items()):
        raise RuntimeError(f"Existing completed result has different identity/configuration; choose a new --output: {success_path.parent}")
    if not Path(csv_path).is_file() or not Path(evidence_path).is_file():
        raise RuntimeError(f"Completed result is missing an output file: {success_path.parent}")
    if old.get("cds_annotations_sha256") != sha256_file(csv_path):
        raise RuntimeError(f"Completed CDS annotation checksum mismatch: {csv_path}")
    if old.get("blast_hits_sha256") != sha256_file(evidence_path):
        raise RuntimeError(f"Completed BLAST evidence checksum mismatch: {evidence_path}")
    return old


@contextmanager
def accession_lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def run_blastx(sequences, temp_dir, blast_db, blast_threads):
    query_path = Path(temp_dir) / "cds_queries.fna"
    output_path = Path(temp_dir) / "blastx.tsv"
    with query_path.open("w", encoding="ascii") as handle:
        for index, sequence in enumerate(sequences):
            handle.write(f">cds{index}\n{sequence}\n")
    command = [
        "blastx", "-query", str(query_path), "-db", str(blast_db),
        "-evalue", "1e-5", "-max_target_seqs", "10",
        "-num_threads", str(blast_threads),
        "-outfmt", "6 " + " ".join(BLAST_FIELDS), "-out", str(output_path),
    ]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        message = (exc.stderr or "").strip()[-2000:]
        raise RuntimeError(f"blastx failed with exit status {exc.returncode}: {message}") from exc
    return output_path


def process_one(md5, accession, reference_hash, config, config_hash, output,
                metadata, SeqIO, blast_threads):
    started = time.time()
    out_dir = output / "plasmids" / accession
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "cds_annotations.csv"
    evidence_path = out_dir / "blast_hits.csv.gz"
    success_path = out_dir / "success.json"
    failure_path = out_dir / "failure.json"
    lock_path = out_dir / "run.lock"
    gbff_path = HELPER.GBFF_DIR / f"{accession}.gbff"
    retained_blast_directory = None

    with accession_lock(lock_path):
        protected_success = success_path.exists()

        def log_event(event, **values):
            with (out_dir / "run.log").open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(dict(time=time.time(), event=event, accession=accession, **values), sort_keys=True) + "\n")

        try:
            log_event("begin")
            input_hash = sha256_file(gbff_path)
            expected = {
                "accession": accession, "md5": md5, "input_sha256": input_hash,
                "reference_sha256": reference_hash, "plasann_version": HELPER.PLASANN_VERSION,
                "plasann_commit": HELPER.PLASANN_COMMIT, "config_sha256": config_hash,
                "config": config,
            }
            if protected_success:
                old = verify_cached_result(success_path, csv_path, evidence_path, expected)
                counts = old["counts"]
                failure_path.unlink(missing_ok=True)
                log_event("cached", cds_total=counts["cds_total"])
                return accession, counts, True

            record = SeqIO.read(gbff_path, "genbank")
            observed_md5 = hashlib.md5(str(record.seq).upper().encode("ascii")).hexdigest()
            if observed_md5 != md5:
                raise ValueError(f"GBFF sequence MD5 mismatch: expected {md5}, observed {observed_md5}")
            base_rows, sequences = HELPER.extract_cds(record, accession, md5)
            if sequences:
                retained_blast_directory = tempfile.mkdtemp(prefix="coverage_blastx_", dir=out_dir)
                blast_path = run_blastx(sequences, retained_blast_directory, HELPER.BLAST_PREFIX, blast_threads)
                by_cds, evidence_rows = parse_hsps(blast_path, accession, md5, base_rows, sequences, metadata)
            else:
                by_cds, evidence_rows = {}, []

            rows = classify_cds_rows(base_rows, by_cds)
            atomic_csv(csv_path, ANNOTATION_COLUMNS, rows)
            write_evidence(evidence_path, evidence_rows)
            cds_total, categories, eligible_categories = HELPER.load_csv_counts(csv_path)
            counts = {
                "cds_total": cds_total,
                "protein_eligible": sum(bool(row["protein_eligible"]) for row in base_rows),
                "category_counts": dict(categories),
                "protein_eligible_category_counts": dict(eligible_categories),
                "coverage_hit_cds": sum(bool(by_cds.get(i)) and any(h["passes_thresholds"] == "true" for h in by_cds[i]) for i in range(len(base_rows))),
                "passing_hsp_total": sum(h["passes_thresholds"] == "true" for h in evidence_rows),
            }
            manifest = dict(
                expected, status="success", counts=counts,
                runner_source_sha256=sha256_file(Path(__file__)),
                cds_annotations_sha256=sha256_file(csv_path),
                blast_hits_sha256=sha256_file(evidence_path),
                retained_blast_directory=retained_blast_directory,
                timing_seconds=round(time.time() - started, 3),
            )
            HELPER.write_json(success_path, manifest)
            failure_path.unlink(missing_ok=True)
            log_event("success", cds_total=cds_total, passing_hsp_total=counts["passing_hsp_total"], seconds=manifest["timing_seconds"])
            print(f"{accession}: {cds_total} CDS, {counts['coverage_hit_cds']} with passing HSPs", flush=True)
            return accession, counts, False
        except Exception as exc:
            if not success_path.exists():
                failure = {
                    "status": "failed", "accession": accession, "md5": md5,
                    "error": f"{type(exc).__name__}: {exc}",
                    "retained_blast_directory": retained_blast_directory,
                    "timing_seconds": round(time.time() - started, 3),
                }
                HELPER.write_json(failure_path, failure)
                log_event("error", error=failure["error"])
            print(f"FAILED {accession}: {type(exc).__name__}: {exc}", flush=True)
            raise


def write_counts(output, completed):
    categories, eligible = Counter(), Counter()
    for counts in completed.values():
        categories.update(counts["category_counts"])
        eligible.update(counts["protein_eligible_category_counts"])
    atomic_csv(output / "category_counts.csv", ["Category", "cds_count", "protein_eligible_count"],
               ({"Category": category, "cds_count": categories[category],
                 "protein_eligible_count": eligible[category]} for category in sorted(categories)))
    return categories


def validate_cpu_request(workers, blast_threads, environ=None):
    environ = os.environ if environ is None else environ
    requested = workers * blast_threads
    allocated = environ.get("SLURM_CPUS_PER_TASK")
    if allocated is not None:
        try:
            allocated = int(allocated)
        except ValueError as exc:
            raise ValueError(f"Invalid SLURM_CPUS_PER_TASK={allocated!r}") from exc
        if requested > allocated:
            raise ValueError(f"workers x blast_threads ({workers} x {blast_threads} = {requested}) exceeds SLURM_CPUS_PER_TASK={allocated}")


def self_test():
    global run_blastx, parse_hsps
    exact = coverage_percent(1, 270, 300)
    reverse_subject = coverage_percent(90, 1, 100)
    assert exact == 90.0 and reverse_subject == 90.0
    base = [{"accession": "a", "md5": "m", "cds_index": 0, "locus_tag": "L", "original_product": "x",
             "original_translation": "KEEP", "protein_eligible": True}]
    meta = {0: {"Gene Name": "A", "Product": "p", "Category": "C"}}
    with tempfile.TemporaryDirectory() as tmp:
        blast = Path(tmp) / "hits.tsv"
        # Exact 80/90/90 boundary, with reverse subject coordinates and nt/aa lengths.
        blast.write_text("cds0\t0\t80\t90\t18\t0\t1e-5\t50\t270\t1\t90\t1\t300\t100\n", encoding="utf-8")
        grouped, hits = parse_hsps(blast, "a", "m", base, ["A" * 300], meta)
        assert hits[0]["passes_thresholds"] == "true" and hits[0]["query_coverage_percent"] == 90.0
        assert hits[0]["subject_coverage_percent"] == 90.0
        assert classify_cds_rows(base, grouped)[0]["Category"] == "C"
        evidence_path = Path(tmp) / "blast_hits.csv.gz"
        write_evidence(evidence_path, hits)
        with gzip.open(evidence_path, "rt", newline="", encoding="utf-8") as handle:
            saved_evidence = list(csv.DictReader(handle))
        assert len(saved_evidence) == 1 and saved_evidence[0]["passes_thresholds"] == "true"

        # Each threshold rejects independently; sub-threshold HSPs never trigger a fallback.
        for pident, qend, send in ((79.99, 270, 1), (80, 269, 1), (80, 270, 10)):
            blast.write_text(f"cds0\t0\t{pident}\t90\t18\t0\t1e-10\t50\t1\t{qend}\t1\t{send}\t300\t100\n", encoding="utf-8")
            grouped, hits = parse_hsps(blast, "a", "m", base, ["A" * 300], meta)
            assert hits[0]["passes_thresholds"] == "false"
            result = classify_cds_rows(base, grouped)[0]
            assert (result["Gene Name"], result["Product"], result["Category"]) == ("ORF", "Open reading frame", "Uncharacterized")

        blast.write_text("cds0\t0\t80\t90\t18\t0\t1.1e-5\t50\t1\t270\t1\t90\t300\t100\n", encoding="utf-8")
        try:
            parse_hsps(blast, "a", "m", base, ["A" * 300], meta)
        except ValueError:
            pass
        else:
            raise AssertionError("an HSP above the fixed e-value was accepted")

        missing_meta = {0: {"Gene Name": None, "Product": None, "Category": "C"}}
        blast.write_text("cds0\t0\t80\t90\t18\t0\t1e-5\t50\t1\t270\t1\t90\t300\t100\n", encoding="utf-8")
        grouped, _ = parse_hsps(blast, "a", "m", base, ["A" * 300], missing_meta)
        missing_row = classify_cds_rows(base, grouped)[0]
        assert (missing_row["Gene Name"], missing_row["Product"], missing_row["Category"]) == ("ORF", "Open reading frame", "C")

        # Per-field majority includes all passing HSPs and keeps the first-seen value on a tie.
        tied = [
            {"passes_thresholds": "true", "Gene Name": "first", "Product": "p1", "Category": "c1", "bitscore": 5,
             "pident_percent": 80, "query_coverage_percent": 90, "subject_coverage_percent": 90},
            {"passes_thresholds": "true", "Gene Name": "second", "Product": "p2", "Category": "c2", "bitscore": 6,
             "pident_percent": 80, "query_coverage_percent": 90, "subject_coverage_percent": 90},
        ]
        tied_row = classify_cds_rows(base, {0: tied})[0]
        assert (tied_row["Gene Name"], tied_row["Product"], tied_row["Category"]) == ("first", "p1", "c1")

        # Cache identity changes invalidate reuse; output corruption is rejected.
        csv_path, evidence_path = Path(tmp) / "a.csv", Path(tmp) / "a.gz"
        csv_path.write_text("csv", encoding="utf-8")
        evidence_path.write_text("hits", encoding="utf-8")
        expected = {"accession": "a", "config_sha256": "one"}
        success = Path(tmp) / "success.json"
        manifest = dict(expected, status="success", cds_annotations_sha256=sha256_file(csv_path),
                        blast_hits_sha256=sha256_file(evidence_path))
        success.write_text(json.dumps(manifest), encoding="utf-8")
        assert verify_cached_result(success, csv_path, evidence_path, expected) == manifest
        try:
            verify_cached_result(success, csv_path, evidence_path, {"accession": "a", "config_sha256": "two"})
        except RuntimeError:
            pass
        else:
            raise AssertionError("cache configuration change was accepted")
        csv_path.write_text("corrupt", encoding="utf-8")
        try:
            verify_cached_result(success, csv_path, evidence_path, expected)
        except RuntimeError:
            pass
        else:
            raise AssertionError("corrupt cached output was accepted")

        # A successful worker keeps its query and raw BLAST output for inspection.
        gbff_dir = Path(tmp) / "gbff"
        gbff_dir.mkdir()
        (gbff_dir / "a.gbff").write_text("fixture", encoding="utf-8")
        worker_md5 = hashlib.md5(b"A").hexdigest()
        worker_rows = [{"accession": "a", "md5": worker_md5, "cds_index": 0,
                        "locus_tag": "L", "protein_eligible": True}]

        class FakeSeqIO:
            @staticmethod
            def read(_path, _format):
                return type("Record", (), {"seq": "A"})()

        def fake_run_blastx(_sequences, temp_dir, _blast_db, _blast_threads):
            (Path(temp_dir) / "cds_queries.fna").write_text(">cds0\nA\n", encoding="ascii")
            blast_path = Path(temp_dir) / "blastx.tsv"
            blast_path.write_text("raw BLAST output\n", encoding="ascii")
            return blast_path

        def fail_temporary_directory(*_args, **_kwargs):
            raise OSError(39, "Directory not empty")

        originals = (HELPER.GBFF_DIR, HELPER.extract_cds, HELPER.load_csv_counts,
                     run_blastx, parse_hsps, tempfile.TemporaryDirectory)
        try:
            HELPER.GBFF_DIR = gbff_dir
            HELPER.extract_cds = lambda *_args: (worker_rows, ["A"])
            HELPER.load_csv_counts = lambda _path: (1, {}, {})
            run_blastx = fake_run_blastx
            parse_hsps = lambda *_args: ({}, [])
            tempfile.TemporaryDirectory = fail_temporary_directory
            process_one(worker_md5, "a", "reference", {"self_test": True}, "config",
                        Path(tmp) / "worker_output", {}, FakeSeqIO, 1)
        finally:
            (HELPER.GBFF_DIR, HELPER.extract_cds, HELPER.load_csv_counts,
             run_blastx, parse_hsps, tempfile.TemporaryDirectory) = originals
        worker_dir = Path(tmp) / "worker_output" / "plasmids" / "a"
        manifest = json.loads((worker_dir / "success.json").read_text(encoding="utf-8"))
        retained = Path(manifest["retained_blast_directory"])
        assert (retained / "cds_queries.fna").is_file() and (retained / "blastx.tsv").is_file()

    sample = list(range(25))
    shards = [HELPER.select_shard(sample, i, 8) for i in range(8)]
    assert sorted(item for shard in shards for item in shard) == sample
    assert all(set(shards[i]).isdisjoint(shards[j]) for i in range(8) for j in range(i))
    assert id_row_index("gnl|BL_ORD_ID|17") == 17
    try:
        validate_cpu_request(10, 4, {"SLURM_CPUS_PER_TASK": "39"})
    except ValueError:
        pass
    else:
        raise AssertionError("CPU over-allocation was accepted")
    validate_cpu_request(10, 4, {"SLURM_CPUS_PER_TASK": "40"})
    print("Self-test passed: coverage boundaries/units, reverse coordinates, consensus, fallback, cache integrity, IDs, sharding, and CPU allocation.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--blast-threads", type=int, default=4)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--limit", type=int, help="Process the first N selected representatives; mark the summary partial")
    parser.add_argument("--accession", action="append", help="Process this representative accession; repeatable and partial")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.workers < 1 or args.blast_threads < 1 or (args.limit is not None and args.limit < 1):
        parser.error("--workers, --blast-threads, and --limit must be positive")
    if args.shard_count < 1 or not 0 <= args.shard_index < args.shard_count:
        parser.error("Require --shard-count >= 1 and 0 <= --shard-index < --shard-count")
    try:
        validate_cpu_request(args.workers, args.blast_threads)
    except ValueError as exc:
        parser.error(str(exc))

    validate_analysis_config()
    reference_hash = HELPER.require_database()
    if reference_hash != EXPECTED_REFERENCE_SHA256:
        raise RuntimeError(f"Frozen PlasAnn reference SHA mismatch: expected {EXPECTED_REFERENCE_SHA256}, got {reference_hash}")
    pd, SeqIO, blast_search, database_download, generate_orf_annotation = HELPER.software_modules()
    del pd
    version = blastx_version()
    config = config_for(reference_hash, source_hashes(blast_search, database_download, generate_orf_annotation), version)
    config_hash = HELPER.json_sha256(config)
    metadata = load_database_metadata(HELPER.DATABASE_CSV)

    full_mapping = HELPER.read_mapping()
    shard_mapping = HELPER.select_shard(full_mapping, args.shard_index, args.shard_count)
    mapping = shard_mapping
    if args.accession:
        requested = set(args.accession)
        known = {accession for _, accession in full_mapping}
        unknown = requested - known
        if unknown:
            parser.error("Unknown representative accession(s): " + ", ".join(sorted(unknown)))
        mapping = [item for item in mapping if item[1] in requested]
    if args.limit is not None:
        mapping = mapping[:args.limit]

    args.output.mkdir(parents=True, exist_ok=True)
    completed, cached_count, failures = {}, 0, []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(process_one, md5, accession, reference_hash, config, config_hash,
                        args.output, metadata, SeqIO, args.blast_threads): accession
            for md5, accession in mapping
        }
        for future in concurrent.futures.as_completed(futures):
            accession = futures[future]
            try:
                _, counts, was_cached = future.result()
                completed[accession] = counts
                cached_count += int(was_cached)
            except Exception as exc:
                failures.append({"accession": accession, "error": f"{type(exc).__name__}: {exc}"})

    summary_output = args.output / "shards" / f"{args.shard_index}_of_{args.shard_count}"
    summary_output.mkdir(parents=True, exist_ok=True)
    categories = write_counts(summary_output, completed)
    partial_scope = args.limit is not None or bool(args.accession)
    summary = {
        "runner_source_sha256": sha256_file(Path(__file__)),
        "status": "failed" if failures else "PARTIAL" if partial_scope else "success",
        "reference_sha256": reference_hash, "plasann_version": HELPER.PLASANN_VERSION,
        "plasann_commit": HELPER.PLASANN_COMMIT, "config": config,
        "config_sha256": config_hash, "workers": args.workers,
        "blast_threads": args.blast_threads, "shard_index": args.shard_index,
        "shard_count": args.shard_count, "plasmids_assigned_to_shard": len(shard_mapping),
        "plasmids_requested": len(mapping), "plasmids_completed": len(completed),
        "plasmids_cached": cached_count, "failures": failures,
        "cds_total": sum(c["cds_total"] for c in completed.values()),
        "protein_eligible_total": sum(c["protein_eligible"] for c in completed.values()),
        "coverage_hit_cds_total": sum(c["coverage_hit_cds"] for c in completed.values()),
        "passing_hsp_total": sum(c["passing_hsp_total"] for c in completed.values()),
        "category_counts": dict(categories),
    }
    HELPER.write_json(summary_output / "run_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
