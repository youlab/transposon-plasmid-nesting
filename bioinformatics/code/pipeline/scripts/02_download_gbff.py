#!/usr/bin/env python3
"""
02_download_gbff.py
Batch-download GenBank files (with sequence + CDS translations) for the
flagged plasmid accessions from NCBI, via Biopython Entrez efetch.

One .gbff per accession in 1_Data/gbff/. Resumable: skips accessions whose
file already exists. Uses chunks with retries and a polite delay.

NCBI requires an email address for E-utilities (used only so NCBI can
contact you about abnormal usage; it is not verified). An API key is
optional but raises the rate limit from 3 to 10 req/s.

Credentials are read only from the NCBI_EMAIL and NCBI_API_KEY environment
variables. Never store a real API key in the repository.

Usage:
    python scripts/02_download_gbff.py [accession_list.txt]
    # or override in the shell:  export NCBI_EMAIL=... NCBI_API_KEY=...

Defaults (relative to the repository root):
    Input:   2_Result/1_Hit_plasmids.txt
    Output:  1_Data/gbff/<accession>.gbff
    Summary: 2_Result/2_Download_assessions_summary.txt
    Failed:  2_Result/2_Failed_assessions.txt
             (re-run this script with that file as argument to retry)

Notes:
    - rettype=gbwithparts gives full sequence; step 3 derives .fna for
      ISEScan from the same file, so only ONE download per accession.
    - ~26k accessions -> a few hours with an API key. Run under nohup/tmux.
"""
import datetime
import os
import sys
import time

from Bio import Entrez

CHUNK = 100            # accessions per efetch request
SLEEP_NO_KEY = 0.4     # ~3 req/s without API key
SLEEP_WITH_KEY = 0.12  # ~10 req/s with API key
RETRIES = 5

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
DATA_DIR = os.path.abspath(os.environ.get(
    "NEE_DATA_DIR", os.path.join(ROOT, "1_Data")))
RESULTS_DIR = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))
OUTDIR = os.path.join(DATA_DIR, "gbff")
SUMMARY_TXT = os.path.join(RESULTS_DIR, "2_Download_assessions_summary.txt")
FAILED_TXT = os.path.join(RESULTS_DIR, "2_Failed_assessions.txt")

Entrez.email = os.environ.get("NCBI_EMAIL", "")
Entrez.api_key = os.environ.get("NCBI_API_KEY") or None
SLEEP = SLEEP_WITH_KEY if Entrez.api_key else SLEEP_NO_KEY


def fetch_chunk(accs):
    ids = ",".join(accs)
    for attempt in range(RETRIES):
        try:
            h = Entrez.efetch(db="nuccore", id=ids,
                              rettype="gbwithparts", retmode="text")
            text = h.read()
            h.close()
            return text
        except Exception as e:
            wait = 5 * (attempt + 1)
            print(f"  retry {attempt + 1}/{RETRIES} after error: {e}; "
                  f"sleeping {wait}s", file=sys.stderr, flush=True)
            time.sleep(wait)
    return None


def split_gbff(text):
    """Split a multi-record GenBank text into per-record strings."""
    records, buf = [], []
    for line in text.splitlines(keepends=True):
        buf.append(line)
        if line.startswith("//"):
            records.append("".join(buf))
            buf = []
    return records


def accession_of(record_text):
    # VERSION line holds the accession.version, e.g. "VERSION     NZ_CP107404.1"
    for line in record_text.splitlines():
        if line.startswith("VERSION"):
            return line.split()[1].strip()
    return None


def main():
    if not Entrez.email:
        sys.exit("ERROR: set NCBI_EMAIL first "
                 "(NCBI requires an email for E-utilities), e.g.\n"
                 "  export NCBI_EMAIL=you@duke.edu")

    acc_file = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        RESULTS_DIR, "1_Hit_plasmids.txt")
    os.makedirs(OUTDIR, exist_ok=True)
    accs = [l.strip() for l in open(acc_file) if l.strip()]
    todo = [a for a in accs
            if not os.path.exists(os.path.join(OUTDIR, a + ".gbff"))]
    n_skipped = len(accs) - len(todo)
    print(f"{len(accs)} accessions total, {n_skipped} already present, "
          f"{len(todo)} to download", flush=True)

    failed = []
    for i in range(0, len(todo), CHUNK):
        chunk = todo[i:i + CHUNK]
        text = fetch_chunk(chunk)
        if text is None:
            failed.extend(chunk)
            print(f"FAILED chunk starting at {i}: {len(chunk)} accessions",
                  file=sys.stderr, flush=True)
            continue
        got = set()
        for rec in split_gbff(text):
            acc = accession_of(rec)
            if acc:
                with open(os.path.join(OUTDIR, acc + ".gbff"), "w") as f:
                    f.write(rec)
                got.add(acc)
        missing = [a for a in chunk if a not in got]
        failed.extend(missing)
        if missing:
            print(f"  {len(missing)} accessions not returned in this chunk: "
                  f"{missing[:5]}", file=sys.stderr, flush=True)
        done = min(i + CHUNK, len(todo))
        print(f"{done}/{len(todo)}", flush=True)
        time.sleep(SLEEP)

    n_gbff = sum(1 for f in os.listdir(OUTDIR) if f.endswith(".gbff"))
    if failed:
        with open(FAILED_TXT, "w") as f:
            f.write("\n".join(failed) + "\n")
        print(f"{len(failed)} failed -> {FAILED_TXT} "
              f"(re-run this script with that file to retry)", flush=True)
    elif os.path.exists(FAILED_TXT):
        os.remove(FAILED_TXT)  # stale list from a previous run; none failed now

    with open(SUMMARY_TXT, "w") as f:
        f.write("\n".join([
            "# Step 2: gbff download summary",
            f"date: {datetime.datetime.now().isoformat(timespec='seconds')}",
            f"input_file: {acc_file}",
            f"input_accessions: {len(accs)}",
            f"already_present_skipped: {n_skipped}",
            f"attempted_this_run: {len(todo)}",
            f"downloaded_this_run: {len(todo) - len(failed)}",
            f"failed_this_run: {len(failed)}",
            f"gbff_files_total: {n_gbff}",
            f"failed_list: {FAILED_TXT if failed else '(none)'}",
            (f"retry: python scripts/02_download_gbff.py {FAILED_TXT}"
             if failed else "retry: (nothing to retry)"),
            "",
        ]))
    print(f"summary -> {SUMMARY_TXT}", flush=True)


if __name__ == "__main__":
    main()
