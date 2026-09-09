#!/usr/bin/env python3
"""
13_strict_analysis.py
Strict-mode (adjacent-only pairing) sensitivity pipeline: protein dedup ->
DIAMOND vs 4 DBs -> enrichment, all on the strict cargo/background sets
(5strict_*.faa from historical step 5). Mirrors historical steps 6-7.

The enrichment table reports BOTH countings side by side (main analysis is
weighted = CDS-level; unweighted is shown only as a sensitivity):
  - weighted   = CDS-level (each cargo/background CDS counts once;
                 identical proteins on different plasmids count each time)
  - unweighted = unique-protein level (identical protein sequences count
                 once regardless of how many plasmids carry them)

Usage:  python scripts/13_strict_analysis.py [threads] [min_pident] [min_qcov]
        defaults: 8 90 60

Input:  2_Result/5strict_Cargo.faa, 2_Result/5strict_Background.faa
        1_Data/database/{card,vfdb,bacmet,tadb}.dmnd
Output: 2_Result/8_Cargo_unique.faa, 8_Cargo_unique_map.csv,
        2_Result/8_Background_unique.faa, 8_Background_unique_map.csv,
        2_Result/8_{Cargo,Background}_vs_{card,vfdb,bacmet,tadb}.tsv,
        2_Result/8_Enrichment_table.csv, 2_Result/8_Summary.txt
"""
import csv
import hashlib
import os
import subprocess
import sys

from scipy.stats import fisher_exact

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
DATA_DIR = os.path.abspath(os.environ.get(
    "NEE_DATA_DIR", os.path.join(ROOT, "1_Data")))
RES = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))
DBDIR = os.path.abspath(os.environ.get(
    "NEE_DATABASE_DIR", os.path.join(DATA_DIR, "database")))
DIAMOND = os.environ.get("DIAMOND_BIN", "diamond")

THREADS = int(sys.argv[1]) if len(sys.argv) > 1 else 8
MIN_PIDENT = float(sys.argv[2]) if len(sys.argv) > 2 else 90.0
MIN_QCOV = float(sys.argv[3]) if len(sys.argv) > 3 else 60.0

CATS = ["card", "vfdb", "bacmet", "tadb"]
LABELS = {"card": "Antibiotic resistance (CARD)",
          "vfdb": "Virulence (VFDB)",
          "bacmet": "Metal/biocide resistance (BacMet)",
          "tadb": "Toxin-antitoxin (TADB)"}


# ---------- 1. protein dedup (same logic as 08_dedup_proteins.py) ----------

def read_fasta(path):
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


def dedup(src, out_faa, out_map):
    by_hash = {}
    n_total = 0
    for header, seq in read_fasta(src):
        n_total += 1
        h = hashlib.md5(seq.encode()).hexdigest()
        if h in by_hash:
            by_hash[h][0] += 1
        else:
            by_hash[h] = [1, header]
    seen = set()
    with open(out_faa, "w") as fo:
        for header, seq in read_fasta(src):
            h = hashlib.md5(seq.encode()).hexdigest()
            if h in seen:
                continue
            seen.add(h)
            fo.write(f">{h}\n{seq}\n")
    with open(out_map, "w") as fm:
        fm.write("md5,n_copies,representative_header\n")
        for h, (n, hdr) in sorted(by_hash.items()):
            fm.write(f"{h},{n},{hdr}\n")
    print(f"  {os.path.basename(src)}: {n_total} -> {len(by_hash)} unique",
          flush=True)
    return by_hash


# ---------- 2. diamond ----------

def run_diamond(setname):
    for cat in CATS:
        o = os.path.join(RES, f"8_{setname}_vs_{cat}.tsv")
        if os.path.exists(o) and os.path.getsize(o) > 0:
            print(f"  skip (exists): {o}", flush=True)
            continue
        print(f"  blastp: {setname} vs {cat}", flush=True)
        subprocess.run(
            [DIAMOND, "blastp", "--very-sensitive", "--threads",
             str(THREADS),
             "-q", os.path.join(RES, f"8_{setname}_unique.faa"),
             "-d", os.path.join(DBDIR, cat),
             "-o", o,
             "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
             "qcovhsp", "evalue", "bitscore", "stitle",
             "-e", "1e-10", "-k", "1"],
            check=True)


# ---------- 3. enrichment (weighted + unweighted) ----------

def load_counts(setname, cat, umap):
    """Return (weighted, unweighted) hit counts passing thresholds."""
    w = u = 0
    with open(os.path.join(RES, f"8_{setname}_vs_{cat}.tsv"),
              errors="replace") as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 7:
                continue
            if float(row[2]) >= MIN_PIDENT and float(row[4]) >= MIN_QCOV:
                w += umap.get(row[0], [1, ""])[0]
                u += 1
    return w, u


def bh_fdr(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        val = min(prev, pvals[i] * m / (m - rank + 1))
        q[i] = val
        prev = val
    return q


def main():
    print("[1/3] protein dedup (strict sets)", flush=True)
    maps = {}
    for setname, src in (("Cargo", "5strict_Cargo.faa"),
                         ("Background", "5strict_Background.faa")):
        maps[setname] = dedup(
            os.path.join(RES, src),
            os.path.join(RES, f"8_{setname}_unique.faa"),
            os.path.join(RES, f"8_{setname}_unique_map.csv"))

    print("[2/3] diamond", flush=True)
    for s in ("Cargo", "Background"):
        run_diamond(s)

    print("[3/3] enrichment (weighted = CDS-level only)", flush=True)
    tot_w = {s: sum(v[0] for v in maps[s].values()) for s in maps}

    rows = []
    for cat in CATS:
        cw, _ = load_counts("Cargo", cat, maps["Cargo"])
        bw, _ = load_counts("Background", cat, maps["Background"])
        odds, p = fisher_exact(
            [[cw, tot_w["Cargo"] - cw], [bw, tot_w["Background"] - bw]])
        rows.append([cat, LABELS[cat],
                     cw, tot_w["Cargo"], f"{cw / tot_w['Cargo']:.4f}",
                     bw, tot_w["Background"], f"{bw / tot_w['Background']:.4f}",
                     f"{odds:.3f}", f"{p:.3e}"])

    q_w = bh_fdr([float(r[9]) for r in rows])
    for r, q in zip(rows, q_w):
        r.append(f"{q:.3e}")

    out_csv = os.path.join(RES, "8_Enrichment_table.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["category", "label", "cargo_hits", "cargo_total",
                    "cargo_frac", "bg_hits", "bg_total", "bg_frac",
                    "odds_ratio", "p_fisher", "q_fdr_bh"])
        w.writerows(rows)

    lines = [
        f"# Step 8: STRICT-mode enrichment "
        f"(pident>={MIN_PIDENT:g}, qcov>={MIN_QCOV:g}, weighted = CDS-level)",
        f"cargo_total: {tot_w['Cargo']}, background_total: "
        f"{tot_w['Background']} (copy-weighted)",
    ]
    for r in rows:
        lines.append(f"{r[1]}: cargo {r[4]} vs bg {r[7]} "
                     f"| OR={r[8]} | P={r[9]} | q={r[10]}")
    summary = "\n".join(lines)
    print(summary)
    with open(os.path.join(RES, "8_Summary.txt"), "w") as f:
        f.write(summary + "\n")

    # figure (same style as step 7)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = [LABELS[c].replace(" (", "\n(") for c in CATS]
    cargo_frac = [float(r[4]) for r in rows]
    bg_frac = [float(r[7]) for r in rows]
    x = range(len(labels))
    w = 0.35
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.bar([i - w / 2 for i in x], cargo_frac, w,
           label="Transposon cargo (strict)", color="#C44E52")
    ax.bar([i + w / 2 for i in x], bg_frac, w,
           label="Background (same plasmids)", color="#8C8C8C")
    ymax = max(max(cargo_frac), max(bg_frac))
    for i, r in enumerate(rows):
        ax.text(i, max(cargo_frac[i], bg_frac[i]) + ymax * 0.02,
                f"OR={r[8]}\nq={r[10]}", ha="center", fontsize=7.5)
    ax.set_ylim(0, ymax * 1.25)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Fraction of proteins with database hit")
    ax.set_title("STRICT cargo vs background "
                 f"(pident>={MIN_PIDENT:g}%, qcov>={MIN_QCOV:g}%)",
                 fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "8_Enrichment_figure.pdf"))
    fig.savefig(os.path.join(RES, "8_Enrichment_figure.png"), dpi=150)
    print(f"figure -> {os.path.join(RES, '8_Enrichment_figure.pdf')} (+ .png)")


if __name__ == "__main__":
    main()
