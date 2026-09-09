#!/usr/bin/env python3
"""
15_is_sequence_jumps.py
Sequence-level analysis of identical IS copies at distinct insertion contexts.

Analysis logic:
  1. Candidate filter: plasmids with >=2 copies of the same IS family
     (family label is ONLY a filter; no identity claim).
  2. Element identity: two copies are "the same transposon" iff their DNA
     sequences are EXACTLY identical, up to orientation (md5 of the copy
     sequence or its reverse complement). No alignment or identity threshold
     is used, and no filter is applied to the ISEScan type field.
  3. Insertion-site comparison: for each identical pair, the 10 bp OUTSIDE
     each end of the element (in the element's own orientation, so inverted
     copies are handled: left-left + right-right after rev-comp where
     needed). Same site = both 10 bp flanks identical (this also excludes
     tandem amplification arrays, whose junctions are identical); different
     site = anything else. Caveat: post-insertion SNP drift in the 10 bp
     counts as "different" (anti-conservative direction, noted in output).
  4. A plasmid is a "multiple-jump plasmid" if any identical-sequence group
     has >=2 distinct contexts. Funnel counts reported at each stage.
  5. Representative example records include accession, coordinates,
     orientation, both 10 bp flanks per copy).

Usage:  python scripts/15_is_sequence_jumps.py
Input:  2_Result/4_ISEScan/fna/*.fna.tsv, 1_Data/fna/*.fna,
        2_Result/3_Seqhash_map.csv
Output: 2_Result/10_Summary.txt, 2_Result/10_Per_plasmid.csv,
        2_Result/10_Jump_examples.txt
"""
import csv
import hashlib
import os
from collections import Counter, defaultdict

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
DATA_DIR = os.path.abspath(os.environ.get(
    "NEE_DATA_DIR", os.path.join(ROOT, "1_Data")))
RES = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))
FNA_DIR = os.path.join(DATA_DIR, "fna")
TSV_DIR = os.path.join(RES, "4_ISEScan", "fna")
MAP_CSV = os.path.join(RES, "3_Seqhash_map.csv")
OUT_SUM = os.path.join(RES, "10_Summary.txt")
OUT_PER = os.path.join(RES, "10_Per_plasmid.csv")
OUT_EX = os.path.join(RES, "10_Jump_examples.txt")

FLANK = 10
MAX_EXAMPLES = 10


def parse_tsv(path):
    """copies: list of (family, begin, end) 1-based, sorted by begin."""
    copies = []
    with open(path) as f:
        for line in f:
            if line.startswith("seqID") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 22:
                continue
            copies.append((p[1], int(p[3]), int(p[4])))
    copies.sort(key=lambda c: c[1])
    return copies


def read_fna(path):
    with open(path) as f:
        return "".join(l.strip() for l in f if not l.startswith(">")).upper()


def circular_slice(seq, start, length):
    L = len(seq)
    return "".join(seq[(start + i) % L] for i in range(length))


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canon_key(s):
    """Orientation-free identity key for an element sequence."""
    return min(s, revcomp(s))


def element_context(seq, begin, end, forward):
    """10 bp outside each end, IN THE ELEMENT'S OWN ORIENTATION:
    returns (left_flank, right_flank) as strings 5'->3' in element frame."""
    up = circular_slice(seq, begin - 1 - FLANK, FLANK)   # upstream of begin
    dn = circular_slice(seq, end, FLANK)                 # downstream of end
    if forward:
        return up, dn
    return revcomp(dn), revcomp(up)


def main():
    md5rep = {}
    with open(MAP_CSV) as f:
        for row in csv.DictReader(f):
            md5rep[row["md5"]] = row["representative_acc"]
    n_total = len(md5rep)

    import glob
    n_stage1 = n_stage2 = n_stage3 = 0
    fam_stage3 = Counter()
    per_rows, examples = [], []
    tsv_files = sorted(glob.glob(os.path.join(TSV_DIR, "*.fna.tsv")))
    for k, tsv in enumerate(tsv_files):
        md5 = os.path.basename(tsv)[:-len(".fna.tsv")]
        copies = parse_tsv(tsv)
        fams = defaultdict(list)
        for fam, b, e in copies:
            fams[fam].append((b, e))
        fams = {f: c for f, c in fams.items() if len(c) >= 2}
        if not fams:
            continue
        n_stage1 += 1
        fna = os.path.join(FNA_DIR, md5 + ".fna")
        if not os.path.exists(fna):
            continue
        seq = read_fna(fna)

        plasmid_has_identical = False
        plasmid_jump_fams = []
        plasmid_maxctx = 0
        for fam, coords in sorted(fams.items()):
            # group copies by exact sequence identity (orientation-free)
            groups = defaultdict(list)
            for b, e in coords:
                s = circular_slice(seq, b - 1, e - b + 1)
                groups[canon_key(s)].append((b, e, s))
            for key, mem in groups.items():
                if len(mem) < 2:
                    continue
                plasmid_has_identical = True
                # contexts in element orientation, keyed per copy
                ctxs = []
                for b, e, s in mem:
                    forward = (s == key)
                    lf, rf = element_context(seq, b, e, forward)
                    ctxs.append((b, e, forward, lf, rf))
                # distinct contexts: exact (left,right) pair strings
                distinct = {}
                for b, e, fw, lf, rf in ctxs:
                    distinct.setdefault((lf, rf), []).append((b, e, fw))
                nctx = len(distinct)
                plasmid_maxctx = max(plasmid_maxctx, nctx)
                if nctx >= 2:
                    plasmid_jump_fams.append(fam)
                    if len(examples) < MAX_EXAMPLES and nctx >= 2:
                        examples.append(
                            (md5rep.get(md5, "NA"), md5, fam, len(key),
                             len(mem), nctx, ctxs, distinct))
        if plasmid_has_identical:
            n_stage2 += 1
        if plasmid_jump_fams:
            n_stage3 += 1
            for f in set(plasmid_jump_fams):
                fam_stage3[f] += 1
        per_rows.append([md5rep.get(md5, "NA"), md5,
                         len(fams), int(plasmid_has_identical),
                         len(set(plasmid_jump_fams)), plasmid_maxctx,
                         ";".join(sorted(set(plasmid_jump_fams)))])
        if (k + 1) % 4000 == 0:
            print(f"{k + 1}/{len(tsv_files)} plasmids scanned", flush=True)

    with open(OUT_PER, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["accession", "md5", "n_families_ge2", "has_identical_pair",
                    "n_jump_families", "max_distinct_contexts", "jump_families"])
        w.writerows(per_rows)

    lines = [
        "# Step 10: sequence-level evidence of multiple jumps",
        f"rule: identical IS sequence (md5, orientation-free) at >=2 distinct",
        f"insertion contexts (outer {FLANK} bp flanks differ)",
        "",
        f"stage 1 - plasmids with >=2 same-family IS copies: {n_stage1}",
        f"stage 2 - ...with >=1 EXACTLY identical copy pair: {n_stage2} "
        f"({n_stage2 / max(n_stage1, 1):.1%} of stage 1)",
        f"stage 3 - ...with the identical element at >=2 DISTINCT insertion "
        f"contexts: {n_stage3} ({n_stage3 / max(n_stage1, 1):.1%} of stage 1, "
        f"{n_stage3 / max(n_total, 1):.1%} of all {n_total:,} plasmids)",
        "",
        "top families among stage-3 plasmids:",
    ]
    lines += [f"  {fam}: {n}" for fam, n in fam_stage3.most_common(15)]
    lines += ["",
              "caveats: (i) post-insertion SNP drift in the 10 bp",
              "flank counts as 'different site' (anti-conservative); (ii) bare",
              "tandem IS-IS pairs count as different sites (adjacent insertion",
              "is still an insertion event; IS26-style targeted insertion);",
              "tandem amplification ARRAYS are excluded (junctions identical)."]

    with open(OUT_SUM, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # examples
    with open(OUT_EX, "w") as f:
        for acc, md5, fam, elen, ncopy, nctx, ctxs, distinct in examples:
            f.write(f"=== {acc} (md5 {md5}) | family {fam} | element length "
                    f"{elen} bp | {ncopy} identical copies | {nctx} distinct "
                    f"insertion contexts\n")
            for i, ((lf, rf), mem) in enumerate(sorted(distinct.items())):
                b, e, fw = mem[0]
                f.write(f"  context {i + 1}: copy at {b}..{e} "
                        f"({'forward' if fw else 'inverted'}"
                        f"{', shared by ' + str(len(mem)) + ' copies' if len(mem) > 1 else ''})\n")
                f.write(f"    left flank : 5'-{lf}-3' (outside, upstream)\n")
                f.write(f"    right flank: 5'-{rf}-3' (outside, downstream)\n")
            if len(examples) > 0:
                f.write("\n")
    print(f"\nexamples -> {OUT_EX}")


if __name__ == "__main__":
    main()
