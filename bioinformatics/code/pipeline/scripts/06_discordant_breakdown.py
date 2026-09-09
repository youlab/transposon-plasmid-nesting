#!/usr/bin/env python3
"""
06_discordant_breakdown.py  (count-discordance QC)
For the plasmids where Rohan's annotation screen found >=1 transposase but
ISEScan found ZERO IS elements (the annotation-positive / sequence-negative
discordant class), extract the regex-matched CDS product strings from their
gbff and classify them:
  - tn3_family_unit_transposon   (outside ISEScan's IS-only scope)
  - other_transposase_or_IS      (plausible real TE annotation; ISEScan miss
                                  or degraded element)
  - unclear_likely_false_positive (regex fired on an unrelated product)

Also validates our reproduction of Rohan's regex: per-plasmid matched-CDS
count vs his transposase_count (same gbff products, same regex).

Defaults (relative to the repository root):
    Input:  2_Result/4_ISEScan/fna/*.fna.tsv
            2_Result/3_Seqhash_map.csv
            1_Data/transposase-counts.csv
            1_Data/gbff/<acc>.gbff
    Output: 2_Result/4c_Discordant_cds.csv
            2_Result/4c_Summary.txt

Usage:  python scripts/06_discordant_breakdown.py
"""
import csv
import glob
import os
import re
from collections import Counter

from Bio import SeqIO

ROOT = os.path.abspath(os.environ.get(
    "NEE_PROJECT_ROOT", os.path.dirname(os.path.dirname(__file__))))
DATA_DIR = os.path.abspath(os.environ.get(
    "NEE_DATA_DIR", os.path.join(ROOT, "1_Data")))
RESULTS_DIR = os.path.abspath(os.environ.get(
    "NEE_RESULTS_DIR", os.path.join(ROOT, "2_Result")))
TSV_DIR = os.path.join(RESULTS_DIR, "4_ISEScan", "fna")
MAP_CSV = os.path.join(RESULTS_DIR, "3_Seqhash_map.csv")
ROHAN_CSV = os.path.join(DATA_DIR, "transposase-counts.csv")
GBFF_DIR = os.path.join(DATA_DIR, "gbff")
OUT_CDS = os.path.join(RESULTS_DIR, "4c_Discordant_cds.csv")
OUT_SUM = os.path.join(RESULTS_DIR, "4c_Summary.txt")

# Rohan's screen regex (exact string from his script, cf. archive §8)
ROHAN_RE = re.compile(
    r"^IS|transpos\S*|insertion|conjugate transposon|Transpos\S*|Tn[0-9]"
    r"|tranposase|Tnp|Ins|ins")

# Tn3-family unit transposons (Tn1/2/3 are the same family; Tn21/Tn501/etc.
# are its best-known members; tnpA/tnpR are its transposase/resolvase genes).
# Numbered Tn like Tn5/Tn7/Tn10 are NOT Tn3-family (they are IS-composites).
TN3_RE = re.compile(
    r"Tn\s?1\b|Tn\s?2\b|Tn\s?3\b|Tn21|Tn501|Tn1721|Tn1331|Tn5393|Tn4400"
    r"|Tn4401|tnpA|tnpR", re.I)

# tokens that make a matched product still plausibly TE-related
TEISH_RE = re.compile(
    r"transpos|insertion|tnp|integrase|resolv|IS[0-9A-Z]|ins[ABEHJ]\b", re.I)


def tsv_copy_count(path):
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith("seqID") or not line.strip():
                continue
            if len(line.split("\t")) >= 22:
                n += 1
    return n


def classify(product):
    if TN3_RE.search(product):
        return "tn3_family_unit_transposon"
    if TEISH_RE.search(product):
        return "other_transposase_or_IS"
    return "unclear_likely_false_positive"


def main():
    # discordant set: md5 with 0-copy ISEScan tsv
    discordant = {}
    for path in glob.glob(os.path.join(TSV_DIR, "*.fna.tsv")):
        if tsv_copy_count(path) == 0:
            md5 = os.path.basename(path)[:-len(".fna.tsv")]
            discordant[md5] = True

    md5rep = {}
    with open(MAP_CSV) as f:
        for row in csv.DictReader(f):
            if row["md5"] in discordant:
                md5rep[row["md5"]] = row["representative_acc"]

    rohan = {}
    with open(ROHAN_CSV) as f:
        r = csv.reader(f)
        next(r)
        for rec in r:
            if rec[2] == "plasmid":
                rohan[rec[1]] = int(rec[6])

    cds_rows = []
    plasmid_cats = {}   # acc -> set of categories
    n_gbff_missing = 0
    n_count_match = 0
    prod_counter = Counter()
    cat_prod = {"tn3_family_unit_transposon": Counter(),
                "other_transposase_or_IS": Counter(),
                "unclear_likely_false_positive": Counter()}

    for md5, acc in sorted(md5rep.items()):
        gbff = os.path.join(GBFF_DIR, acc + ".gbff")
        if not os.path.exists(gbff):
            n_gbff_missing += 1
            continue
        rec = SeqIO.read(gbff, "genbank")
        n_hit = 0
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            product = feat.qualifiers.get("product", [""])[0]
            m = ROHAN_RE.search(product)
            if not m:
                continue
            cat = classify(product)
            n_hit += 1
            plasmid_cats.setdefault(acc, set()).add(cat)
            locus = feat.qualifiers.get("locus_tag", ["NA"])[0]
            cds_rows.append([acc, md5, locus, product, m.group(0), cat,
                             rohan.get(acc, "")])
            prod_counter[product] += 1
            cat_prod[cat][product] += 1
        if rohan.get(acc) is not None and n_hit == rohan[acc]:
            n_count_match += 1

    with open(OUT_CDS, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["accession", "md5", "locus_tag", "product",
                    "regex_match", "category", "rohan_transposase_count"])
        w.writerows(cds_rows)

    # plasmid-level category: tn3 > other > unclear; plasmids with NO matched
    # CDS in the current gbff are reported separately (annotation drift since
    # Rohan's screen snapshot).
    pl_cat = {"tn3_family_unit_transposon": 0,
              "other_transposase_or_IS": 0,
              "unclear_likely_false_positive": 0,
              "no_match_in_current_gbff": 0}
    for acc in md5rep.values():
        cats = plasmid_cats.get(acc, set())
        for c in ("tn3_family_unit_transposon", "other_transposase_or_IS",
                  "unclear_likely_false_positive"):
            if c in cats:
                pl_cat[c] += 1
                break
        else:
            pl_cat["no_match_in_current_gbff"] += 1

    cds_cat = Counter(r[5] for r in cds_rows)
    n_pl = len(md5rep)
    lines = [
        "# 4c: breakdown of annotation-positive / ISEScan-zero plasmids",
        f"discordant_plasmids_unique: {n_pl}  (gbff missing: {n_gbff_missing})",
        f"total_regex_matched_cds: {len(cds_rows)}",
        f"rohan_count_reproduced_exactly: {n_count_match}/{n_pl} plasmids "
        f"(matched CDS == his transposase_count)",
        "",
        "## plasmid-level (each plasmid assigned to its most specific "
        "category: tn3 > other > unclear)",
    ]
    lines += [f"  {k}: {v} ({v / max(n_pl, 1):.1%})" for k, v in pl_cat.items()]
    lines += ["", "## CDS-level"]
    lines += [f"  {k}: {cds_cat.get(k, 0)}" for k in pl_cat
              if k != "no_match_in_current_gbff"]
    for k in cat_prod:
        lines += ["", f"## top products: {k}"]
        lines += [f"  {n}\t{p}" for p, n in cat_prod[k].most_common(10)]
    lines += ["", "## top 20 products overall"]
    lines += [f"  {n}\t{p}" for p, n in prod_counter.most_common(20)]

    with open(OUT_SUM, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines[:14]))
    print(f"...\ncds_csv -> {OUT_CDS}\nsummary -> {OUT_SUM}")


if __name__ == "__main__":
    main()
