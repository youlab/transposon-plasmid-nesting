#!/usr/bin/env python3
"""Run the count and Fisher/BH backbone analysis without plotting.

Inputs: bundled data/19_cds_annotation_membership.csv.gz by default, or an
annotation directory with --from-annotation. Outputs go to results/ by default.
Use --plot explicitly to run the separate plotting entry point.
"""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "results"
DEFAULT_MEMBERSHIP = ROOT / "data" / "19_cds_annotation_membership.csv.gz"
DEFAULT_CONFIG = ROOT / "code" / "pipeline" / "config" / "analysis.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--membership-file", type=Path,
                        help="override the membership input selected by the package config")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--input-dir", type=Path,
                        default=DEFAULT_OUTPUT / "19_plasann_coverage_annotation",
                        help="annotation directory used with --from-annotation")
    parser.add_argument("--from-annotation", action="store_true",
                        help="rebuild membership by joining existing per-CDS annotations")
    parser.add_argument("--plot", action="store_true",
                        help="run script 23 to create figures after statistics")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    count_command = [sys.executable, str(SCRIPTS / "20_filter_and_count_cds.py"),
                     "--output-dir", str(args.output_dir), "--config", str(args.config)]
    if args.from_annotation:
        count_command.extend(("--from-annotation", "--input-dir", str(args.input_dir)))
    elif args.membership_file is not None:
        count_command.extend(("--membership-file", str(args.membership_file)))
    subprocess.run(count_command, check=True)

    stats_command = [sys.executable, str(SCRIPTS / "21_calculate_backbone_enrichment.py"),
                     "--input-dir", str(args.output_dir), "--output-dir", str(args.output_dir),
                     "--config", str(args.config)]
    subprocess.run(stats_command, check=True)

    if args.plot:
        subprocess.run([sys.executable, str(SCRIPTS / "23_plot_backbone_enrichment.py"),
                        "--input-dir", str(args.output_dir),
                        "--output-dir", str(args.output_dir / "figures"),
                        "--config", str(args.config)], check=True)


if __name__ == "__main__":
    main()
