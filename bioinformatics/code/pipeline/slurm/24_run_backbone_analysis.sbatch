#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --job-name=nee_backbone_analysis
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=2-00:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${NEE_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}
PIPELINE="$ROOT/code/pipeline"
DATA_DIR=${NEE_DATA_DIR:-$ROOT/work/raw}
RESULTS_DIR=${NEE_RESULTS_DIR:-$ROOT/work/processed}
export NEE_PROJECT_ROOT="$ROOT" NEE_DATA_DIR="$DATA_DIR" NEE_RESULTS_DIR="$RESULTS_DIR"
PYTHON_BIN=${PYTHON_BIN:-python}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
OUTPUT_DIR="$ROOT/results"
mkdir -p "$OUTPUT_DIR"
"$PYTHON_BIN" "$PIPELINE/scripts/24_run_backbone_analysis.py" --output-dir "$OUTPUT_DIR"
