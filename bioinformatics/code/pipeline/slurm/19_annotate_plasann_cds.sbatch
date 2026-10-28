#!/usr/bin/env bash
#SBATCH --array=0-7%8
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G
#SBATCH --job-name=nee_plasann_cds
#SBATCH --output=slurm-%x-%A_%a.out
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
OUT=${NEE_PLASANN_ANNOTATION_DIR:-$ROOT/results/19_plasann_coverage_annotation}
export MPLBACKEND=Agg MPLCONFIGDIR="$RESULTS_DIR/matplotlib_plasann_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$OUT" "$MPLCONFIGDIR"
"$PYTHON_BIN" "$PIPELINE/scripts/19_annotate_plasann_cds.py" \
    --output "$OUT" --workers 10 --blast-threads 4 \
    --shard-index "$SLURM_ARRAY_TASK_ID" --shard-count 8
