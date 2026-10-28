#!/usr/bin/env bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --job-name=nee_strict_sensitivity
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=6:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${NEE_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}
PIPELINE="$ROOT/code/pipeline"
DATA_DIR=${NEE_DATA_DIR:-$ROOT/work/raw}
RESULTS_DIR=${NEE_RESULTS_DIR:-$ROOT/work/processed}
export NEE_PROJECT_ROOT="$ROOT" NEE_DATA_DIR="$DATA_DIR" NEE_RESULTS_DIR="$RESULTS_DIR"
PYTHON_BIN=${PYTHON_BIN:-python}

"$PYTHON_BIN" "$PIPELINE/scripts/13_strict_analysis.py" \
    "${SLURM_CPUS_PER_TASK:-1}" 90 60
"$PYTHON_BIN" "$PIPELINE/scripts/14_strict_counts.py"
