#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=nee_count_qc
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=4:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${NEE_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}
PIPELINE="$ROOT/code/pipeline"
DATA_DIR=${NEE_DATA_DIR:-$ROOT/work/raw}
RESULTS_DIR=${NEE_RESULTS_DIR:-$ROOT/work/processed}
export NEE_PROJECT_ROOT="$ROOT" NEE_DATA_DIR="$DATA_DIR" NEE_RESULTS_DIR="$RESULTS_DIR"
PYTHON_BIN=${PYTHON_BIN:-python}

"$PYTHON_BIN" "$PIPELINE/scripts/05_compare_transposase_counts.py"
"$PYTHON_BIN" "$PIPELINE/scripts/06_discordant_breakdown.py"
