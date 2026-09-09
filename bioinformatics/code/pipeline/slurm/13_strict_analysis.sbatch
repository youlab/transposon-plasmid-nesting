#!/usr/bin/env bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --job-name=nee_strict_sensitivity
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=6:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${NEE_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}
PYTHON_BIN=${PYTHON_BIN:-python}
export NEE_PROJECT_ROOT="$ROOT"

"$PYTHON_BIN" "$ROOT/scripts/13_strict_analysis.py" \
    "${SLURM_CPUS_PER_TASK:-1}" 90 60
"$PYTHON_BIN" "$ROOT/scripts/14_strict_counts.py"
