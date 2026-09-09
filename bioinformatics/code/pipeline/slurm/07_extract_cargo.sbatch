#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=nee_extract_cargo
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=6:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${NEE_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}
PYTHON_BIN=${PYTHON_BIN:-python}
export NEE_PROJECT_ROOT="$ROOT"

"$PYTHON_BIN" "$ROOT/scripts/07_extract_cargo.py" 25000
"$PYTHON_BIN" "$ROOT/scripts/07_extract_cargo.py" 25000 --strict
