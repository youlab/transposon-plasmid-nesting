#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --job-name=nee_download_gbff
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=2-00:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${NEE_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}
PYTHON_BIN=${PYTHON_BIN:-python}
export NEE_PROJECT_ROOT="$ROOT"

: "${NCBI_EMAIL:?Set NCBI_EMAIL before submitting this job}"
"$PYTHON_BIN" "$ROOT/scripts/02_download_gbff.py"
