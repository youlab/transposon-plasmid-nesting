#!/usr/bin/env bash
#SBATCH --cpus-per-task=44
#SBATCH --mem=32G
#SBATCH --job-name=nee_isescan
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

ISESCAN_BIN=${ISESCAN_BIN:-isescan.py}
FNA_DIR="$RESULTS_DIR/fna"
OUTDIR="$RESULTS_DIR/04_isescan"
export OUTDIR ISESCAN_BIN
mkdir -p "$OUTDIR/tmp"
n_fna=$(find "$FNA_DIR" -maxdepth 1 -name '*.fna' | wc -l)
if [[ "$n_fna" -eq 0 ]]; then
    echo "ERROR: no FNA files in $FNA_DIR" >&2
    exit 1
fi
find "$FNA_DIR" -maxdepth 1 -name '*.fna' -print0 | \
    xargs -0 -P "${SLURM_CPUS_PER_TASK:-1}" -I{} bash -c '
        set -euo pipefail
        f="$1"
        base=$(basename "$f")
        [[ -s "$OUTDIR/fna/${base}.tsv" ]] && exit 0
        cd "$OUTDIR/tmp"
        "$ISESCAN_BIN" --seqfile "$f" --output "$OUTDIR" --nthread 1 \
            > "$OUTDIR/${base}.log" 2>&1
    ' _ {}
n_tsv=$(find "$OUTDIR/fna" -maxdepth 1 -name '*.tsv' | wc -l)
echo "ISEScan finished: $n_tsv / $n_fna TSV files"
[[ "$n_tsv" -eq "$n_fna" ]] || exit 1
