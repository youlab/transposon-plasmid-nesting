#!/usr/bin/env bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --job-name=nee_functional_annotation
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=1-00:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${NEE_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}
DATA_DIR=${NEE_DATA_DIR:-$ROOT/1_Data}
RESULTS_DIR=${NEE_RESULTS_DIR:-$ROOT/2_Result}
DB_DIR=${NEE_DATABASE_DIR:-$DATA_DIR/database}
PYTHON_BIN=${PYTHON_BIN:-python}
DIAMOND_BIN=${DIAMOND_BIN:-diamond}
AMRFINDER_BIN=${AMRFINDER_BIN:-amrfinder}
THREADS=${SLURM_CPUS_PER_TASK:-1}
export NEE_PROJECT_ROOT="$ROOT"
export NEE_DATA_DIR="$DATA_DIR"
export NEE_RESULTS_DIR="$RESULTS_DIR"
export NEE_DATABASE_DIR="$DB_DIR"

CARD_FASTA=${CARD_FASTA:-$DB_DIR/card_protein_homolog_model.fasta}
VFDB_FASTA=${VFDB_FASTA:-$DB_DIR/VFDB_setB_pro.fas}
BACMET_FASTA=${BACMET_FASTA:-$DB_DIR/BacMet2_EXP.fasta}
TADB_FASTA=${TADB_FASTA:-$DB_DIR/TADB2_protein.fasta}

declare -A FASTAS=(
    [card]="$CARD_FASTA"
    [vfdb]="$VFDB_FASTA"
    [bacmet]="$BACMET_FASTA"
    [tadb]="$TADB_FASTA"
)

"$PYTHON_BIN" "$ROOT/scripts/08_dedup_proteins.py"

for name in card vfdb bacmet tadb; do
    fasta=${FASTAS[$name]}
    database="$DB_DIR/$name"
    [[ -s "$fasta" ]] || { echo "Missing database FASTA: $fasta" >&2; exit 1; }
    if [[ ! -s "${database}.dmnd" || "$fasta" -nt "${database}.dmnd" ]]; then
        "$DIAMOND_BIN" makedb --in "$fasta" -d "$database"
    fi
done

for setname in Cargo Background; do
    for name in card vfdb bacmet tadb; do
        output="$RESULTS_DIR/6_${setname}_vs_${name}.tsv"
        if [[ ! -s "$output" ]]; then
            "$DIAMOND_BIN" blastp --very-sensitive --threads "$THREADS" \
                -q "$RESULTS_DIR/6_${setname}_unique.faa" \
                -d "$DB_DIR/$name" -o "$output" \
                --outfmt 6 qseqid sseqid pident length qcovhsp evalue bitscore stitle \
                -e 1e-10 -k 1
        fi
    done
done

AMRFINDER_THREADS=${AMRFINDER_THREADS:-8}
for setname in Cargo Background; do
    output="$RESULTS_DIR/6_AMRF_${setname}.tsv"
    [[ -s "$output" ]] && continue
    chunk_dir="$RESULTS_DIR/6_AMRF_${setname}_chunks"
    mkdir -p "$chunk_dir"
    split -l 100000 -d -a 3 --additional-suffix=.faa \
        "$RESULTS_DIR/6_${setname}_unique.faa" "$chunk_dir/chunk_"
    first=1
    : > "$output"
    for chunk in "$chunk_dir"/chunk_*.faa; do
        chunk_output="${chunk%.faa}.tsv"
        amrfinder_args=(-p "$chunk" --threads "$AMRFINDER_THREADS" -o "$chunk_output")
        if [[ -n "${AMRFINDER_DB:-}" ]]; then
            amrfinder_args+=(-d "$AMRFINDER_DB")
        fi
        "$AMRFINDER_BIN" "${amrfinder_args[@]}"
        if [[ "$first" -eq 1 ]]; then
            cat "$chunk_output" >> "$output"
            first=0
        else
            tail -n +2 "$chunk_output" >> "$output"
        fi
    done
    rm -rf "$chunk_dir"
done

"$PYTHON_BIN" "$ROOT/scripts/10_function_inventory.py" 90 60
"$PYTHON_BIN" "$ROOT/scripts/11_enrichment.py" 90 60
"$PYTHON_BIN" "$ROOT/scripts/11_enrichment.py" 80 70 _p80q70
"$PYTHON_BIN" "$ROOT/scripts/12_amrf_enrichment.py"
