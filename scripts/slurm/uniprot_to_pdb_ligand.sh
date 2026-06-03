#!/bin/bash
#SBATCH --job-name=uniprot_ligands
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=112000
#SBATCH --time=12:00:00
#SBATCH --output=logs/uniprot_ligands-%j.out
#SBATCH --error=logs/uniprot_ligands-%j.out



if [ "$#" -ne 3 ]; then
  echo "Usage: sbatch $0 <input_txt> <output_csv> <blacklist_txt>"
  exit 1
fi

INPUT_TXT="$1"
OUTPUT_CSV="$2"
BLACKLIST_TXT="$3"

echo "Starting UniProt ligand mapping on $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "INPUT_TXT=${INPUT_TXT}"
echo "OUTPUT_CSV=${OUTPUT_CSV}"
echo "BLACKLIST_TXT=${BLACKLIST_TXT}"

# --- environment ---
REPO_ROOT=$(git rev-parse --show-toplevel)

# --- thread control ---
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --- activate conda env ---
set +u
source ~/.bashrc
set -u
conda activate MaSIF

echo "Running uniprot_to_pdb_ligand.py"
python "${REPO_ROOT}/scripts/python/uniprot_to_pdb_ligand.py" \
  --input "${INPUT_TXT}" \
  --output "${OUTPUT_CSV}" \
  --blacklist "${BLACKLIST_TXT}" 
