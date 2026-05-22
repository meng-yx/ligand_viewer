#!/bin/bash
#SBATCH --job-name=jupyter
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=56000
#SBATCH --time=04:00:00
#SBATCH --output=logs/jupyter-%j.out
#SBATCH --error=logs/jupyter-%j.out

echo "Starting Jupyter job on $(hostname)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"

# --- environment ---
REPO_ROOT=$(pwd)

# --- thread control ---
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --- pick a port ---
PORT=$(python -c "
import socket
s=socket.socket()
s.bind(('',0))
print(s.getsockname()[1])
s.close()
")

# --- activate conda env and launch Jupyter ---
source ~/.bashrc
conda activate MaSIF

jupyter notebook \
    --no-browser \
    --ip=0.0.0.0 \
    --port=${PORT} \
    --ServerApp.allow_remote_access=True \
    --notebook-dir=$REPO_ROOT