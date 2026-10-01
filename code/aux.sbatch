#!/bin/bash
#SBATCH --job-name=am_aux
#SBATCH --partition=shortrun
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/aux_%A_%a.out
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
N=37800; C=4; PER=$((N / C + 1))
export UNIT_START=$((SLURM_ARRAY_TASK_ID * PER)) UNIT_END=$(((SLURM_ARRAY_TASK_ID + 1) * PER))
export WORKERS=16 RESULTS_DIR=results_aux CHUNK_ID=$SLURM_ARRAY_TASK_ID
cd "${SLURM_SUBMIT_DIR:-.}"
Rscript sim_aux.R
