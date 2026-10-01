#!/bin/bash
#SBATCH --job-name=am_rer
#SBATCH --partition=shortrun
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=0-9
#SBATCH --output=logs/rer_%A_%a.out
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
N=40000; C=10; PER=$((N / C))
export UNIT_START=$((SLURM_ARRAY_TASK_ID * PER)) UNIT_END=$(((SLURM_ARRAY_TASK_ID + 1) * PER))
export WORKERS=16 RESULTS_DIR=results_rer CHUNK_ID=$SLURM_ARRAY_TASK_ID
cd "${SLURM_SUBMIT_DIR:-.}"
Rscript sim_rer.R
