#!/bin/bash
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
module load conda_R/4.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript run.R eem task "$SLURM_ARRAY_TASK_ID"
