#!/bin/bash
set -eo pipefail
test -f validation/VERIFIED
test -f validation/SMOKE_COMPLETE
sbatch --parsable --job-name=kennedy_drcmd --partition="${SIM_PARTITION:-shared}" \
 --array="${SIM_TASKS:-1-320}%${SIM_CONCURRENT_TASKS:-24}" --cpus-per-task=1 \
 --mem="${SIM_MEMORY:-3G}" --time="${SIM_WALLTIME:-02:00:00}" \
 --exclude="${SIM_EXCLUDE:-compute-124}" \
 --output=results/logs/run_%A_%a.out --error=results/logs/run_%A_%a.err scripts/slurm.sh
