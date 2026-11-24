#!/bin/bash
# Run from the project root. No original method is fitted by this array.
set -eo pipefail
test -f results/drcmd_eem_v1/verification/VERIFIED
selection="${1:-all}"
case "$selection" in
 all) task_range=1-320 ;;
 1) task_range=1-120 ;;
 2) task_range=121-240 ;;
 3) task_range=241-320 ;;
 *) echo 'Usage: bash scripts/eem/submit.sh all|1|2|3' >&2; exit 2 ;;
esac
sbatch --parsable --job-name=drcmd_eem --partition="${SIM_PARTITION:-shared}" \
 --array="${SIM_TASKS:-$task_range}%${SIM_CONCURRENT_TASKS:-24}" --cpus-per-task=1 \
 --mem="${SIM_MEMORY:-3G}" --time="${SIM_WALLTIME:-02:00:00}" \
 --exclude="${SIM_EXCLUDE:-compute-124}" \
 --output=results/drcmd_eem_v1/logs/run_%A_%a.out \
 --error=results/drcmd_eem_v1/logs/run_%A_%a.err scripts/eem/slurm.sh
