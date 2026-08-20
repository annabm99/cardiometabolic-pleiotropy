#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: run_array_item.sh INPUT_LIST COMMAND [ARG ...]" >&2
    exit 1
fi

INPUT_LIST=$1
shift

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is not set; run this script through sbatch --array." >&2
    exit 1
fi

CVP_ARRAY_ITEM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_LIST")
export CVP_ARRAY_ITEM

if [[ -z "$CVP_ARRAY_ITEM" ]]; then
    echo "No input item found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "Array task ${SLURM_ARRAY_TASK_ID}: ${CVP_ARRAY_ITEM}"

exec "$@"
