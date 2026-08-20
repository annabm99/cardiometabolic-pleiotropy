#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: submit_array.sh INPUT_LIST COMMAND [ARG ...]" >&2
    echo "Example:" >&2
    echo "  submit_array.sh inputs.txt bash -lc 'python script.py \"$CVP_ARRAY_ITEM\" /path/to/out'" >&2
    exit 1
fi

INPUT_LIST=$1
shift

if [[ ! -f "$INPUT_LIST" ]]; then
    echo "Input list not found: $INPUT_LIST" >&2
    exit 1
fi

JOBS_COUNT=$(grep -cve '^[[:space:]]*$' "$INPUT_LIST")

if [[ "$JOBS_COUNT" -lt 1 ]]; then
    echo "Input list is empty: $INPUT_LIST" >&2
    exit 1
fi

echo "Submitting ${JOBS_COUNT} Slurm array jobs"

sbatch \
    --array=1-"${JOBS_COUNT}" \
    "$(dirname "$0")/run_array_item.sh" \
    "$INPUT_LIST" \
    "$@"
