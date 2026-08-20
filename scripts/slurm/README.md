# Slurm Execution

The original analyses were run on a Slurm-based computing cluster. Several stages were parallelized using job arrays, with one phenotype, phenotype pair or analysis directory processed per task.

Generic helpers are provided here to illustrate this execution pattern:

- `submit_array.sh`
- `run_array_item.sh`

The selected input for each array task is made available through
`CVP_ARRAY_ITEM`.

Example:

```bash
bash scripts/slurm/submit_array.sh input_list.txt \
  python path/to/script.py '$CVP_ARRAY_ITEM'