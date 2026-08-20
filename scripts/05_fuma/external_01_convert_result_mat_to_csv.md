# External Step 01: Convert PleioFDR Result MAT To CSV

PleioFDR's `python_convert/fdrmat2csv.py` utility was used as third-party code
and is not redistributed here.

Input from previous stage:

- Pair-specific PleioFDR `result.mat` files from stage 03.

Original command pattern:

```bash
python python_convert/fdrmat2csv.py \
  --mat <pair_result_dir>/result.mat \
  --ref <9545380.ref> \
  --out <pair>-FdrResults.csv
```

Output required by next external step:

- CSV-formatted PleioFDR FDR results for PLINK clumping.

Note:

- The historical wrapper also called `ClumpPrep.py`, but that script is not
  present in this local copy and should be checked against the original cluster
  project if exact clump-prep reproduction is needed.
