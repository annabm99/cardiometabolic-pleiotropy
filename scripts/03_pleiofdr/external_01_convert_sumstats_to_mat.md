# External Step 01: Convert Sumstats To PleioFDR MAT

PleioFDR's `python_convert/sumstats.py` utility was used as third-party code
and is not redistributed here.

Input from previous stage:

- Filtered GWAS summary statistics from
  `scripts/01_download_formatting/04_filter_to_ldsc_snps.py`.

Original command pattern:

```bash
python python_convert/sumstats.py mat \
  --sumstats <filtered_gwas> \
  --ref <9545380.ref> \
  --out <phenotype>.mat \
  --force
```

Output required by next included script:

- One PleioFDR `.mat` file per phenotype for
  `02_make_pleiofdr_pair_list.py`.
