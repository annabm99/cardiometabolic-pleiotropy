# External Step 04: Run PleioFDR In MATLAB

The PleioFDR MATLAB implementation, including `runme.m`, was used as
third-party software and is not redistributed here.

Input from previous steps:

- Pair list from `02_make_pleiofdr_pair_list.py`.
- Pair-specific config files created with `03_modify_pleiofdr_config.py`.
- PleioFDR `.mat` files from `external_01_convert_sumstats_to_mat.md`.

Original command pattern:

```bash
matlab -nodisplay -nosplash -nodesktop \
  -r "config='<config-file>'; run('runme.m'); exit;"
```

Output required by next included script:

- Pair-specific PleioFDR result files for
  `05_extract_pleiotropic_snps.py`.
