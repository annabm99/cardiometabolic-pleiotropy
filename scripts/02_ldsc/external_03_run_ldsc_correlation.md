# External Step 03: LDSC Genetic Correlation

LDSC `ldsc.py` was used as third-party software and is not redistributed here.

Input from previous step:

- Pair list from `01_make_ldsc_pair_list.py`.
- LDSC-munged summary statistics from stage 01.

Original command pattern:

```bash
ldsc.py \
  --rg <trait1.sumstats.gz,trait2.sumstats.gz> \
  --ref-ld-chr <eur_w_ld_chr/> \
  --w-ld-chr <eur_w_ld_chr/> \
  --out <phenotype1>_vs_<phenotype2>-GenCorr
```

Output required by next included script:

- LDSC genetic-correlation `.log` files for `05_parse_correlation_logs.py`.
