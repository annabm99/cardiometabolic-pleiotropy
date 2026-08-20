# External Step 02: LDSC Heritability

LDSC `ldsc.py` was used as third-party software and is not redistributed here.

Input from previous step:

- LDSC-munged summary statistics from stage 01.

Original command pattern:

```bash
ldsc.py \
  --h2 <munged_sumstats> \
  --ref-ld-chr <eur_w_ld_chr/> \
  --w-ld-chr <eur_w_ld_chr/> \
  --out <heritability_output_prefix>
```

Output required by next included script:

- LDSC heritability `.log` files for `04_parse_heritability_logs.py`.
