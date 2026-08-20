# External Step 03: LDSC Munging

LDSC `munge_sumstats.py` was used as third-party software and is not
redistributed in this code archive.

Input from previous step:

- Formatted GWAS summary statistics from `02_format_sumstats.py`.

Original command pattern:

```bash
munge_sumstats.py \
  --sumstats <formatted_gwas> \
  --out <munged_output_prefix> \
  --merge-alleles <1kgPhase3_SNPA1A2.ref.gz> \
  --snp RSID \
  --p PVAL \
  --N-col N \
  --maf-min 0.05
```

Disease phenotypes additionally used:

```bash
--N-cas-col N_CASES --N-con-col N_CONTROLS
```

Output required by next included script:

- LDSC-munged summary statistics, paired with the corresponding formatted file
  for `04_filter_to_ldsc_snps.py`.
