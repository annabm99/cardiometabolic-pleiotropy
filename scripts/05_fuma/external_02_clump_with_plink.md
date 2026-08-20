# External Step 02: Clump PleioFDR Results With PLINK

PLINK 1.9 was used through PleioFDR's `sumstats.py clump` wrapper. PLINK,
PleioFDR's wrapper, and genotype reference files are external dependencies and
are not redistributed here.

Input from previous step:

- CSV-formatted PleioFDR FDR results from
  `external_01_convert_result_mat_to_csv.md`.

Original command pattern:

```bash
python python_convert/sumstats.py clump \
  --sumstats <pair>-FdrResults.csv \
  --out <pair>-FdrClumped \
  --bfile-chr <plink_reference_chr@> \
  --clump-field FDR \
  --clump-snp-field SNP \
  --clump-p1 0.05 \
  --plink plink
```

Output required by next included script:

- `<pair>-FdrClumped.snps.csv`
- `<pair>-FdrClumped.loci.csv`

These are consumed by `03_split_clumped_loci_by_direction.py`.
