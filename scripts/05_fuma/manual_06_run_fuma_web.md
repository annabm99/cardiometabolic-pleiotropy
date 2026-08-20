# Manual Step 06: Run FUMA

FUMA analyses were run outside this repository. The FUMA web service/software
and generated FUMA outputs are not redistributed here.

Inputs from previous included scripts:

- Direction-specific lead and candidate SNP files from
  `04_prepare_fuma_pair_uploads.sh`.
- Full GWAS summary statistics prepared by
  `05_prepare_full_sumstats_for_fuma.py`.

Procedure:

- Upload the appropriate positive and negative directionality inputs for each
  no-HDL/LDL phenotype pair to FUMA.
- Run SNP2GENE/GENE2FUNC analyses using the same FUMA settings as the original
  publication analysis.
- Download pair-specific FUMA output folders.

Outputs required by next included scripts:

- For annotation summaries: `leadSNPs.txt`, `snps.txt`, and `genes.txt`.
- For convergence tables: `GenomicRiskLoci.txt`, `leadSNPs.txt`, `snps.txt`,
  and `genes.txt`.
- For enrichment: `GS.txt` files from GENE2FUNC.
