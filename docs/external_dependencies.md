# External Dependencies

The analysis relies on software and reference datasets that are not distributed with this repository.

## Software

- **LDSC** — SNP heritability and genetic correlation
- **PleioFDR** — condFDR/conjFDR analyses
- **PLINK 1.9** — LD-based clumping
- **MATLAB** — execution of PleioFDR
- **FUMA** — SNP2GENE and GENE2FUNC analyses

The Python scripts additionally use standard scientific packages including
`pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `scikit-learn`, `openpyxl`, and `pyliftover`.

## Reference data

External reference resources include:

- 1000 Genomes Project Phase 3 variant and allele references
- European LD scores for LDSC
- PleioFDR reference files
- PLINK genotype reference files for clumping

GWAS summary statistics and these reference datasets are not redistributed in this repository.

## PleioFDR utilities

The analysis also uses utilities distributed with PleioFDR, including `sumstats.py` and `fdrmat2csv.py`. These are treated as third-party dependencies and are not included here.

## FUMA

FUMA analyses were performed externally using the FUMA web platform. Scripts in this repository prepare FUMA inputs and process outputs such as:

- `GenomicRiskLoci.txt`
- `leadSNPs.txt`
- `snps.txt`
- `genes.txt`
- `GS.txt`

More detailed commands for external steps are documented in the corresponding `external_*.md` files within each pipeline stage.