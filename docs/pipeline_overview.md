# Pipeline Overview

The analysis is organised into six main analytical stages, followed by a supplementary-table generation stage.

| Stage | Analysis | Main output |
| --- | --- | --- |
| `01_download_formatting` | GWAS formatting, harmonization and QC | Filtered summary statistics |
| `02_ldsc` | SNP heritability and genetic correlation | LDSC estimates |
| `03_pleiofdr` | condFDR/conjFDR pleiotropy analysis | Significant pleiotropic loci |
| `04_directionality` | Concordant/discordant effect classification | Directionality summaries |
| `05_fuma` | Clumping, convergence analysis and FUMA preparation | Mapped loci and genes |
| `06_enrichment` | Functional annotation and pathway enrichment | Functional summaries |
| `07_supplement` | Assembly of publication supplementary tables | Supplementary tables |

Some steps rely on external software rather than scripts included in this
repository. These transitions are marked by `external_*.md` or `manual_*.md`
files within the relevant stage.

For software requirements and reference datasets, see
[`external_dependencies.md`](external_dependencies.md).
