# Cardiometabolic Pleiotropy Analysis

Analysis code accompanying the manuscript:

**“[MANUSCRIPT TITLE]”**

**Analysis code:** Anna Basquet Muniente

This repository contains the code used to investigate shared genetic architecture
and pleiotropy across four cardiometabolic diseases—coronary artery disease,
hypertension, stroke, and type 2 diabetes—and eight quantitative risk factors:
BMI, waist circumference, systolic and diastolic blood pressure, fasting glucose,
triglycerides, HDL cholesterol, and LDL cholesterol.

## Analysis workflow

The pipeline comprises six main stages:

1. **GWAS formatting and quality control**
2. **SNP heritability and genetic correlation** using LDSC
3. **Pleiotropic locus discovery** using condFDR/conjFDR
4. **Pleiotropic directionality** based on concordant and discordant effect signs
5. **Multi-disease convergence and functional mapping** using FUMA
6. **Functional enrichment and pathway analysis**

Scripts are organized by analysis stage and numbered according to their
execution or dependency order.

## Repository structure

```text
.
├── config/
├── docs/
│   ├── external_dependencies.md
│   └── pipeline_overview.md
└── scripts/
    ├── 01_download_formatting/
    ├── 02_ldsc/
    ├── 03_pleiofdr/
    ├── 04_directionality/
    ├── 05_fuma/
    ├── 06_enrichment/
    └── slurm/