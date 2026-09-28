# Cardiometabolic Pleiotropy Analysis

Analysis code accompanying a manuscript investigating shared genetic architecture
and pleiotropy across four cardiometabolic diseases—coronary artery disease,
hypertension, stroke, and type 2 diabetes—and eight quantitative risk factors:
BMI, waist circumference, systolic and diastolic blood pressure, random
(non-fasting) glucose, triglycerides, HDL cholesterol, and LDL cholesterol.

**Analysis code:** Anna Basquet Muniente

## Analysis workflow

The pipeline comprises six main analytical stages:

1. **GWAS formatting and quality control**
2. **SNP heritability and genetic correlation** using LDSC
3. **Pleiotropic locus discovery** using condFDR/conjFDR
4. **Pleiotropic directionality** based on concordant and discordant effect signs
5. **Multi-disease convergence and functional mapping** using FUMA
6. **Functional enrichment and pathway analysis**

Publication supplementary tables are generated in a separate final stage.

Scripts are organised by analysis stage and numbered according to their
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
    ├── 07_supplement/
    └── slurm/
```

## Configuration

Analysis scripts use environment variables to define project and external-data
locations. The principal variable is:

```bash
export CVP_PROJECT_DIR=/path/to/cardiovascular_pleiotropies
```

Additional stage-specific environment variables are used where required for
reference data, software installations, output directories, or external inputs.
Placeholder paths such as `/path/to/...` should be replaced with local paths
before execution.

## Dependencies

The workflow uses Python scientific-computing packages together with external
software including LDSC, PleioFDR, PLINK 1.9, MATLAB, and FUMA.

Software requirements, reference datasets, and externally executed steps are
described in [`docs/external_dependencies.md`](docs/external_dependencies.md).
An overview of the pipeline is provided in
[`docs/pipeline_overview.md`](docs/pipeline_overview.md).

## Data availability

GWAS summary statistics and external reference datasets are not redistributed
with this repository. Their sources are described in the accompanying
manuscript and publication tables.