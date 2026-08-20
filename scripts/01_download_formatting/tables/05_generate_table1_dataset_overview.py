#!/usr/bin/env python3

import gzip
import os
import pandas as pd
from pathlib import Path

#################################################
# PATHS
#################################################

PROJECT_DIR = Path(
    os.environ.get("CVP_PROJECT_DIR", "/path/to/cardiovascular_pleiotropies")
)

FINAL_SUMSTATS = (
    PROJECT_DIR /
    "0-Download/5-Filtered/FinalSumstats"
)

OUTPUT_FILE = (
    PROJECT_DIR /
    "FinalTables/Table1_dataset_overview.csv"
)

#################################################
# PHENOTYPE METADATA
#################################################

metadata = {
    "CAD_d": {
        "category": "Disease",
        "abbr": "CAD",
        "domain": "Clinical disease",
        "phenotype": "Coronary artery disease (CAD)",
        "source": "CARDIoGRAMplusC4D",
        "ancestry": "European",
        "N": 86995,
        "Cases": 22233,
        "Controls": 64762
    },

    "T2D_d": {
        "category": "Disease",
        "abbr": "T2D",
        "phenotype": "Type 2 diabetes (T2D)",
        "domain": "Clinical disease",
        "source": "DIAGRAM",
        "ancestry": "European",
        "N": 441894,
        "Cases": 18197,
        "Controls": 423697
    },

    "HT_d": {
        "category": "Disease",
        "abbr": "HT",
        "phenotype": "Hypertension (HT)",
        "domain": "Clinical disease",
        "source": "GERA",
        "ancestry": "European",
        "N": 28391,
        "Cases": 28246,
        "Controls": None
    },

    "STR_d": {
        "category": "Disease",
        "abbr": "STR",
        "phenotype": "Stroke (STR)",   
        "domain": "Clinical disease",
        "source": "MEGASTROKE",
        "ancestry": "European",
        "N": 446696,
        "Cases": 40585,
        "Controls": 406111
    },

    "BMI_t": {
        "category": "Quantitative trait",
        "abbr": "BMI",
        "domain": "Anthropometric Trait",
        "phenotype": "Body mass index (BMI)",
        "source": "Neale Lab UK Biobank",
        "ancestry": "European",
        "N": 359983
    },

    "WC_t": {
        "category": "Quantitative trait",
        "abbr": "WC",
        "domain": "Anthropometric Trait",
        "phenotype": "Waist circumference (WC)",
        "source": "Neale Lab UK Biobank",
        "ancestry": "European",
        "N": 360564
    },

    "SBP_t": {
        "category": "Quantitative trait",
        "abbr": "SBP",
        "domain": "Blood pressure measure",
        "phenotype": "Systolic blood pressure (SBP)",
        "source": "Neale Lab UK Biobank",
        "ancestry": "European",
        "N": 340159
    },

    "DBP_t": {
        "category": "Quantitative trait",
        "abbr": "DBP",
        "domain": "Blood pressure measure",
        "phenotype": "Diastolic blood pressure",
        "source": "Neale Lab UK Biobank",
        "ancestry": "European",
        "N": 340162
    },

    "TGL_t": {
        "category": "Quantitative trait",
        "abbr": "TG",
        "domain": "Metabolic marker",
        "phenotype": "Triglycerides",
        "source": "Neale Lab UK Biobank",
        "ancestry": "European",
        "N": 343991
    },

    "LDL_t": {
        "category": "Quantitative trait",
        "abbr": "LDL",
        "domain": "Metabolic marker",
        "phenotype": "LDL cholesterol (LDL)",
        "source": "Neale Lab UK Biobank",
        "ancestry": "European",
        "N": 343621
    },

    "HDL_t": {
        "category": "Quantitative trait",
        "abbr": "HDL",
        "phenotype": "HDL cholesterol (HDL)",
        "domain": "Metabolic marker",
        "source": "Neale Lab UK Biobank",
        "ancestry": "European",
        "N": 315133
    },

    "FG_t": {
        "category": "Quantitative trait",
        "abbr": "FG",
        "domain": "Metabolic marker",
        "phenotype": "Fasting glucose (FG)",
        "source": "Neale Lab UK Biobank",
        "ancestry": "European",
        "N": 314914
    }
}

#################################################
# COUNT VARIANTS
#################################################

def count_variants(gwas_file):

    with gzip.open(gwas_file, "rt") as f:

        # subtract header
        n = sum(1 for _ in f) - 1

    return n

#################################################
# BUILD TABLE
#################################################

rows = []

# IMPORTANT:
# your files are named like:
# CAD_d-filtered.tsv.gz

files = sorted(FINAL_SUMSTATS.glob("*-filtered.tsv.gz"))

print("\nFound {} GWAS files\n".format(len(files)))

for file in files:

    phenotype = file.name.replace("-filtered.tsv.gz", "")

    print("Processing {}".format(phenotype))

    if phenotype not in metadata:
        print("WARNING: {} not in metadata".format(phenotype))
        continue

    info = metadata[phenotype]

    n_variants = count_variants(file)

    rows.append({
        "Phenotype category": info["category"],
        "Phenotype domain": info["domain"],
        "Phenotype": info["phenotype"],
        "Source / consortium": info["source"],
        "Ancestry": info["ancestry"],
        "Sample size (N)": info["N"],
        "N Cases": info.get("Cases"),
        "N Controls": info.get("Controls"),
        "N variants post-harmonization": n_variants
    })

#################################################
# EXPORT
#################################################

if len(rows) == 0:
    raise ValueError(
        "No rows generated. Check filenames and paths."
    )

table1 = pd.DataFrame(rows)

table1 = table1.sort_values(
    by=["Phenotype category", "Phenotype"]
)

# create output directory if needed
OUTPUT_DIR = PROJECT_DIR / "FinalTables"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

#################################################
# EXPORT CSV
#################################################

csv_file = OUTPUT_DIR / "Table1_dataset_overview.csv"

table1.to_csv(csv_file, index=False)

#################################################
# EXPORT EXCEL
#################################################

excel_file = OUTPUT_DIR / "Table1_dataset_overview.xlsx"

table1.to_excel(
    excel_file,
    index=False,
    sheet_name="Table1",
    engine="openpyxl"
)

#################################################
# EXPORT LATEX
#################################################

latex_file = OUTPUT_DIR / "Table1_dataset_overview.tex"

latex_table = table1.to_latex(
    index=False,
    escape=False,
    longtable=False
)

with open(latex_file, "w") as f:
    f.write(latex_table)

#################################################
# DONE
#################################################

print("\nTable 1 generated successfully\n")

print("CSV:")
print(csv_file)

print("\nExcel:")
print(excel_file)

print("\nLaTeX:")
print(latex_file)
