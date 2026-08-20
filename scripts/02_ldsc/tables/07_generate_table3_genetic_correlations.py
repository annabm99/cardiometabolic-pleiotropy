#!/usr/bin/env python3

import os
import pandas as pd
from statsmodels.stats.multitest import multipletests
from openpyxl import load_workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter


# =========================
# Input / Output
# =========================

PROJECT_DIR = os.environ.get(
    "CVP_PROJECT_DIR",
    "/path/to/cardiovascular_pleiotropies"
)

INPUT_FILE = os.environ.get(
    "CVP_LDSC_CORRELATION_DF",
    f"{PROJECT_DIR}/1-LDSC/2-Correlation/JointData/PairsDf.gz"
)

OUTPUT_DIR = os.environ.get(
    "CVP_FINAL_TABLES_DIR",
    f"{PROJECT_DIR}/FinalTables"
)
BASE_NAME = "Table3_GeneticCorrelations"

OUTPUT_CSV = f"{OUTPUT_DIR}/{BASE_NAME}.csv"
OUTPUT_XLSX = f"{OUTPUT_DIR}/{BASE_NAME}.xlsx"
OUTPUT_LATEX = f"{OUTPUT_DIR}/{BASE_NAME}.tex"


# =========================
# Phenotype name mapping
# =========================

PHENO_MAP = {
    "CAD_d": "Coronary artery disease",
    "T2D_d": "Type 2 diabetes",
    "HT_d": "Hypertension",
    "STR_d": "Stroke",
    "BMI_t": "Body mass index",
    "WC_t": "Waist circumference",
    "SBP_t": "Systolic blood pressure",
    "DBP_t": "Diastolic blood pressure",
    "FG_t": "Fasting glucose",
    "TGL_t": "Triglycerides",
    "HDL_t": "HDL cholesterol",
    "LDL_t": "LDL cholesterol"
}


# =========================
# Load data
# =========================

print("Loading pairwise genetic correlation table...")


df = pd.read_csv(INPUT_FILE, sep="\t")


# =========================
# Recalculate multiple testing corrections
# =========================

print("Applying multiple testing correction...")

pvals = df["p-value"]

fdr_corrected = multipletests(pvals, method="fdr_bh")[1]
bonf_corrected = multipletests(pvals, method="bonferroni")[1]


df["p_FDR"] = fdr_corrected

df["p_Bonferroni"] = bonf_corrected


# =========================
# Rename phenotypes
# =========================

print("Mapping phenotype names...")


df["Phenotype1"] = df["Phenotype1"].map(PHENO_MAP)
df["Phenotype2"] = df["Phenotype2"].map(PHENO_MAP)

# =========================
# Add significance columns
# =========================

print("Adding significance columns...")


df["FDR significant"] = df["p_FDR"] < 0.05

df["Bonferroni significant"] = df["p_Bonferroni"] < 0.05


# =========================
# Rename columns for publication
# =========================

print("Formatting manuscript-ready columns...")


df = df.rename(columns={
    "Phenotype1": "Phenotype 1",
    "Phenotype2": "Phenotype 2",
    "rg": "Genetic correlation (rg)",
    "se": "Standard error",
    "p-value": "P-value",
    "p_FDR": "FDR-adjusted p-value",
    "p_Bonferroni": "Bonferroni-adjusted p-value"
})


# =========================
# Sort by significance
# =========================

print("Sorting by P-value...")


df = df.sort_values("P-value")


# =========================
# Round numeric columns
# =========================

print("Rounding numeric columns...")

round_cols = [
    "Genetic correlation (rg)",
    "Standard error"
]

for col in round_cols:
    df[col] = df[col].round(4)


# =========================
# Export TSV
# =========================

print(f"Writing TSV: {OUTPUT_CSV}")


df.to_csv(OUTPUT_CSV, index=False)


# =========================
# Export LaTeX
# =========================

print(f"Writing LaTeX: {OUTPUT_LATEX}")


latex_df = df.copy()

latex_df.to_latex(
    OUTPUT_LATEX,
    index=False,
    escape=False,
    float_format="%.4g",
    longtable=True
)


# =========================
# Export Excel
# =========================

print(f"Writing Excel: {OUTPUT_XLSX}")


df.to_excel(OUTPUT_XLSX, index=False)


# =========================
# Excel formatting
# =========================

print("Applying Excel formatting...")


wb = load_workbook(OUTPUT_XLSX)
ws = wb.active


# Freeze top row
ws.freeze_panes = "A2"


# Autofilter
ws.auto_filter.ref = ws.dimensions


# Bold header
for cell in ws[1]:
    cell.font = Font(bold=True)


# Adjust column widths
for column_cells in ws.columns:
    length = max(len(str(cell.value)) if cell.value is not None else 0 for cell in column_cells)
    adjusted_width = min(length + 4, 50)
    column_letter = get_column_letter(column_cells[0].column)
    ws.column_dimensions[column_letter].width = adjusted_width


# Scientific notation for p-values
pval_columns = [
    "E",  # P-value
    "F",  # FDR-adjusted
    "G"   # Bonferroni-adjusted
]

for col in pval_columns:
    for cell in ws[col][1:]:
        cell.number_format = '0.00E+00'


# Format rg and SE
numeric_columns = ["C", "D"]

for col in numeric_columns:
    for cell in ws[col][1:]:
        cell.number_format = '0.0000'


wb.save(OUTPUT_XLSX)


print("Done!")
print(f"Generated files:")
print(f" - {OUTPUT_CSV}")
print(f" - {OUTPUT_XLSX}")
print(f" - {OUTPUT_LATEX}")
