#!/usr/bin/env python3

"""
04_generate_table8_pathways.py

Purpose
-------
Generate final Supplementary Table 8:
Nonredundant FUMA pathway enrichment results.

Exports:
- CSV
- XLSX
- LaTeX

Input
-----
all_fuma_pathways_nonredundant.tsv

Output
------
Table8_FUMA_Pathway_Enrichment.csv
Table8_FUMA_Pathway_Enrichment.xlsx
Table8_FUMA_Pathway_Enrichment.tex

Configuration
-------------
Set CVP_NONREDUNDANT_PATHWAYS to the nonredundant enrichment table and
CVP_FINAL_TABLES_DIR to the output directory. If unset, placeholder paths under
CVP_PROJECT_DIR are used.
"""

import pandas as pd
import os
from pathlib import Path
import sys

# ============================================================
# INPUTS
# ============================================================

PROJECT_DIR = Path(
    os.environ.get("CVP_PROJECT_DIR", "/path/to/cardiovascular_pleiotropies")
)

INPUT_FILE = Path(
    os.environ.get(
        "CVP_NONREDUNDANT_PATHWAYS",
        str(PROJECT_DIR / "5-EnrichmentAnalysis/3-Collapse/all_fuma_pathways_nonredundant.tsv")
    )
)
OUTDIR = Path(
    os.environ.get("CVP_FINAL_TABLES_DIR", str(PROJECT_DIR / "FinalTables"))
)

OUTDIR.mkdir(parents=True, exist_ok=True)

BASENAME = "Table8_FUMA_Pathway_Enrichment"

CSV_OUT = OUTDIR / f"{BASENAME}.csv"
XLSX_OUT = OUTDIR / f"{BASENAME}.xlsx"
LATEX_OUT = OUTDIR / f"{BASENAME}.tex"

# ============================================================
# LOAD
# ============================================================

print("Loading nonredundant enrichment dataset...")

df = pd.read_csv(INPUT_FILE, sep="\t")

df = df[
    ~df["phenotype_pair"].str.contains(
        "_d-HDL_|_d-LDL_",
        regex=True
    )
].copy()

print(f"Rows loaded: {len(df):,}")

# ============================================================
# SORT
# ============================================================

df = df.sort_values(
    by=[
        "phenotype_pair",
        "pleiotropy_direction",
        "fdr",
        "p_value"
    ]
)

# ============================================================
# FORMAT NUMERIC COLUMNS
# ============================================================

df["p_value"] = df["p_value"].apply(
    lambda x: f"{x:.3e}"
)

df["fdr"] = df["fdr"].apply(
    lambda x: f"{x:.3e}"
)

# ============================================================
# CLEAN COLUMN NAMES
# ============================================================

df = df.rename(columns={
    "phenotype_pair": "Phenotype_Pair",
    "pleiotropy_direction": "Pleiotropy_Direction",
    "pathway_name": "Pathway",
    "category": "Category",
    "n_overlap": "N_Overlap",
    "p_value": "P_Value",
    "fdr": "FDR"
})

# ============================================================
# EXPORT CSV
# ============================================================

df.to_csv(
    CSV_OUT,
    index=False
)

print(f"CSV saved: {CSV_OUT}")

# ============================================================
# EXPORT XLSX
# ============================================================

try:
    with pd.ExcelWriter(XLSX_OUT, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Table8")
    print(f"XLSX saved: {XLSX_OUT}")
except ImportError:
    print("WARNING: openpyxl not installed; skipping XLSX export.")

# ============================================================
# EXPORT LATEX
# ============================================================

latex_table = df.to_latex(
    index=False,
    escape=False,
    longtable=True
)

with open(LATEX_OUT, "w") as f:
    f.write(latex_table)

print(f"LaTeX saved: {LATEX_OUT}")

# ============================================================
# SUMMARY
# ============================================================

print("\n====================================================")
print("TABLE 8 GENERATED")
print("====================================================")
print(f"Rows exported: {len(df):,}")
print(f"CSV:   {CSV_OUT}")
print(f"XLSX:  {XLSX_OUT}")
print(f"LaTeX: {LATEX_OUT}")
print("====================================================")
