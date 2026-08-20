#!/usr/bin/env python3

"""
05_generate_table9_recurrent_pathways.py

Purpose
-------
Generate final main manuscript Table 9:
Recurrent nonredundant pathway enrichments across
phenotype pairs, separated by pleiotropy direction.

This script summarizes pathway recurrence across analyses
to identify the dominant biological themes associated with:
- positive pleiotropy
- negative pleiotropy

Exports:
- CSV
- XLSX
- LaTeX

Input
-----
all_fuma_pathways_nonredundant.tsv

Output
------
Table9_Top_Pathway_Enrichment.csv
Table9_Top_Pathway_Enrichment.xlsx
Table9_Top_Pathway_Enrichment.tex

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

BASENAME = "Table9_Top_Pathway_Enrichment"

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
# SUMMARIZE PATHWAY RECURRENCE
# ============================================================

summary = (
    df.groupby(
        [
            "pleiotropy_direction",
            "pathway_name",
            "category"
        ]
    )
    .agg(
        N_Phenotype_Pairs=("phenotype_pair", "nunique"),
        Mean_N_Overlap=("n_overlap", "mean"),
        Best_P_Value=("p_value", "min"),
        Best_FDR=("fdr", "min")
    )
    .reset_index()
)

# ============================================================
# REPRESENTATIVE PHENOTYPE PAIRS
# ============================================================

disease_summary = (
    df.groupby(
        [
            "pleiotropy_direction",
            "pathway_name"
        ]
    )["phenotype_pair"]
    .apply(
        lambda x: sorted(
            set(
                pair.split("_d-")[0]
                for pair in x
            )
        )
    )
    .reset_index(name="Diseases_Involved")
)

disease_summary["N_Diseases"] = (
    disease_summary["Diseases_Involved"]
    .apply(len)
)

disease_summary["Diseases_Involved"] = (
    disease_summary["Diseases_Involved"]
    .apply(lambda x: "; ".join(x))
)

summary = summary.merge(
    disease_summary,
    on=[
        "pleiotropy_direction",
        "pathway_name"
    ],
    how="left"
)


# ============================================================
# SORT WITHIN DIRECTION
# ============================================================

summary = summary.sort_values(
    by=[
        "pleiotropy_direction",
        "N_Diseases",
        "N_Phenotype_Pairs",
        "Best_FDR"
    ],
    ascending=[True, False, False, True]
)

# ============================================================
# KEEP ALL PATHWAYS
# ============================================================

final_df = summary.copy()

# ============================================================
# FORMAT NUMERIC COLUMNS
# ============================================================

final_df["Mean_N_Overlap"] = final_df["Mean_N_Overlap"].round(2)

final_df["Best_P_Value"] = final_df["Best_P_Value"].apply(
    lambda x: f"{x:.3e}"
)

final_df["Best_FDR"] = final_df["Best_FDR"].apply(
    lambda x: f"{x:.3e}"
)

# ============================================================
# CLEAN COLUMN NAMES
# ============================================================

final_df = final_df.rename(columns={
    "pleiotropy_direction": "Pleiotropy_Direction",
    "pathway_name": "Pathway",
    "category": "Category"
})

# ============================================================
# FINAL COLUMN ORDER
# ============================================================

final_df = final_df[
    [
        "Pleiotropy_Direction",
        "Pathway",
        "Category",
        "N_Diseases",
        "Diseases_Involved",
        "N_Phenotype_Pairs",
        "Mean_N_Overlap",
        "Best_P_Value",
        "Best_FDR"
    ]
]

# ============================================================
# EXPORT CSV
# ============================================================

final_df.to_csv(
    CSV_OUT,
    index=False
)

print(f"CSV saved: {CSV_OUT}")

# ============================================================
# EXPORT XLSX
# ============================================================

with pd.ExcelWriter(XLSX_OUT, engine="openpyxl") as writer:

    final_df.to_excel(
        writer,
        sheet_name="Table9",
        index=False
    )

print(f"XLSX saved: {XLSX_OUT}")

# ============================================================
# EXPORT LATEX
# ============================================================

latex_table = final_df.to_latex(
    index=False,
    escape=False,
    longtable=False
)

with open(LATEX_OUT, "w") as f:
    f.write(latex_table)

print(f"LaTeX saved: {LATEX_OUT}")

# ============================================================
# SUMMARY
# ============================================================

print("\n====================================================")
print("TABLE 9 GENERATED")
print("====================================================")
print(f"Rows exported: {len(final_df):,}")
print(f"CSV:   {CSV_OUT}")
print(f"XLSX:  {XLSX_OUT}")
print(f"LaTeX: {LATEX_OUT}")
print("====================================================")
