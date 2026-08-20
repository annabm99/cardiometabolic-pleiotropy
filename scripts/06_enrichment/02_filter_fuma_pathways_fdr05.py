#!/usr/bin/env python3

"""
02_filter_fuma_pathways_fdr05.py

Purpose
-------
Filter the master FUMA pathway enrichment dataset
to retain statistically significant enrichments.

This creates the canonical downstream dataset used for:
- Table 8 (full supplementary enrichment table)
- Table 9 (summary/top pathways)

Input
-----
all_fuma_pathways.tsv

Output
------
all_fuma_pathways_fdr05.tsv

Filtering
---------
- FDR < 0.05
- Remove duplicated rows
- Remove missing values in critical columns

Usage
-----
python 02_filter_fuma_pathways_fdr05.py INPUT_FILE OUTPUT_DIR

Example
-------
python 02_filter_fuma_pathways_fdr05.py \
    /path/to/all_fuma_pathways.tsv \
    /path/to/5-Enrichment/2-Filtered
"""

import pandas as pd
from pathlib import Path
import sys

# ============================================================
# INPUTS
# ============================================================

INPUT_FILE = Path(sys.argv[1].rstrip())
OUTDIR = Path(sys.argv[2].rstrip())

OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "all_fuma_pathways_fdr05.tsv"

# ============================================================
# LOAD
# ============================================================

print("Loading master enrichment dataset...")

df = pd.read_csv(INPUT_FILE, sep="\t")

print(f"Initial rows: {len(df):,}")

# ============================================================
# BASIC CLEANING
# ============================================================

required_cols = [
    "phenotype_pair",
    "pleiotropy_direction",
    "pathway_name",
    "category",
    "n_overlap",
    "p_value",
    "fdr"
]

missing = [c for c in required_cols if c not in df.columns]

if len(missing) > 0:
    raise ValueError(f"Missing required columns: {missing}")

# Remove rows with missing critical values
df = df.dropna(subset=required_cols)

print(f"After NA removal: {len(df):,}")

# ============================================================
# ENSURE NUMERIC TYPES
# ============================================================

df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
df["fdr"] = pd.to_numeric(df["fdr"], errors="coerce")
df["n_overlap"] = pd.to_numeric(df["n_overlap"], errors="coerce")

# Remove rows becoming NA after conversion
df = df.dropna(subset=["p_value", "fdr", "n_overlap"])

# ============================================================
# FILTER SIGNIFICANT RESULTS
# ============================================================

df_sig = df[df["fdr"] < 0.05].copy()

print(f"Significant rows (FDR < 0.05): {len(df_sig):,}")

# ============================================================
# REMOVE DUPLICATES
# ============================================================

before = len(df_sig)

df_sig = df_sig.drop_duplicates()

after = len(df_sig)

print(f"Removed duplicates: {before - after:,}")

# ============================================================
# SORT
# ============================================================

df_sig = df_sig.sort_values(
    by=[
        "phenotype_pair",
        "pleiotropy_direction",
        "fdr",
        "p_value"
    ]
)

# ============================================================
# SAVE
# ============================================================

df_sig.to_csv(
    OUTFILE,
    sep="\t",
    index=False
)

# ============================================================
# SUMMARY
# ============================================================

print("\n====================================================")
print("DONE")
print("====================================================")
print(f"Final significant rows: {len(df_sig):,}")
print(f"Output saved to:")
print(OUTFILE)
print("====================================================")
