#!/usr/bin/env python3

"""
03_collapse_fuma_pathways.py

Purpose
-------
Collapse and curate significant FUMA pathway enrichments
to generate a cleaner, publication-oriented enrichment dataset.

This script:
- Removes exact duplicates
- Prioritizes pathway databases
- Removes less informative categories
- Collapses duplicated pathways across databases
- Keeps the best (lowest FDR) representative pathway

Input
-----
all_fuma_pathways_fdr05.tsv

Output
------
all_fuma_pathways_nonredundant.tsv

Usage
-----
python 03_collapse_fuma_pathways.py INPUT_FILE OUTPUT_DIR

Example
-------
python 03_collapse_fuma_pathways.py \
    /path/to/all_fuma_pathways_fdr05.tsv \
    /path/to/5-Enrichment/3-Collapsed
"""

import pandas as pd
from pathlib import Path
import sys
import re

# ============================================================
# INPUTS
# ============================================================

INPUT_FILE = Path(sys.argv[1].rstrip())
OUTDIR = Path(sys.argv[2].rstrip())

OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "all_fuma_pathways_nonredundant.tsv"

# ============================================================
# LOAD
# ============================================================

print("Loading significant enrichment dataset...")

df = pd.read_csv(INPUT_FILE, sep="\t")

print(f"Initial rows: {len(df):,}")

# ============================================================
# REMOVE EXACT DUPLICATES
# ============================================================

before = len(df)

df = df.drop_duplicates()

after = len(df)

print(f"Removed exact duplicates: {before - after:,}")

# ============================================================
# REMOVE LESS INFORMATIVE CATEGORIES
# ============================================================

REMOVE_CATEGORIES = [
    "Positional_gene_sets",
    "TF_targets"
]

before = len(df)

df = df[~df["category"].isin(REMOVE_CATEGORIES)].copy()

after = len(df)

print(f"Removed low-interpretability categories: {before - after:,}")

# ============================================================
# DATABASE PRIORITY
# ============================================================

DATABASE_PRIORITY = {
    "KEGG": 1,
    "Reactome": 2,
    "GO_bp": 3,
    "GO_mf": 4,
    "GO_cc": 5,
    "Canonical_Pathways": 6,
    "GWAScatalog": 7
}

df["db_priority"] = df["category"].map(DATABASE_PRIORITY)

# Unknown categories get lowest priority
df["db_priority"] = df["db_priority"].fillna(999)

# ============================================================
# CLEAN PATHWAY NAMES
# ============================================================

def clean_pathway_name(name):

    name = str(name)

    # Remove common prefixes
    prefixes = [
        "KEGG_",
        "REACTOME_",
        "WP_",
        "GOBP_",
        "GO_BP_",
        "HALLMARK_"
    ]

    for prefix in prefixes:
        if name.startswith(prefix):
            name = name.replace(prefix, "", 1)

    # Normalize separators
    name = name.replace("_", " ")

    # Uppercase consistency
    name = re.sub(r"\s+", " ", name)

    return name.strip().lower()

df["clean_pathway"] = df["pathway_name"].apply(clean_pathway_name)

# ============================================================
# SORT FOR BEST REPRESENTATIVE SELECTION
# ============================================================

df = df.sort_values(
    by=[
        "phenotype_pair",
        "pleiotropy_direction",
        "clean_pathway",
        "db_priority",
        "fdr",
        "p_value"
    ]
)

# ============================================================
# COLLAPSE REDUNDANT PATHWAYS
# ============================================================

before = len(df)

df_nonredundant = df.drop_duplicates(
    subset=[
        "phenotype_pair",
        "pleiotropy_direction",
        "clean_pathway"
    ],
    keep="first"
)

after = len(df_nonredundant)

print(f"Collapsed redundant pathways: {before - after:,}")

# ============================================================
# FINAL SORT
# ============================================================

df_nonredundant = df_nonredundant.sort_values(
    by=[
        "phenotype_pair",
        "pleiotropy_direction",
        "fdr",
        "p_value"
    ]
)

# ============================================================
# FINAL COLUMNS
# ============================================================

df_nonredundant = df_nonredundant[
    [
        "phenotype_pair",
        "pleiotropy_direction",
        "pathway_name",
        "category",
        "n_overlap",
        "p_value",
        "fdr"
    ]
]

# ============================================================
# SAVE
# ============================================================

df_nonredundant.to_csv(
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
print(f"Final nonredundant rows: {len(df_nonredundant):,}")
print(f"Output saved to:")
print(OUTFILE)
print("====================================================")
