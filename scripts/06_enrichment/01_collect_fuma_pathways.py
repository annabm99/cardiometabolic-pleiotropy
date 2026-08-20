#!/usr/bin/env python3

"""
01_collect_fuma_pathways.py

Purpose
-------
Recursively scan FUMA GENE2FUNC output directories and collect all
gene-set enrichment results (GS.txt) into a unified master table.

This script:
- Detects positive and negative pleiotropy analyses
- Extracts phenotype pair metadata
- Reads GS.txt files
- Standardizes columns
- Saves a combined intermediate dataset

Expected structure
------------------
4-FUMA/5-FUMAOut/a-Pleiotropic/
    ├── CAD_d-BMI_t/
    │   ├── CAD_d-BMI_t-PosFunc/GS.txt
    │   └── CAD_d-BMI_t-NegFunc/GS.txt
    ├── HT_d-FG_t/
    │   ├── ...
    │   └── ...
    └── ...

Output
------
5-EnrichmentAnalysis/intermediate/all_fuma_pathways.tsv

Usage
-----
python 01_collect_fuma_pathways.py FUMA_ROOT OUTDIR

"""

import pandas as pd
from pathlib import Path
import re
import sys

# ============================================================
# PATHS
# ============================================================

# EDIT THESE IF NEEDED
FUMA_ROOT = Path(sys.argv[1].rstrip())
OUTDIR = Path(sys.argv[2].rstrip())

OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "all_fuma_pathways.tsv"

# Final downstream convergence/FUMA/enrichment analyses exclude cholesterol
# traits, although HDL and LDL were retained in the genome-wide PleioFDR stage.
EXCLUDED_TRAITS = {"HDL_t", "LDL_t"}

# ============================================================
# FUNCTIONS
# ============================================================

def infer_direction(folder_name):
    """
    Determine pleiotropy direction from folder name.
    """
    if "PosFunc" in folder_name:
        return "positive"
    elif "NegFunc" in folder_name:
        return "negative"
    else:
        return "unknown"


def infer_pheno_pair(folder_name):
    """
    Extract phenotype pair name from folder name.

    Example:
    CAD_d-BMI_t-PosFunc -> CAD_d-BMI_t
    """
    folder_name = folder_name.replace("-PosFunc", "")
    folder_name = folder_name.replace("-NegFunc", "")
    return folder_name


# ============================================================
# MAIN
# ============================================================

all_dfs = []

gs_files = sorted(FUMA_ROOT.rglob("GS.txt"))

if len(gs_files) == 0:
    sys.exit(f"ERROR: No GS.txt files found under:\n{FUMA_ROOT}")

print(f"Found {len(gs_files)} GS.txt files")

for gs_file in gs_files:

    try:
        parent_folder = gs_file.parent.name

        direction = infer_direction(parent_folder)
        phenotype_pair = infer_pheno_pair(parent_folder)

        if any(trait in phenotype_pair for trait in EXCLUDED_TRAITS):
            print(f"Skipping excluded cholesterol pair: {phenotype_pair}")
            continue

        print(f"Processing: {phenotype_pair} ({direction})")

        # ----------------------------------------------------
        # READ FILE
        # ----------------------------------------------------

        df = pd.read_csv(
            gs_file,
            sep="\t",
            engine="python"
        )

        # Explicit FUMA GS.txt columns
        df.columns = [
            "category",
            "pathway_name",
            "geneset_size",
            "n_overlap",
            "p_value",
            "fdr",
            "overlapping_genes",
            "link"
        ]


        # ----------------------------------------------------
        # ADD METADATA
        # ----------------------------------------------------

        df["phenotype_pair"] = phenotype_pair
        df["pleiotropy_direction"] = direction

        # ----------------------------------------------------
        # KEEP FINAL COLUMNS
        # ----------------------------------------------------

        df = df[
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

        all_dfs.append(df)

    except Exception as e:
        print(f"WARNING: Failed to process {gs_file}")
        print(e)

# ============================================================
# CONCATENATE
# ============================================================

if len(all_dfs) == 0:
    sys.exit("ERROR: No valid GS.txt files processed")

master_df = pd.concat(all_dfs, ignore_index=True)

# ============================================================
# SAVE
# ============================================================

master_df.to_csv(OUTFILE, sep="\t", index=False)

print("\n====================================================")
print("DONE")
print("====================================================")
print(f"Rows collected: {len(master_df):,}")
print(f"Output saved to:\n{OUTFILE}")
print("====================================================")
