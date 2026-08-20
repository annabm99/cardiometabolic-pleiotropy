#!/usr/bin/env python3

"""
Generate Figure 3 input datasets from Table4_PleioLoci.csv

Outputs:
--------
1. Table11_PleioCounts.csv
    One row per locus × disease × trait

2. Fig3_DiseaseBurden.csv
    Unique pleiotropic loci per disease and direction

3. Fig3_MatrixCounts.csv
    Disease × Trait matrix counts for plotting

Author: Anna Basquet project
"""

import os
import re
import pandas as pd

# ============================================================
# PATHS
# ============================================================

PROJECT_DIR = os.environ.get(
    "CVP_PROJECT_DIR",
    "/path/to/cardiovascular_pleiotropies"
)

TABLE4_PATH = os.environ.get(
    "CVP_TABLE4_FOR_FIG3",
    f"{PROJECT_DIR}/FinalTables/Table4_PleioLoci.csv"
)
# Figure 3 summarizes genome-wide PleioFDR results, where HDL and LDL were
# retained. Set CVP_TABLE4_FOR_FIG3 to a no-HDL/LDL table only if intentionally
# reproducing the downstream convergence branch instead.

OUTPUT_DIR = os.environ.get(
    "CVP_FINAL_TABLES_DIR",
    f"{PROJECT_DIR}/FinalTables"
)

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# OUTPUT FILES
# ============================================================

TABLE11_OUT = os.path.join(
    OUTPUT_DIR,
    "Table11_PleioCounts.csv"
)

BURDEN_OUT = os.path.join(
    OUTPUT_DIR,
    "Fig3_DiseaseBurden.csv"
)

MATRIX_OUT = os.path.join(
    OUTPUT_DIR,
    "Fig3_MatrixCounts.csv"
)

# ============================================================
# EXPECTED VALUES
# ============================================================

VALID_DISEASES = {
    "CAD",
    "HT",
    "STR",
    "T2D"
}

VALID_TRAITS = {
    "BMI",
    "WC",
    "SBP",
    "DBP",
    "FG",
    "HDL",
    "LDL",
    "TGL"
}

DISEASE_ORDER = [
    "CAD",
    "HT",
    "STR",
    "T2D"
]

TRAIT_ORDER = [
    "BMI",
    "WC",
    "SBP",
    "DBP",
    "FG",
    "HDL",
    "LDL",
    "TGL"
]

# ============================================================
# LOAD TABLE 4
# ============================================================

print(f"\nLoading Table 4 source: {TABLE4_PATH}")

t4 = pd.read_csv(TABLE4_PATH)

print(f"Rows loaded: {len(t4):,}")

# ============================================================
# PARSE PHENOTYPE PAIRS
# ============================================================

print("\nParsing Phenotype_Pairs ...")

records = []

pair_pattern = re.compile(
    r"^(CAD|HT|STR|T2D)_d-(BMI|WC|SBP|DBP|FG|HDL|LDL|TGL)_t$"
)

for _, row in t4.iterrows():

    locus = row["MergedLocusID"]
    direction = row["Pleiotropy_Type"]

    pairs = str(row["Phenotype_Pairs"]).split("|")

    # --------------------------------------------------------
    # Deduplicate within locus
    # --------------------------------------------------------

    pairs = sorted(set(pairs))

    for pair in pairs:

        pair = pair.strip()

        match = pair_pattern.match(pair)

        ignored_pairs = []

        if match is None:
            ignored_pairs.append(pair)
            continue

        disease = match.group(1)
        trait = match.group(2)

        records.append(
            {
                "MergedLocusID": locus,
                "Disease": disease,
                "Trait": trait,
                "Pleiotropy_Type": direction
            }
        )

# ============================================================
# TABLE 11
# ============================================================

table11 = pd.DataFrame(records)

table11 = table11.drop_duplicates()

table11.to_csv(
    TABLE11_OUT,
    index=False
)

print(
    f"\nSaved Table11_PleioCounts.csv "
    f"({len(table11):,} rows)"
)

# ============================================================
# FIGURE 3A
# DISEASE BURDEN
# ============================================================

print("\nGenerating disease burden table ...")

disease_burden = (
    table11
    .groupby(
        ["Disease", "Pleiotropy_Type"]
    )["MergedLocusID"]
    .nunique()
    .unstack(fill_value=0)
)

# Ensure all diseases present

disease_burden = (
    disease_burden
    .reindex(DISEASE_ORDER)
    .fillna(0)
)

# Ensure both columns exist

for col in ["Concordant", "Discordant"]:
    if col not in disease_burden.columns:
        disease_burden[col] = 0

disease_burden["Total"] = (
    disease_burden["Concordant"]
    + disease_burden["Discordant"]
)

disease_burden = (
    disease_burden
    .reset_index()
)

disease_burden.columns.name = None

disease_burden.to_csv(
    BURDEN_OUT,
    index=False
)

print(
    f"Saved Fig3_DiseaseBurden.csv "
    f"({len(disease_burden)} rows)"
)

# ============================================================
# FIGURE 3B
# MATRIX COUNTS
# ============================================================

print("\nGenerating matrix counts table ...")

matrix_counts = (
    table11
    .groupby(
        ["Disease", "Trait", "Pleiotropy_Type"]
    )["MergedLocusID"]
    .nunique()
    .unstack(fill_value=0)
)

# Ensure columns exist

for col in ["Concordant", "Discordant"]:
    if col not in matrix_counts.columns:
        matrix_counts[col] = 0

# ------------------------------------------------------------
# Create full 4 × 8 matrix
# ------------------------------------------------------------

full_grid = pd.MultiIndex.from_product(
    [DISEASE_ORDER, TRAIT_ORDER],
    names=["Disease", "Trait"]
)

matrix_counts = (
    matrix_counts
    .reindex(full_grid, fill_value=0)
    .reset_index()
)

matrix_counts.columns.name = None

matrix_counts.to_csv(
    MATRIX_OUT,
    index=False
)

print(
    f"Saved Fig3_MatrixCounts.csv "
    f"({len(matrix_counts)} rows)"
)

# ============================================================
# QC SUMMARY
# ============================================================

print("\nIgnored pair types:")

for p in sorted(set(ignored_pairs)):
    print(f"  {p}")

print("\n" + "=" * 60)
print("QC SUMMARY")
print("=" * 60)

print(
    f"\nUnique pleiotropic loci: "
    f"{table11['MergedLocusID'].nunique():,}"
)

print(
    f"Table11 rows: "
    f"{len(table11):,}"
)

print("\nDisease burden:")

print(
    disease_burden[
        [
            "Disease",
            "Concordant",
            "Discordant",
            "Total"
        ]
    ]
)

print("\nMatrix dimensions:")

print(matrix_counts.shape)

print("\nExpected matrix size: 32 rows")

if len(matrix_counts) != 32:
    print(
        "WARNING: Matrix does not contain "
        "exactly 32 disease-trait combinations."
    )

print("\nOutputs written to:")

print(TABLE11_OUT)
print(BURDEN_OUT)
print(MATRIX_OUT)

print("\nDone.")
