#!/usr/bin/env python3

"""
Generate Supplementary Table S3: pleiotropic loci.

Primary analysis:
- 4 diseases × 6 cardiometabolic traits
- HDL and LDL excluded from the primary convergence definition
- HDL/LDL overlap retained in separate sensitivity columns
- Dis-Dis loci excluded upstream

Input:
    FinalTables/Table4_PleioLoci_noHDL_LDL.csv

Output:
    Supplementary_S3_PleiotropicLoci.csv
"""

import os
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_DIR = Path(
    os.environ.get(
        "CVP_PROJECT_DIR",
        "/path/to/cardiovascular_pleiotropies"
    )
)

FINAL_TABLES_DIR = Path(
    os.environ.get(
        "CVP_FINAL_TABLES_DIR",
        str(PROJECT_DIR / "FinalTables")
    )
)

OUT_DIR = Path(
    os.environ.get(
        "CVP_SUPPLEMENT_TABLES_DIR",
        str(FINAL_TABLES_DIR / "Supplement")
    )
)

OUT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FILE = FINAL_TABLES_DIR / "Table4_PleioLoci_noHDL_LDL.csv"
OUTPUT_FILE = OUT_DIR / "Supplementary_S3_PleiotropicLoci.csv"


# ---------------------------------------------------------------------
# Load canonical convergence table
# ---------------------------------------------------------------------

df = pd.read_csv(INPUT_FILE)

required = [
    "MergedLocusID",
    "Chromosome",
    "Merged_Start",
    "Merged_End",
    "Representative_Lead_SNP",
    "Best_P",
    "Pleiotropy_Type",
    "Diseases_Involved",
    "Number_of_Diseases",
    "Phenotype_Pairs",
    "Num_Original_Loci",
    "Cholesterol_Concordant",
    "Cholesterol_Discordant",
    "Mapped_Genes",
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise ValueError(
        "Missing required input columns: "
        + ", ".join(missing)
    )


# ---------------------------------------------------------------------
# Construct supplementary table
# ---------------------------------------------------------------------

s3 = pd.DataFrame(
    {
        "Locus ID": df["MergedLocusID"],
        "Chr": df["Chromosome"],
        "Start (bp)": df["Merged_Start"],
        "End (bp)": df["Merged_End"],
        "Lead SNP": df["Representative_Lead_SNP"],
        "Best P-Value": df["Best_P"],
        "Pleiotropy Type": df["Pleiotropy_Type"],
        "Diseases Involved": df["Diseases_Involved"],
        "N Diseases": df["Number_of_Diseases"],

        "Concordant Pairs": df["Phenotype_Pairs"].where(
            df["Pleiotropy_Type"].eq("Concordant"),
            ""
        ),

        "Discordant Pairs": df["Phenotype_Pairs"].where(
            df["Pleiotropy_Type"].eq("Discordant"),
            ""
        ),

        "Cholesterol Concordant Pairs": df["Cholesterol_Concordant"],
        "Cholesterol Discordant Pairs": df["Cholesterol_Discordant"],

        "N Original Loci": df["Num_Original_Loci"],
        "Mapped Genes": df["Mapped_Genes"],
    }
)


# ---------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------

print("=" * 60)
print("SUPPLEMENTARY TABLE S3 QC")
print("=" * 60)

print(f"Input: {INPUT_FILE}")
print(f"Rows: {len(s3):,}")

# Expected corrected primary convergence result
assert len(s3) == 1542, (
    f"Expected 1,542 loci, found {len(s3)}"
)

assert s3["Locus ID"].notna().all(), (
    "Missing locus IDs detected"
)

assert s3["Locus ID"].is_unique, (
    "Duplicated locus IDs detected"
)

direction_counts = (
    s3["Pleiotropy Type"]
    .value_counts(dropna=False)
    .to_dict()
)

print("Pleiotropy types:", direction_counts)

assert direction_counts.get("Concordant", 0) == 1279, (
    "Expected 1,279 concordant loci"
)

assert direction_counts.get("Discordant", 0) == 263, (
    "Expected 263 discordant loci"
)

# No Dis-Dis contamination
text = s3.astype(str)

assert not text.apply(
    lambda col: col.str.contains("Dis-Dis", regex=False).any()
).any(), "Dis-Dis contamination detected"

# HDL/LDL must not occur in the primary pair columns
primary_pairs = (
    s3["Concordant Pairs"].fillna("").astype(str)
    + ";"
    + s3["Discordant Pairs"].fillna("").astype(str)
)

assert not primary_pairs.str.contains(
    r"(?:^|[-_;|])HDL(?:$|[-_;|])",
    regex=True
).any(), "HDL detected in primary phenotype pairs"

assert not primary_pairs.str.contains(
    r"(?:^|[-_;|])LDL(?:$|[-_;|])",
    regex=True
).any(), "LDL detected in primary phenotype pairs"

# Direction-specific pair columns should be mutually exclusive
assert (
    s3.loc[
        s3["Pleiotropy Type"].eq("Concordant"),
        "Discordant Pairs"
    ]
    .fillna("")
    .eq("")
    .all()
), "Discordant pairs populated for concordant loci"

assert (
    s3.loc[
        s3["Pleiotropy Type"].eq("Discordant"),
        "Concordant Pairs"
    ]
    .fillna("")
    .eq("")
    .all()
), "Concordant pairs populated for discordant loci"


# ---------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------

s3.to_csv(
    OUTPUT_FILE,
    index=False
)

print()
print("QC PASSED")
print(f"Output: {OUTPUT_FILE}")
print("=" * 60)