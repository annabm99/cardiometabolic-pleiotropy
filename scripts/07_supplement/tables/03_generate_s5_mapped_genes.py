#!/usr/bin/env python3

"""
Generate Supplementary Table S5: all FUMA-mapped pleiotropic genes.

Input:
    FinalTables/Table5_GeneAnnotationSummary_noHDL_LDL.csv

This table uses the corrected Ensembl-level gene aggregation from the
primary six-trait convergence analysis. HDL/LDL evidence is retained in
separate sensitivity columns.

Important:
- one row per Ensembl ID
- Dis-Dis is excluded upstream
- repeated symbols such as Y_RNA and snoU13 remain separate Ensembl genes
- source FUMA/Ensembl nomenclature is preserved
"""

import os
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(
    os.environ.get(
        "CVP_PROJECT_DIR",
        "/path/to/cardiovascular_pleiotropies",
    )
)

INPUT = (
    PROJECT_DIR
    / "FinalTables"
    / "Table5_GeneAnnotationSummary_noHDL_LDL.csv"
)

OUT_DIR = Path(
    os.environ.get(
        "CVP_SUPPLEMENT_DIR",
        PROJECT_DIR / "FinalTables",
    )
)

OUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT = OUT_DIR / "Supplementary_S5_MappedGenes.csv"


###############################################################################
# LOAD
###############################################################################

df = pd.read_csv(INPUT)

required = {
    "symbol",
    "Ensembl_ID",
    "Gene_Type",
    "Chromosome",
    "Gene_Start",
    "Gene_End",
    "N_Observations",
    "N_Merged_Loci",
    "MergedLocusIDs",
    "Diseases",
    "N_Diseases",
    "N_Phenotype_Pairs",
    "Concordant_Pairs",
    "Discordant_Pairs",
    "Cholesterol_Concordant_Pairs",
    "Cholesterol_Discordant_Pairs",
    "pLI",
    "ncRVIS",
}

missing = required - set(df.columns)

if missing:
    raise RuntimeError(
        "Missing required columns: "
        + ", ".join(sorted(missing))
    )


###############################################################################
# QC SOURCE
###############################################################################

if df["Ensembl_ID"].duplicated().any():
    dup = df.loc[
        df["Ensembl_ID"].duplicated(False),
        ["symbol", "Ensembl_ID"],
    ]

    raise RuntimeError(
        "Duplicate Ensembl IDs detected:\n"
        + dup.to_string(index=False)
    )

if (df["N_Merged_Loci"] < 1).any():
    raise RuntimeError(
        "Mapped genes with N_Merged_Loci < 1 detected"
    )

pair_cols = [
    "Concordant_Pairs",
    "Discordant_Pairs",
    "Cholesterol_Concordant_Pairs",
    "Cholesterol_Discordant_Pairs",
]

contains_disdis = (
    df[pair_cols]
    .fillna("")
    .astype(str)
    .apply(
        lambda col:
        col.str.contains(
            "Dis-Dis",
            regex=False,
        )
    )
    .any(axis=1)
)

if contains_disdis.any():
    raise RuntimeError(
        "Dis-Dis contamination detected in S5 source"
    )


###############################################################################
# CURATED NOMENCLATURE NOTE
###############################################################################

annotation_notes = {
    "ENSG00000270316": (
        "Source FUMA symbol C10orf32-ASMT and source type "
        "protein_coding retained for reproducibility; current official "
        "symbol is BORCS7-ASMT, classified as an ncRNA read-through / "
        "NMD candidate."
    )
}

df["Annotation_Note"] = (
    df["Ensembl_ID"]
    .map(annotation_notes)
    .fillna("")
)


###############################################################################
# FORMAT S5
###############################################################################

s5 = pd.DataFrame({

    "Gene Symbol":
        df["symbol"],

    "Ensembl ID":
        df["Ensembl_ID"],

    "Gene Type":
        df["Gene_Type"],

    "Chr":
        df["Chromosome"],

    "Gene Start (bp)":
        df["Gene_Start"],

    "Gene End (bp)":
        df["Gene_End"],

    "N Observations":
        df["N_Observations"],

    "N Merged Loci":
        df["N_Merged_Loci"],

    "Merged Locus IDs":
        df["MergedLocusIDs"],

    "N Diseases":
        df["N_Diseases"],

    "Diseases":
        df["Diseases"],

    "N Phenotype Pairs":
        df["N_Phenotype_Pairs"],

    "Concordant Pairs":
        df["Concordant_Pairs"],

    "Discordant Pairs":
        df["Discordant_Pairs"],

    "Cholesterol Concordant Pairs":
        df["Cholesterol_Concordant_Pairs"],

    "Cholesterol Discordant Pairs":
        df["Cholesterol_Discordant_Pairs"],

    "pLI":
        df["pLI"],

    "ncRVIS":
        df["ncRVIS"],

    "Annotation Note":
        df["Annotation_Note"],
})


###############################################################################
# SORT / RANK
###############################################################################

s5 = (
    s5.sort_values(
        [
            "N Observations",
            "N Merged Loci",
            "N Diseases",
            "N Phenotype Pairs",
            "Gene Symbol",
            "Ensembl ID",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            True,
            True,
        ],
        kind="mergesort",
    )
    .reset_index(drop=True)
)

s5.insert(
    0,
    "Rank",
    range(1, len(s5) + 1),
)


###############################################################################
# FINAL QC
###############################################################################

assert len(s5) == 4806
assert s5["Ensembl ID"].nunique() == 4806
assert s5["Ensembl ID"].duplicated().sum() == 0
assert (s5["N Merged Loci"] >= 1).all()
assert (s5["N Diseases"] == 4).sum() == 12

rbm6 = s5[
    s5["Ensembl ID"] == "ENSG00000004534"
]

assert len(rbm6) == 1
assert int(rbm6.iloc[0]["N Diseases"]) == 3
assert rbm6.iloc[0]["Diseases"] == "CAD|HT|T2D"

borcs = s5[
    s5["Ensembl ID"] == "ENSG00000270316"
]

assert len(borcs) == 1
assert borcs.iloc[0]["Gene Symbol"] == "C10orf32-ASMT"
assert borcs.iloc[0]["Annotation Note"] != ""


###############################################################################
# SAVE
###############################################################################

s5.to_csv(
    OUTPUT,
    index=False,
)

print("=" * 70)
print("SUPPLEMENTARY TABLE S5")
print("=" * 70)

print("Rows:", len(s5))
print(
    "Unique Ensembl IDs:",
    s5["Ensembl ID"].nunique(),
)
print(
    "Four-disease genes:",
    (s5["N Diseases"] == 4).sum(),
)
print(
    "Rows containing Dis-Dis:",
    s5[
        [
            "Concordant Pairs",
            "Discordant Pairs",
            "Cholesterol Concordant Pairs",
            "Cholesterol Discordant Pairs",
        ]
    ]
    .fillna("")
    .astype(str)
    .apply(
        lambda col:
        col.str.contains(
            "Dis-Dis",
            regex=False,
        )
    )
    .any(axis=1)
    .sum(),
)

print()
print("Top 10 genes by observation count:")
print(
    s5.head(10)[
        [
            "Rank",
            "Gene Symbol",
            "Ensembl ID",
            "N Observations",
            "N Merged Loci",
            "N Diseases",
            "Diseases",
        ]
    ].to_string(index=False)
)

print()
print("RBM6:")
print(
    rbm6[
        [
            "Gene Symbol",
            "Ensembl ID",
            "N Observations",
            "N Merged Loci",
            "N Diseases",
            "Diseases",
        ]
    ].to_string(index=False)
)

print()
print("Saved:")
print(OUTPUT)
print()
print("DONE")
