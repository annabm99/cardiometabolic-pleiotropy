#!/usr/bin/env python3

"""
Generate Supplementary Table S6 components.

S6A
    Disease-specific association strength (|Z|) for concordant and
    discordant pleiotropic signals after disease-level LD pruning,
    including the sensitivity analysis excluding mixed-context variants.

S6B
    Individual mixed-discordant disease-SNP observations.

S6C
    LD-pruning retention / QC summary.

S6D
    FUMA variant-level annotation summary.

S6E
    FUMA mapped-gene-level annotation summary.

This script reformats validated analysis outputs into Supplementary
Table S6 components.
"""

import os
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(
    os.environ.get(
        "CVP_PROJECT_DIR",
        "/data/samanthafs/scratch/lab_anavarro/"
        "anbasquet/cardiovascular_pleiotropies",
    )
)

FINAL = PROJECT_DIR / "FinalTables"

OUT_DIR = Path(
    os.environ.get(
        "CVP_SUPPLEMENT_DIR",
        FINAL,
    )
)

OUT_DIR.mkdir(parents=True, exist_ok=True)

###############################################################################
# INPUTS
###############################################################################

ASSOCIATION_DIR = (
    PROJECT_DIR
    / "3-Directionality"
    / "3-AssociationStrength"
)

ASSOCIATION_FILE = (
    ASSOCIATION_DIR
    / "AssociationStrength_LDPruned_Summary.csv"
)

CASES_FILE = (
    ASSOCIATION_DIR
    / "MixedDirectionality_Cases.csv"
)

LD_QC_FILE = (
    ASSOCIATION_DIR
    / "AssociationStrength_LDPruning_QC.csv"
)

FUMA_FILE = FINAL / "Table7_FUMA_AnnotationSummary.csv"

###############################################################################
# S6A — LD-PRUNED ASSOCIATION-STRENGTH ANALYSIS
###############################################################################

assoc = pd.read_csv(
    ASSOCIATION_FILE
)

analysis_labels = {
    "Primary_LD_pruned":
        "Primary: LD-pruned",
    "Sensitivity_LD_pruned_no_mixed":
        "Sensitivity: LD-pruned, mixed-context excluded",
}

assoc["Analysis"] = (
    assoc["Analysis"]
    .map(analysis_labels)
)

if assoc["Analysis"].isna().any():
    raise RuntimeError(
        "Unexpected analysis label in association-strength results"
    )

s6a = pd.DataFrame({

    "Analysis":
        assoc["Analysis"],

    "Disease":
        assoc["Disease"],

    "N Concordant Signals":
        assoc["N_concordant"],

    "N Discordant Signals":
        assoc["N_discordant"],

    "Median |Z| Concordant":
        assoc["Concordant_median_abs_z"],

    "Median |Z| Discordant":
        assoc["Discordant_median_abs_z"],

    "Difference in Median |Z| (Discordant - Concordant)":
        assoc[
            "Delta_median_abs_z_discordant_minus_concordant"
        ],

    "Mean |Z| Concordant":
        assoc["Concordant_mean_abs_z"],

    "Mean |Z| Discordant":
        assoc["Discordant_mean_abs_z"],

    "Mann-Whitney U":
        assoc["Mann_Whitney_U"],

    "P-value":
        assoc["P_value"],

    "FDR-adjusted P-value":
        assoc["P_FDR"],

    "Bonferroni-adjusted P-value":
        assoc["P_Bonferroni"],

    "Rank-Biserial (Discordant > Concordant)":
        assoc[
            "Rank_biserial_discordant_gt_concordant"
        ],
})

s6a["Significant after FDR"] = (
    s6a["FDR-adjusted P-value"] < 0.05
).map({
    True: "Yes",
    False: "No",
})

###############################################################################
# S6B — MIXED-DISCORDANT CASES
###############################################################################

cases = pd.read_csv(CASES_FILE)

required_case_cols = {
    "SNP",
    "Disease",
    "signed_z",
    "abs_z",
    "n_positive_pairs",
    "n_negative_pairs",
    "n_total_pairs",
    "compared_traits",
}

missing = required_case_cols - set(cases.columns)

if missing:
    raise RuntimeError(
        "Missing mixed-case columns: "
        + ", ".join(sorted(missing))
    )

# A mixed case must contain at least one same-sign and one opposite-sign pair.
if not (
    (cases["n_positive_pairs"] >= 1)
    & (cases["n_negative_pairs"] >= 1)
).all():
    raise RuntimeError(
        "A purported mixed-discordant case lacks both direction classes"
    )

disease_labels = {
    "CAD_d": "Coronary Artery Disease",
    "HT_d": "Hypertension",
    "STR_d": "Stroke",
    "T2D_d": "Type 2 Diabetes",
}

s6b = pd.DataFrame({
    "SNP":
        cases["SNP"],

    "Disease":
        cases["Disease"].map(disease_labels),

    "Disease Code":
        cases["Disease"],

    "Disease Z":
        cases["signed_z"],

    "|Disease Z|":
        cases["abs_z"],

    "N Concordant Pairs":
        cases["n_positive_pairs"],

    "N Discordant Pairs":
        cases["n_negative_pairs"],

    "N Total Pairs":
        cases["n_total_pairs"],

    "Compared Traits":
        cases["compared_traits"],

    "Directionality Class":
        "Mixed discordant",
})

s6b = (
    s6b.sort_values(
        ["Disease", "SNP"],
        kind="mergesort",
    )
    .reset_index(drop=True)
)

###############################################################################
# S6C — LD-PRUNING QC
###############################################################################

ld_qc = pd.read_csv(
    LD_QC_FILE
)

s6c = pd.DataFrame({

    "Disease":
        ld_qc["Disease"],

    "Disease-SNP Observations Before LD Pruning":
        ld_qc["Original_N"],

    "Independent Signals Retained: Primary":
        ld_qc["Primary_LD_pruned_N"],

    "Independent Signals Retained: Mixed-context Excluded":
        ld_qc["No_mixed_LD_pruned_N"],
})

s6c[
    "Primary Retained (%)"
] = (
    100
    * s6c[
        "Independent Signals Retained: Primary"
    ]
    / s6c[
        "Disease-SNP Observations Before LD Pruning"
    ]
)

s6c[
    "Mixed-excluded Retained (%)"
] = (
    100
    * s6c[
        "Independent Signals Retained: Mixed-context Excluded"
    ]
    / s6c[
        "Disease-SNP Observations Before LD Pruning"
    ]
)

###############################################################################
# S6D / S6E — FUMA ANNOTATION LEVELS
###############################################################################

fuma = pd.read_csv(FUMA_FILE)

required_fuma_cols = {
    "Layer",
    "Annotation",
    "Statistic",
    "Concordant",
    "Discordant",
}

missing = required_fuma_cols - set(fuma.columns)

if missing:
    raise RuntimeError(
        "Missing FUMA summary columns: "
        + ", ".join(sorted(missing))
    )

s6d = (
    fuma.loc[
        fuma["Layer"] == "Variant-level",
        [
            "Annotation",
            "Statistic",
            "Concordant",
            "Discordant",
        ],
    ]
    .copy()
    .reset_index(drop=True)
)

s6d.insert(
    0,
    "Analysis Level",
    "Variant-level",
)

s6e = (
    fuma.loc[
        fuma["Layer"] == "Mapped-gene-level",
        [
            "Annotation",
            "Statistic",
            "Concordant",
            "Discordant",
        ],
    ]
    .copy()
    .reset_index(drop=True)
)

s6e.insert(
    0,
    "Analysis Level",
    "Mapped-gene-level",
)

###############################################################################
# QC
###############################################################################

# S6A contains four diseases for the primary analysis and the same four
# diseases for the mixed-context-excluded sensitivity analysis.
assert len(s6a) == 8

assert set(s6a["Analysis"]) == {
    "Primary: LD-pruned",
    "Sensitivity: LD-pruned, mixed-context excluded",
}

assert set(s6a["Disease"]) == {
    "Coronary Artery Disease",
    "Hypertension",
    "Stroke",
    "Type 2 Diabetes",
}

# No disease-specific comparison survives FDR correction.
assert (
    s6a["Significant after FDR"]
    .eq("No")
    .all()
)

# The established mixed-context sensitivity set contains 66
# disease-SNP observations.
assert len(s6b) == 66

assert (
    (s6b["N Concordant Pairs"] >= 1)
    & (s6b["N Discordant Pairs"] >= 1)
).all()

# S6C contains one LD-pruning QC row per disease.
assert len(s6c) == 4

assert set(s6c["Disease"]) == {
    "Coronary Artery Disease",
    "Hypertension",
    "Stroke",
    "Type 2 Diabetes",
}

# LD pruning cannot increase the number of retained observations.
assert (
    s6c["Independent Signals Retained: Primary"]
    <= s6c["Disease-SNP Observations Before LD Pruning"]
).all()

assert (
    s6c[
        "Independent Signals Retained: Mixed-context Excluded"
    ]
    <= s6c["Disease-SNP Observations Before LD Pruning"]
).all()

# FUMA annotation layers remain separated correctly.
assert set(s6d["Analysis Level"]) == {"Variant-level"}

assert set(s6e["Analysis Level"]) == {"Mapped-gene-level"}

assert {"CADD score", "RegulomeDB score"} <= set(
    s6d["Annotation"]
)

assert {"pLI", "ncRVIS"} <= set(
    s6e["Annotation"]
)

assert "pLI" not in set(s6d["Annotation"])
assert "ncRVIS" not in set(s6d["Annotation"])

assert "CADD score" not in set(s6e["Annotation"])
assert "RegulomeDB score" not in set(s6e["Annotation"])


###############################################################################
# SAVE
###############################################################################

outputs = {

    "Supplementary_S6A_AssociationStrength.csv":
        s6a,

    "Supplementary_S6B_MixedDiscordantCases.csv":
        s6b,

    "Supplementary_S6C_LDPruningQC.csv":
        s6c,

    "Supplementary_S6D_VariantAnnotations.csv":
        s6d,

    "Supplementary_S6E_GeneAnnotations.csv":
        s6e,
}

for filename, table in outputs.items():

    table.to_csv(
        OUT_DIR / filename,
        index=False,
    )


###############################################################################
# REPORT
###############################################################################

print("=" * 70)
print("SUPPLEMENTARY TABLE S6")
print("=" * 70)

print("S6A association-strength rows:", len(s6a))
print("S6B mixed-context cases:", len(s6b))
print("S6C LD-pruning QC rows:", len(s6c))
print("S6D variant-annotation rows:", len(s6d))
print("S6E gene-annotation rows:", len(s6e))

print()

print("Primary LD-pruned association-strength analysis:")

print(
    s6a.loc[
        s6a["Analysis"] == "Primary: LD-pruned",
        [
            "Disease",
            "N Concordant Signals",
            "N Discordant Signals",
            "Median |Z| Concordant",
            "Median |Z| Discordant",
            "Difference in Median |Z| (Discordant - Concordant)",
            "P-value",
            "FDR-adjusted P-value",
            "Significant after FDR",
        ],
    ].to_string(index=False)
)

print()

print("Mixed-context-excluded sensitivity:")

print(
    s6a.loc[
        s6a["Analysis"]
        == "Sensitivity: LD-pruned, mixed-context excluded",
        [
            "Disease",
            "N Concordant Signals",
            "N Discordant Signals",
            "Median |Z| Concordant",
            "Median |Z| Discordant",
            "P-value",
            "FDR-adjusted P-value",
        ],
    ].to_string(index=False)
)

print()

print("LD-pruning retention:")

print(
    s6c.to_string(
        index=False
    )
)

print()

print("Annotation separation:")

print(
    "Variant-level:",
    sorted(s6d["Annotation"].unique())
)

print(
    "Mapped-gene-level:",
    sorted(s6e["Annotation"].unique())
)

print()

print("Saved:")

for filename in outputs:
    print(OUT_DIR / filename)

print()
print("DONE")