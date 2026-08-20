#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# =============================================================================
# PATHS
# =============================================================================

PROJECT_DIR = os.environ.get(
    "CVP_PROJECT_DIR",
    "/path/to/cardiovascular_pleiotropies"
)

FINAL_TABLES_DIR = os.environ.get(
    "CVP_FINAL_TABLES_DIR",
    f"{PROJECT_DIR}/FinalTables"
)

INPUT_FILE = os.environ.get(
    "CVP_TABLE9_PATHWAYS",
    f"{FINAL_TABLES_DIR}/Table9_Top_Pathway_Enrichment.csv"
)

OUTPUT_MAPPING = os.environ.get(
    "CVP_TABLE10_MAPPING",
    f"{FINAL_TABLES_DIR}/Table10_PathwayThemeMapping.csv"
)

OUTPUT_SUMMARY = os.environ.get(
    "CVP_TABLE10_SUMMARY",
    f"{FINAL_TABLES_DIR}/Table10_ThemeSummary.csv"
)

# =============================================================================
# LOAD
# =============================================================================

df = pd.read_csv(INPUT_FILE)

# =============================================================================
# FILTER
# =============================================================================

df = df[
    (df["Category"] == "GWAScatalog") &
    (df["N_Diseases"] >= 3)
].copy()

print(f"Rows after filtering: {len(df)}")

# =============================================================================
# THEME ASSIGNMENT
# =============================================================================

def assign_theme(pathway):

    p = pathway.lower()

    # -------------------------------------------------------------------------
    # Cardiovascular Disease
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "coronary artery disease",
        "coronary heart disease",
        "myocardial infarction",
        "ischemic stroke",
        "cardioembolic stroke",
        "lacunar stroke",
        "intracranial aneurysm",
        "intracerebral hemorrhage",
        "stroke",
        "aortic root",
        "electrocardiogram",
        "intracranial",
        "haemorrhoidal disease"
    ]):
        return "Cardiovascular Disease"

    # -------------------------------------------------------------------------
    # Blood Pressure & Hemodynamics
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "blood pressure",
        "systolic blood pressure",
        "diastolic blood pressure",
        "mean arterial pressure",
        "pulse pressure",
        "hypertension",
        "resting heart rate"
    ]):
        return "Blood Pressure & Hemodynamics"

    # -------------------------------------------------------------------------
    # Adiposity & Body Composition
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "body mass index",
        "waist",
        "hip circumference",
        "hip index",
        "weight",
        "body fat distribution",
        "appendicular lean mass",
        "body shape",
        "adult body size",
        "mineral density"
    ]):
        return "Adiposity & Body Composition"

    # -------------------------------------------------------------------------
    # Lipids & Lipoproteins
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "hdl",
        "ldl",
        "cholesterol",
        "triglyceride",
        "apolipoprotein"
    ]):
        return "Lipids & Lipoproteins"

    # -------------------------------------------------------------------------
    # Glycaemic & Metabolic
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "type 2 diabetes",
        "fasting insulin",
        "glycated hemoglobin",
        "metabolic syndrome",
        "leptin",
        "serum metabolite"
    ]):
        return "Glycaemic & Metabolic"

    # -------------------------------------------------------------------------
    # Hematological
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "platelet",
        "reticulocyte",
        "hematocrit",
        "hemoglobin",
        "white blood cell",
        "myeloid",
        "neutrophil",
        "monocyte",
        "eosinophil",
        "red blood cell",
        "corpuscular",
        "hematological",
        "red cell"
    ]):
        return "Hematological"

    # -------------------------------------------------------------------------
    # Neurological / Psychiatric / Cognitive
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "brain morphology",
        "cortical",
        "subcortical",
        "schizophrenia",
        "autism",
        "anorexia",
        "mood instability",
        "intelligence",
        "cognitive",
        "alzheimer",
        "white matter",
        "hearing difficulty",
        "cbt",
        "sleep duration",
        "multisite chronic pain"
    ]):
        return "Neurological & Psychiatric"

    # -------------------------------------------------------------------------
    # Immune / Inflammatory
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "crohn",
        "ulcerative colitis",
        "inflammatory bowel disease",
        "asthma",
        "vitiligo",
        "c-reactive protein",
        "endometriosis"
    ]):
        return "Immune & Inflammatory"

    # -------------------------------------------------------------------------
    # Renal
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "glomerular filtration",
        "chronic kidney disease",
        "blood urea nitrogen",
        "albumin-to-creatinine",
        "urate",
        "uric acid"
    ]):
        return "Renal Function"

    # -------------------------------------------------------------------------
    # Cancer
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "cancer",
        "glioma",
        "glioblastoma",
        "carcinoma"
    ]):
        return "Cancer"

    # -------------------------------------------------------------------------
    # Ophthalmological
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "glaucoma",
        "intraocular",
        "refractive error"
    ]):
        return "Ophthalmological"

    # -------------------------------------------------------------------------
    # Longevity & Aging
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "parental longevity",
        "parental lifespan",
        "attained age",
        "age at death"
    ]):
        return "Longevity & Aging"

    # -------------------------------------------------------------------------
    # Lifestyle & Behaviour
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "smoking",
        "coffee",
        "fruit consumption",
        "morning person",
        "gym",
        "religious group",
        "sexual intercourse",
        "household income"
    ]):
        return "Lifestyle & Behaviour"

    # -------------------------------------------------------------------------
    # Medication Use
    # -------------------------------------------------------------------------

    if "medication use" in p:
        return "Medication Use"

    # -------------------------------------------------------------------------
    # Liver / Biomarkers
    # -------------------------------------------------------------------------

    if any(x in p for x in [
        "aminotransferase",
        "alkaline phosphatase",
        "non-albumin protein",
        "spleen volume",
        "arsenic metabolism",
        "alanine transaminase",
        "liver enzyme"
    ]):
        return "Liver & Biomarkers"

    return "UNCLASSIFIED"

# =============================================================================
# APPLY
# =============================================================================

df["Theme"] = df["Pathway"].apply(assign_theme)

# =============================================================================
# SAVE PATHWAY MAPPING
# =============================================================================

mapping = (
    df[
        [
            "Pathway",
            "Theme",
            "Pleiotropy_Direction",
            "N_Diseases",
            "N_Phenotype_Pairs"
        ]
    ]
    .drop_duplicates()
)

mapping.to_csv(
    OUTPUT_MAPPING,
    index=False
)

# =============================================================================
# THEME SUMMARY
# =============================================================================

summary = (
    df.groupby(
        ["Theme", "Pleiotropy_Direction"],
        as_index=False
    )
    .agg(
        Total_Pairs=("N_Phenotype_Pairs", "sum"),
        N_Unique_Pathways=("Pathway", "nunique")
    )
)

summary["Mean_Pairs_Per_Pathway"] = (
    summary["Total_Pairs"] /
    summary["N_Unique_Pathways"]
)

summary["Relative_Contribution"] = (
    summary.groupby("Pleiotropy_Direction")["Total_Pairs"]
    .transform(lambda x: 100 * x / x.sum())
)

# =============================================================================
# POSITIVE VS NEGATIVE ENRICHMENT TEST
# =============================================================================

wide = (
    summary.pivot(
        index="Theme",
        columns="Pleiotropy_Direction",
        values="Total_Pairs"
    )
    .fillna(0)
)

wide.columns.name = None

positive_total = wide["positive"].sum()
negative_total = wide["negative"].sum()

results = []

for theme, row in wide.iterrows():

    pos_theme = row["positive"]
    neg_theme = row["negative"]

    pos_other = positive_total - pos_theme
    neg_other = negative_total - neg_theme

    contingency = [
        [pos_theme, neg_theme],
        [pos_other, neg_other]
    ]

    odds_ratio, pvalue = fisher_exact(contingency)

    fold_change = (
        (pos_theme / positive_total) /
        (neg_theme / negative_total)
    )

    results.append({
        "Theme": theme,
        "Positive_Pairs": pos_theme,
        "Negative_Pairs": neg_theme,
        "Fold_Change": fold_change,
        "Odds_Ratio": odds_ratio,
        "P_Value": pvalue
    })

stats_df = pd.DataFrame(results)

stats_df["FDR"] = multipletests(
    stats_df["P_Value"],
    method="fdr_bh"
)[1]

def significance_star(fdr):

    if fdr < 0.001:
        return "***"
    elif fdr < 0.01:
        return "**"
    elif fdr < 0.05:
        return "*"
    else:
        return ""

stats_df["Significance"] = (
    stats_df["FDR"]
    .apply(significance_star)
)

# =============================================================================
# MERGE STATISTICS INTO SUMMARY TABLE
# =============================================================================

summary = summary.merge(
    stats_df[
        [
            "Theme",
            "Fold_Change",
            "Odds_Ratio",
            "P_Value",
            "FDR",
            "Significance"
        ]
    ],
    on="Theme",
    how="left"
)

# =============================================================================
# POSITIVE vs NEGATIVE CONTRAST
# =============================================================================

positive_map = (
    summary.loc[
        summary["Pleiotropy_Direction"] == "positive",
        ["Theme", "Relative_Contribution"]
    ]
    .set_index("Theme")
    ["Relative_Contribution"]
)

negative_map = (
    summary.loc[
        summary["Pleiotropy_Direction"] == "negative",
        ["Theme", "Relative_Contribution"]
    ]
    .set_index("Theme")
    ["Relative_Contribution"]
)

summary["Positive_Percent"] = (
    summary["Theme"]
    .map(positive_map)
)

summary["Negative_Percent"] = (
    summary["Theme"]
    .map(negative_map)
)

summary["Fold_Change"] = (
    summary["Positive_Percent"] /
    summary["Negative_Percent"]
)

summary = summary.sort_values(
    "Total_Pairs",
    ascending=False
)

summary.to_csv(
    OUTPUT_SUMMARY,
    index=False
)
# ---------------------------------------------------------
# Average recurrence per pathway
# ---------------------------------------------------------

summary["Mean_Pairs_Per_Pathway"] = (
    summary["Total_Pairs"] /
    summary["N_Unique_Pathways"]
)

summary = summary.sort_values(
    "Total_Pairs",
    ascending=False
)

# =============================================================================
# QC
# =============================================================================

print("\nTheme summary:\n")

print(
    summary.sort_values(
        "Total_Pairs",
        ascending=False
    ).to_string(index=False)
)

print("\nUnclassified pathways:")

unclassified = (
    df.loc[
        df["Theme"] == "UNCLASSIFIED",
        "Pathway"
    ]
    .drop_duplicates()
)

print(unclassified.to_string(index=False))

print(f"\nSaved: {OUTPUT_MAPPING}")
print(f"Saved: {OUTPUT_SUMMARY}")
