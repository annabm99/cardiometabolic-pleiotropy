#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

# =============================================================================
# FILE PATHS
# =============================================================================

PROJECT_DIR = Path(
    os.environ.get("CVP_PROJECT_DIR", "/path/to/cardiovascular_pleiotropies")
)

TABLE1 = PROJECT_DIR / "FinalTables/Table1_dataset_overview.csv"
TABLE2 = PROJECT_DIR / "FinalTables/Table2_HeritabilityEstimates.csv"
TABLE3 = PROJECT_DIR / "FinalTables/Table3_GeneticCorrelations.csv"

OUTDIR = PROJECT_DIR / "Figures"
OUTDIR.mkdir(exist_ok=True)

# =============================================================================
# SETTINGS
# =============================================================================

sns.set_style("white")

POSITIVE_COLOR = "#E41A1C"   # red
NEGATIVE_COLOR = "#0072B2"   # blue

TEXT = "#333333"

TITLE_SIZE = 18
LABEL_SIZE = 14
TICK_SIZE = 14

# -----------------------------------------------------------------------------
# DOMAIN COLORS
# -----------------------------------------------------------------------------

DOMAIN_COLORS = {
    "Cardiometabolic diseases": "#1B9E77",
    "Anthropometric traits": "#D95F02",
    "Blood pressure traits": "#7570B3",
    "Metabolic biomarkers": "#E7298A",
}

# =============================================================================
# PHENOTYPE DEFINITIONS
# =============================================================================

ordered_codes = [
    "CAD_d",
    "HT_d",
    "STR_d",
    "T2D_d",
    "BMI_t",
    "WC_t",
    "DBP_t",
    "SBP_t",
    "FG_t",
    "HDL_t",
    "LDL_t",
    "TGL_t"
]

code_to_name = {
    "CAD_d": "Coronary Artery Disease",
    "HT_d": "Hypertension",
    "STR_d": "Stroke",
    "T2D_d": "Type 2 Diabetes",
    "BMI_t": "Body Mass Index",
    "WC_t": "Waist Circumference",
    "DBP_t": "Diastolic Blood Pressure",
    "SBP_t": "Systolic Blood Pressure",
    "FG_t": "Fasting Glucose",
    "HDL_t": "HDL Cholesterol",
    "LDL_t": "LDL Cholesterol",
    "TGL_t": "Triglycerides"
}

domain_map = {
    "CAD_d": "Clinical disease",
    "HT_d": "Clinical disease",
    "STR_d": "Clinical disease",
    "T2D_d": "Clinical disease",

    "BMI_t": "Anthropometric Trait",
    "WC_t": "Anthropometric Trait",

    "DBP_t": "Blood pressure measure",
    "SBP_t": "Blood pressure measure",

    "FG_t": "Metabolic marker",
    "HDL_t": "Metabolic marker",
    "LDL_t": "Metabolic marker",
    "TGL_t": "Metabolic marker"
}

domain_order = [
    "Clinical disease",
    "Anthropometric Trait",
    "Blood pressure measure",
    "Metabolic marker"
]

phenotype_type = {
    "CAD_d": "Disease",
    "HT_d": "Disease",
    "STR_d": "Disease",
    "T2D_d": "Disease",

    "BMI_t": "Trait",
    "WC_t": "Trait",
    "DBP_t": "Trait",
    "SBP_t": "Trait",
    "FG_t": "Trait",
    "HDL_t": "Trait",
    "LDL_t": "Trait",
    "TGL_t": "Trait"
}

# =============================================================================
# PANEL A
# HERITABILITY DOT PLOT
# =============================================================================

h2 = pd.read_csv(TABLE2)

h2.columns = h2.columns.str.strip()

h2 = h2.set_index("Phenotype")
h2 = h2.loc[ordered_codes].reset_index()


fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(len(h2))

ax.errorbar(
    x,
    h2["h2"],
    yerr=h2["SE"],
    fmt="none",
    ecolor="black",
    capsize=3,
    linewidth=1.2,
    zorder=1
)

# diseases
disease_mask = h2["Phenotype"].map(
    lambda x: phenotype_type[x] == "Disease"
)

# traits
trait_mask = ~disease_mask

ax.scatter(
    x[disease_mask],
    h2.loc[disease_mask, "h2"],
    s=150,
    marker="^",
    color="#333333",
    edgecolor="black",
    linewidth=0.8,
    zorder=2
)

ax.scatter(
    x[trait_mask],
    h2.loc[trait_mask, "h2"],
    s=120,
    marker="o",
    color="#333333",
    edgecolor="black",
    linewidth=0.8,
    zorder=2
)

from matplotlib.lines import Line2D

legend_elements = [
    Line2D(
        [0], [0],
        marker="^",
        color="w",
        markerfacecolor="#333333",
        markeredgecolor="black",
        markersize=10,
        label="Disease"
    ),
    Line2D(
        [0], [0],
        marker="o",
        color="w",
        markerfacecolor="#333333",
        markeredgecolor="black",
        markersize=10,
        label="Quantitative trait"
    )
]

ax.legend(
    handles=legend_elements,
    frameon=True,
    loc="upper left",
    fontsize=14
)

ax.set_xticks(x)

ax.set_xticklabels(
    [code_to_name[p] for p in h2["Phenotype"]],
    rotation=45,
    ha="right",
    fontsize=TICK_SIZE
)

ax.set_ylabel(
    r"SNP heritability ($h^2$)",
    fontsize=LABEL_SIZE
)

ax.set_xlabel("")
ax.set_title(
    "Heritability Estimates",
    fontsize=TITLE_SIZE,
    pad=15
)

sns.despine()

plt.tight_layout()

plt.savefig(
    OUTDIR / "Figure1A_Heritability.svg",
    bbox_inches="tight"
)

plt.savefig(
    OUTDIR / "Figure1A_Heritability.png",
    dpi=600,
    bbox_inches="tight"
)

plt.close()

# =============================================================================
# PANEL B
# GENETIC CORRELATION HEATMAP
# =============================================================================

rg = pd.read_csv(TABLE3)

rg.columns = rg.columns.str.strip()

heatmap_order = [
    "Coronary artery disease",
    "Hypertension",
    "Stroke",
    "Type 2 diabetes",
    "Body mass index",
    "Waist circumference",
    "Diastolic blood pressure",
    "Systolic blood pressure",
    "Fasting glucose",
    "HDL cholesterol",
    "LDL cholesterol",
    "Triglycerides"
]

# =============================================================================
# BUILD MATRICES
# =============================================================================

rg_matrix = pd.DataFrame(
    np.nan,
    index=heatmap_order,
    columns=heatmap_order
)

sig_matrix = pd.DataFrame(
    np.nan,
    index=heatmap_order,
    columns=heatmap_order
)

for _, row in rg.iterrows():

    p1 = row["Phenotype 1"]
    p2 = row["Phenotype 2"]

    rg_matrix.loc[p1, p2] = row["Genetic correlation (rg)"]
    rg_matrix.loc[p2, p1] = row["Genetic correlation (rg)"]

    sig_matrix.loc[p1, p2] = row["FDR-adjusted p-value"]
    sig_matrix.loc[p2, p1] = row["FDR-adjusted p-value"]

for p in heatmap_order:
    rg_matrix.loc[p, p] = 1
    sig_matrix.loc[p, p] = 0

# labels
labels = []

for pheno in heatmap_order:

    if pheno in [
        "Coronary artery disease",
        "Hypertension",
        "Stroke",
        "Type 2 diabetes"
    ]:
        labels.append(f"▲ {pheno}")
    else:
        labels.append(f"● {pheno}")

# significance stars
annot = np.empty(rg_matrix.shape, dtype=object)

for i in range(rg_matrix.shape[0]):
    for j in range(rg_matrix.shape[1]):

        annot[i, j] = (
            ""
            if i == j
            else "*" if sig_matrix.iloc[i, j] < 0.05
            else ""
        )

# -----------------------------------------------------------------------------
# HEATMAP
# -----------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(11, 9))

cmap = sns.diverging_palette(
    240,
    10,
    as_cmap=True
)

sns.heatmap(
    rg_matrix,
    cmap=cmap,
    vmin=-1,
    vmax=1,
    center=0,
    square=True,
    linewidths=0.5,
    linecolor="white",
    annot=annot,
    fmt="",
    cbar_kws={"label": r"$r_g$"},
    ax=ax
)

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=14)
cbar.set_label(r"$r_g$", fontsize=16)

ax.set_xticklabels(
    labels,
    rotation=45,
    ha="right",
    fontsize=TICK_SIZE
)

ax.set_yticklabels(
    labels,
    rotation=0,
    fontsize=TICK_SIZE
)

ax.set_title(
    "Genetic Correlation Matrix",
    fontsize=TITLE_SIZE,
    pad=15
)

# -----------------------------------------------------------------------------
# DOMAIN SEPARATORS
# -----------------------------------------------------------------------------

domain_sizes = [4, 2, 2, 4]

running = 0

for n in domain_sizes[:-1]:

    running += n

    ax.axhline(
        running,
        color="black",
        linewidth=1,
        alpha=0.6
    )

    ax.axvline(
        running,
        color="black",
        linewidth=1,
        alpha=0.6
    )

plt.tight_layout()

plt.savefig(
    OUTDIR / "Figure1B_GeneticCorrelationHeatmap.svg",
    bbox_inches="tight"
)

plt.savefig(
    OUTDIR / "Figure1B_GeneticCorrelationHeatmap.png",
    dpi=600,
    bbox_inches="tight"
)

plt.close()

print("Figure 1 exported successfully.")
