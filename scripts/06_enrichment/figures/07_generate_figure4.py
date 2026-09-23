#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

PROJECT_DIR = Path(
    os.environ.get("CVP_PROJECT_DIR", "/path/to/cardiovascular_pleiotropies")
)
FINAL_TABLES_DIR = Path(
    os.environ.get("CVP_FINAL_TABLES_DIR", str(PROJECT_DIR / "FinalTables"))
)

THEME_TABLE = Path(
    os.environ.get(
        "CVP_THEME_TABLE",
        str(
            FINAL_TABLES_DIR
            / "Table10_ThemePermutation.csv"
        )
    )
)

FIGURES_DIR = Path(
    os.environ.get("CVP_FIGURES_DIR", str(PROJECT_DIR / "Figures"))
)

ASSOCIATION_STRENGTH_DIR = Path(
    os.environ.get(
        "CVP_ASSOCIATION_STRENGTH_DIR",
        str(
            PROJECT_DIR
            / "3-Directionality"
            / "3-AssociationStrength"
        )
    )
)

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# SETTINGS
# =============================================================================

TEXT = "#333333"

TITLE_SIZE = 24
LABEL_SIZE = 20
TICK_SIZE = 20

sns.set_style("white")

# =============================================================================
# LOAD DATA
# =============================================================================

df = pd.read_csv(FINAL_TABLES_DIR / "Table7_FUMA_AnnotationSummary.csv")

# Extract percentages from:
# "864 (46.7%)"

annotation_rows = [
    "UTR3",
    "UTR5",
    "downstream",
    "exonic",
    "intergenic",
    "intronic",
    "ncRNA_exonic",
    "ncRNA_intronic",
    "upstream"
]

df = df[
    (df["Layer"] == "Variant-level")
    & (df["Statistic"] == "Count (% of unique lead SNPs)")
    & (df["Annotation"].isin(annotation_rows))
].copy()

for col in ["Concordant", "Discordant"]:
    df[col] = (
        df[col]
        .str.extract(r"\(([\d\.]+)%\)", expand=False)
        .astype(float)
    )

# =============================================================================
# CLEAN LABELS
# =============================================================================

label_map = {
    "UTR3": "3' UTR",
    "UTR5": "5' UTR",
    "downstream": "Downstream",
    "exonic": "Exonic",
    "intergenic": "Intergenic",
    "intronic": "Intronic",
    "ncRNA_exonic": "ncRNA Exonic",
    "ncRNA_intronic": "ncRNA Intronic",
    "upstream": "Upstream"
}

df["Label"] = (
    df["Annotation"]
    .map(label_map)
    .fillna(df["Annotation"])
)
# =============================================================================
# ORDER SEGMENTS
# =============================================================================

df["MeanPct"] = (
    df["Concordant"] +
    df["Discordant"]
) / 2

df = df.sort_values(
    "MeanPct",
    ascending=False
)

# =============================================================================
# GREYSCALE PALETTE
# =============================================================================

colors = [
    "#2F2F2F",
    "#4A4A4A",
    "#666666",
    "#808080",
    "#999999",
    "#B3B3B3",
    "#CCCCCC",
    "#DDDDDD",
    "#EEEEEE"
]

color_map = dict(
    zip(df["Label"], colors[:len(df)])
)

# =============================================================================
# FIGURE
# =============================================================================

fig, ax = plt.subplots(
    figsize=(10, 3)
)

# =============================================================================
# STACKED BARS
# =============================================================================

groups = ["Concordant", "Discordant"]
y_positions = [1, 0]

for group, y in zip(groups, y_positions):

    left = 0

    for _, row in df.iterrows():

        value = row[group]

        ax.barh(
            y=y,
            width=value,
            left=left,
            height=0.55,
            color=color_map[row["Label"]],
            edgecolor="white",
            linewidth=1
        )

        # Label only major segments

        if value >= 8:

            text_color = (
                "white"
                if row["Label"] in [
                    "Intronic",
                    "Intergenic",
                    "ncRNA Intronic"
                ]
                else TEXT
            )

            short_label = row["Label"]

            if short_label == "ncRNA Intronic":
                short_label = "ncRNA"

            ax.text(
                left + value/2,
                y,
                short_label,
                ha="center",
                va="center",
                fontsize=16,
                color=text_color,
                fontweight="bold"
            )

        left += value

# =============================================================================
# AXES
# =============================================================================

ax.set_xlim(0, 100)

ax.set_xticks([0, 25, 50, 75, 100])

ax.set_xticklabels(
    ["0%", "25%", "50%", "75%", "100%"],
    fontsize=TICK_SIZE,
    color=TEXT
)

ax.set_yticks(y_positions)

ax.set_yticklabels(
    groups,
    fontsize=LABEL_SIZE,
    color=TEXT,
    fontweight="bold"
)

ax.set_xlabel(
    "Percentage of unique lead SNPs",
    fontsize=16,
    color=TEXT
)

ax.set_title(
    "A. Functional Annotation Profiles of Concordant and Discordant Lead Variants",
    fontsize=TITLE_SIZE,
    fontweight="bold",
    color=TEXT,
    pad=20
)

# =============================================================================
# CLEANUP
# =============================================================================

ax.tick_params(colors=TEXT)

sns.despine(
    left=True,
    bottom=False
)

plt.tight_layout()

# =============================================================================
# SAVE
# =============================================================================

plt.savefig(
    FIGURES_DIR / "Figure4A_FunctionalAnnotationComposition.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    FIGURES_DIR / "Figure4A_FunctionalAnnotationComposition.svg",
    bbox_inches="tight"
)

# =============================================================================
# PANEL B
# Biological Theme Heatmap
# =============================================================================

summary = pd.read_csv(
    THEME_TABLE
)

# =============================================================================
# RETAIN ALL 15 PREDEFINED THEMES
# =============================================================================

plot_df = summary.copy()

plot_df = (
    plot_df[
        [
            "Theme",
            "Positive_Percent",
            "Negative_Percent",
            "Delta_Percent"
        ]
    ]
    .drop_duplicates()
)

plot_df["Theme"] = (
    plot_df["Theme"]
    .str.replace(r"\\n", "\n", regex=True)
)

if len(plot_df) != 15:

    raise RuntimeError(
        "QC FAIL: expected 15 GWAS Catalog themes, "
        f"found {len(plot_df)}"
    )

# =============================================================================
# ORDER BY TOTAL IMPORTANCE
# =============================================================================

plot_df["Total"] = (
    plot_df["Positive_Percent"] +
    plot_df["Negative_Percent"]
)

plot_df = (
    plot_df
    .sort_values(
        "Total",
        ascending=False
    )
)

# =============================================================================
# MATRIX
# =============================================================================

heatmap_df = (
    plot_df[
        [
            "Theme",
            "Positive_Percent",
            "Negative_Percent",
            "Delta_Percent"
        ]
    ]
    .set_index("Theme")
)

# =============================================================================
# DESCRIPTIVE DIFFERENCE
# =============================================================================

max_abs_delta = (
    heatmap_df["Delta_Percent"]
    .abs()
    .max()
)

# =============================================================================
# CUSTOM MATRIX
# =============================================================================

from matplotlib.patches import Rectangle
from matplotlib import cm

fig, ax = plt.subplots(
    figsize=(10, 10)
)

plt.subplots_adjust(left=0.45, right=0.98, bottom=0.12)

n_rows = len(heatmap_df)

red_cmap = cm.Reds
blue_cmap = cm.Blues

max_pct = max(
    heatmap_df["Positive_Percent"].max(),
    heatmap_df["Negative_Percent"].max(),
)

for i, (theme, row) in enumerate(heatmap_df.iterrows()):

    y = n_rows - i - 1

    pos = row["Positive_Percent"]
    neg = row["Negative_Percent"]
    delta = row["Delta_Percent"]

    # Positive
    ax.add_patch(
        Rectangle(
            (0, y),
            1,
            1,
            facecolor=red_cmap(
                0.15 + 0.85 * pos / max_pct
            ),
            edgecolor="white",
            linewidth=2
        )
    )

    # Negative
    ax.add_patch(
        Rectangle(
            (1, y),
            1,
            1,
            facecolor=blue_cmap(
                0.15 + 0.85 * neg/max_pct
            ),
            edgecolor="white",
            linewidth=2
        )
    )

    delta_intensity = (
        abs(delta) / max_abs_delta
        if max_abs_delta > 0
        else 0
    )

    grey_value = (
        0.95
        - 0.75 * delta_intensity
    )

    ax.add_patch(
        Rectangle(
            (2, y),
            1,
            1,
            facecolor=(
                grey_value,
                grey_value,
                grey_value
            ),
            edgecolor="lightgrey",
            linewidth=1.5
        )
    )


    # Numbers

    pos_color = (
        "white"
        if pos > 0.65 * max_pct
        else TEXT
    )

    ax.text(
        0.5,
        y + 0.5,
        f"{pos:.1f}",
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold",
        color=pos_color
    )

    neg_color = (
        "white"
        if neg > 0.65 * max_pct
        else TEXT
    )

    ax.text(
        1.5,
        y + 0.5,
        f"{neg:.1f}",
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold",
        color=neg_color
    )

    delta_label = f"{delta:+.1f}"

    delta_text_color = (
        "white"
        if grey_value < 0.5
        else TEXT
    )

    ax.text(
        2.5,
        y + 0.5,
        delta_label,
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold",
        color=delta_text_color
    )

ax.set_xlim(0, 3.05)
ax.set_ylim(0, n_rows)

ax.set_xticks([0.5, 1.5, 2.5])

ax.set_xticklabels(
    [
        "Concordant\n(%)",
        "Discordant\n(%)",
        "Δ percentage\npoints"
    ],
    fontsize=16,
    fontweight="bold"
)

ax.tick_params(axis="x", length=0, pad=8)

ax.set_yticks(
    np.arange(n_rows) + 0.5
)

display_labels = []

for theme in heatmap_df.index[::-1]:
    display_labels.append(theme)

ax.set_yticklabels(
    display_labels,
    fontsize=20
)

ax.tick_params(
    length=0
)

for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()

ax.set_title(
    "B. Recurrent GWAS Catalog Trait-Enrichment Profiles",
    fontsize=22,
    fontweight="bold",
    color=TEXT,
    pad=20
)

# =============================================================================
# SAVE
# =============================================================================

plt.savefig(
    FIGURES_DIR / "Figure4B_BiologicalThemesHeatmap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    FIGURES_DIR / "Figure4B_BiologicalThemesHeatmap.svg",
    bbox_inches="tight"
)

plt.show()

# =============================================================================
# PANEL C
# Disease-specific association strength by pleiotropic direction
# =============================================================================

TEXT = "#333333"

# =============================================================================
# LOAD DATA
# =============================================================================

raw = pd.read_csv(
    ASSOCIATION_STRENGTH_DIR
    / "AssociationStrength_LDPruned_Primary_SNPs.csv.gz",
    sep="\t",
    compression="gzip"
)

stats = pd.read_csv(
    FINAL_TABLES_DIR
    / "Table6_AssociationStrengthDirectionality_noChol.csv"
)

# Figure 4C shows the primary LD-pruned analysis.
stats = stats[
    stats["Analysis"] == "Primary: LD-pruned"
].copy()

# =============================================================================
# HARMONIZE DISEASE NAMES
# =============================================================================

raw["Disease"] = (
    raw["Disease"]
    .str.replace("_d", "", regex=False)
)

stats["Disease"] = stats["Disease"].replace({
    "Coronary Artery Disease": "CAD",
    "Hypertension": "HT",
    "Stroke": "STR",
    "Type 2 Diabetes": "T2D"
})

# =============================================================================
# ORDER
# =============================================================================

disease_order = [
    "CAD",
    "HT",
    "STR",
    "T2D"
]

disease_labels = {
    "CAD": "Coronary\nArtery Disease",
    "HT": "Hypertension",
    "STR": "Stroke",
    "T2D": "Type 2\nDiabetes"
}

# =============================================================================
# FIGURE
# =============================================================================

fig, ax = plt.subplots(
    figsize=(10, 8)
)

# =============================================================================
# VIOLINS
# =============================================================================

for i, disease in enumerate(disease_order):

    y = len(disease_order) - i - 1

    negative = raw[
        (raw["Disease"] == disease)
        &
        (raw["classification"] == "Negative")
    ]["abs_z"]

    positive = raw[
        (raw["Disease"] == disease)
        &
        (raw["classification"] == "Positive")
    ]["abs_z"]

    # --------------------------------------------------
    # Negative (left)
    # --------------------------------------------------

    if len(negative) > 1:

        vp = ax.violinplot(
            -negative,
            positions=[y],
            vert=False,
            widths=1,
            showmeans=False,
            showmedians=False,
            showextrema=False
        )

        for body in vp["bodies"]:

            body.set_facecolor("#377EB8")
            body.set_edgecolor("#377EB8")
            body.set_alpha(0.50)

    # --------------------------------------------------
    # Positive (right)
    # --------------------------------------------------

    if len(positive) > 1:

        vp = ax.violinplot(
            positive,
            positions=[y],
            vert=False,
            widths=1,
            showmeans=False,
            showmedians=False,
            showextrema=False
        )

        for body in vp["bodies"]:

            body.set_facecolor("#E41A1C")
            body.set_edgecolor("#E41A1C")
            body.set_alpha(0.60)

# =============================================================================
# DUMBBELLS
# =============================================================================

for i, disease in enumerate(disease_order):

    y = len(disease_order) - i - 1

    row = stats[
        stats["Disease"] == disease
    ].iloc[0]

    n_neg = int(row["N discordant signals"])
    n_pos = int(row["N concordant signals"])

    neg_med = -row["Discordant median |Z|"]
    pos_med = row["Concordant median |Z|"]

    delta = abs(neg_med) - pos_med

    ax.plot(
        [neg_med, pos_med],
        [y, y],
        color="#333333",
        linewidth=3,
        zorder=5
    )

    ax.scatter(
        neg_med,
        y,
        s=160,
        color=TEXT,
        edgecolor="white",
        linewidth=1.5,
        zorder=6
    )

    ax.scatter(
        pos_med,
        y,
        s=160,
        color=TEXT,
        edgecolor="white",
        linewidth=1.5,
        zorder=6
    )

    # --------------------------------------------------
    # MEDIAN LABELS
    # --------------------------------------------------

    ax.text(
        neg_med - 0.35,
        y,
        f"{abs(neg_med):.2f}",
        ha="right",
        va="center",
        fontsize=20,
        color=TEXT,
        fontweight="bold"
    )

    ax.text(
        pos_med + 0.35,
        y,
        f"{pos_med:.2f}",
        ha="left",
        va="center",
        fontsize=20,
        color=TEXT,
        fontweight="bold"
    )

# --------------------------------------------------
# SIGNIFICANCE
# --------------------------------------------------

    # --------------------------------------------------
    # SIGNIFICANCE
    # --------------------------------------------------

    p = row["FDR-adjusted P-value"]

    if p < 0.001:
        stars = "***"
    elif p < 0.01:
        stars = "**"
    elif p < 0.05:
        stars = "*"
    else:
        stars = ""

    if stars:

        ax.text(
            0,
            y + 0.2,
            stars,
            ha="center",
            va="bottom",
            fontsize=30,
            fontweight="bold",
            color=TEXT
        )

        ax.text(
            0,
            y + 0.1,
            f"Δ={delta:.2f}",
            ha="center",
            va="bottom",
            fontsize=20,
            fontweight="bold",
            color=TEXT
        )

# =============================================================================
# AXES
# =============================================================================

ax.axvline(
    0,
    color="#BDBDBD",
    linestyle="--",
    linewidth=1
)

max_z = np.ceil(
    np.percentile(
        raw["abs_z"],
        99
    )
)

ax.set_xlim(
    -max_z,
    max_z
)

ticks = np.arange(
    -max_z,
    max_z + 0.1,
    2
)

ax.set_xticks(ticks)

ax.set_xticklabels(
    [str(int(abs(x))) for x in ticks],
    fontsize=18
)

ax.set_yticks(
    np.arange(len(disease_order))
)

ax.set_yticklabels(
    [disease_labels[d] for d in disease_order[::-1]],
    fontsize=24,
    fontweight="bold"
)

ax.set_xlabel(
    "Disease-specific association strength (|Z|)",
    fontsize=20,
    fontweight="bold"
)

ax.set_ylabel("")

ax.set_title(
    "C. Disease-specific Association Strength of Pleiotropic Signals",
    fontsize=28,
    fontweight="bold",
    pad=5
)

# =============================================================================
# DIRECTION LABELS
# =============================================================================

ax.set_ylim(-0.7, len(disease_order)-0.2 + 0.5)

# =============================================================================
# CLEANUP
# =============================================================================

sns.despine(
    left=True,
    bottom=False
)

ax.tick_params(
    labelsize=18
)

plt.tight_layout()

# =============================================================================
# SAVE
# =============================================================================

plt.savefig(
    FIGURES_DIR / "Figure4C_AssociationStrength.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    FIGURES_DIR / "Figure4C_AssociationStrength.svg",
    bbox_inches="tight"
)

plt.show()
