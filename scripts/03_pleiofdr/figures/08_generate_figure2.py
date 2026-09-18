#!/usr/bin/env python3

"""
Figure 2 – Pleiotropic Architecture Matrix

Inputs
------
Fig2_DiseaseBurden.csv
Fig2_MatrixCounts.csv

Outputs
-------
Figure2_PleioMatrix.svg
Figure2_PleioMatrix.png

Figure2_legend.svg
Figure2_legend.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

# ============================================================
# PATHS
# ============================================================

PROJECT_DIR = os.environ.get(
    "CVP_PROJECT_DIR",
    "/path/to/cardiovascular_pleiotropies"
)

INPUT_DIR = os.environ.get(
    "CVP_FINAL_TABLES_DIR",
    f"{PROJECT_DIR}/FinalTables"
)

OUTPUT_DIR = os.environ.get(
    "CVP_FIGURES_DIR",
    f"{PROJECT_DIR}/Figures"
)

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)

BURDEN_PATH = os.path.join(
    INPUT_DIR,
    "Fig2_DiseaseBurden.csv"
)

MATRIX_PATH = os.path.join(
    INPUT_DIR,
    "Fig2_MatrixCounts.csv"
)

# =============================================================================
# CONSTANTS
# =============================================================================

RED = "#E41A1C"
BLUE = "#0072B2"

DISEASES = ["CAD", "HT", "STR", "T2D"]

DISEASE_LABELS = {
    "CAD": "Coronary\nArtery Disease",
    "HT": "Hypertension",
    "STR": "Stroke",
    "T2D": "Type 2 Diabetes"
}

TRAIT_ORDER = [
    "BMI",
    "WC",
    "SBP",
    "DBP",
    "FG",
    "TGL",
    "HDL",
    "LDL"
]

TRAIT_LABELS = {
    "BMI": "Body Mass\nIndex",
    "WC": "Waist\nCircumference",
    "SBP": "Systolic Blood\nPressure",
    "DBP": "Diastolic Blood\nPressure",
    "FG": "Fasting\nGlucose",
    "TGL": "Triglycerides",
    "HDL": "HDL\nCholesterol",
    "LDL": "LDL\nCholesterol"
}

# =============================================================================
# AXIS TRANSFORMATIONS
# =============================================================================

def transform_cad(x):
    return x

def transform_str(x):
    return x

HT_GAP = 12

def transform_ht(x):
    x = np.asarray(x)

    return np.where(
        x <= 50,
        x,
        50 + HT_GAP + (x - 150)
    )

T2D_GAP1 = 12
T2D_GAP2 = 12

def transform_t2d(x):

    x = np.asarray(x)

    result = np.zeros_like(x, dtype=float)

    # Segment 1: 0–125
    mask1 = x <= 125
    result[mask1] = x[mask1]

    # Segment 2: 300–325
    mask2 = (x > 125) & (x <= 325)
    result[mask2] = (
        125
        + T2D_GAP1
        + (x[mask2] - 300)
    )

    # Segment 3: 375–390
    mask3 = x > 325
    result[mask3] = (
        125
        + T2D_GAP1
        + 25
        + T2D_GAP2
        + (x[mask3] - 375)
    )

    return result

def add_axis_break(ax, x_pos, dx=2.0):

    trans = ax.get_xaxis_transform()

    ax.add_patch(
        plt.Rectangle(
            (x_pos+3.2 - 3*dx, -0.01),
            3*dx,
            0.02,
            transform=trans,
            facecolor="white",
            edgecolor="none",
            clip_on=False,
            zorder=10
        )
    )

    ax.plot(
        [x_pos - dx, x_pos],
        [-0.02, 0.02],
        transform=trans,
        color="0.5",
        lw=1.2,
        clip_on=False,
        zorder=11
    )

    ax.plot(
        [x_pos + 0.4*dx, x_pos + 1.4*dx],
        [-0.02, 0.02],
        transform=trans,
        color="0.5",
        lw=1.2,
        clip_on=False,
        zorder=11
    )

TRANSFORM = {
    "CAD": transform_cad,
    "HT": transform_ht,
    "STR": transform_str,
    "T2D": transform_t2d
}

RAW_MAX = {
    "CAD": 59+5,
    "HT": 168+5,
    "STR": 42+5,
    "T2D": 390+5
}

DISPLAY_MAX = {
    d: float(TRANSFORM[d](RAW_MAX[d]))
    for d in DISEASES
}

#f proportional panel widths
WIDTH_RATIOS = [DISPLAY_MAX[d] for d in DISEASES]

TICKS = {

    "CAD": [0, 25, 50, 75],

    "STR": [0, 25, 50],

    "HT": [0, 25, 50, 150, 175],

    "T2D": [
    0, 25, 50, 75, 100, 125,
    300, 325,
    375, 400
]
}

# =============================================================================
# LOAD DATA
# =============================================================================

burden = pd.read_csv(BURDEN_PATH)
matrix = pd.read_csv(MATRIX_PATH)

# =============================================================================
# FIGURE
# =============================================================================

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 11
})

fig = plt.figure(figsize=(14, 8))

gs = GridSpec(
    nrows=2,
    ncols=4,
    figure=fig,
    height_ratios=[1.1, 5],
    width_ratios=WIDTH_RATIOS,
    hspace=0.15,
    wspace=0.08
)

# =============================================================================
# BOTTOM ROW : MATRIX
# =============================================================================

bottom_axes = []

for i, disease in enumerate(DISEASES):

    if i == 0:
        ax = fig.add_subplot(gs[1, i])
    else:
        ax = fig.add_subplot(gs[1, i], sharey=bottom_axes[0])

    bottom_axes.append(ax)

    dsub = matrix.loc[matrix["Disease"] == disease].copy()
    dsub["Trait"] = pd.Categorical(
        dsub["Trait"],
        categories=TRAIT_ORDER,
        ordered=True
    )
    dsub = dsub.sort_values("Trait")

    trans = TRANSFORM[disease]

    y = np.array([
        0,
        1,
        3,   # gap after WC
        4,
        6,   # gap after DBP
        7,
        8,
        9
    ])
    concordant = trans(dsub["Concordant"].values)
    discordant = trans(dsub["Discordant"].values)

    ax.barh(
        y - 0.18,
        concordant,
        height=0.32,
        color=RED,
        zorder=3,
        alpha=0.9
    )

    ax.barh(
        y + 0.18,
        discordant,
        height=0.32,
        color=BLUE,
        zorder=3,
        alpha=0.9
    )

    ax.set_xlim(0, DISPLAY_MAX[disease])

    if disease == "HT":
        add_axis_break(ax, 56)

    if disease == "T2D":
        add_axis_break(ax, 131)
        add_axis_break(ax, 168)

    tick_positions = trans(np.array(TICKS[disease]))
    tick_labels = [str(x) for x in TICKS[disease]]

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)

    ax.grid(
        axis="x",
        color="0.85",
        linewidth=0.8,
        zorder=0
    )

    if i == 0:

        ax.set_yticks(y)
        ax.set_yticklabels(
            [TRAIT_LABELS[t] for t in TRAIT_ORDER]
        )

    else:

        ax.tick_params(
            axis="y",
            left=False,
            labelleft=False
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# reverse trait order so BMI at top
bottom_axes[0].invert_yaxis()

ax0 = bottom_axes[0]

line_x = -0.76
label_x = -0.81

# ==========================
# Anthropometric
# ==========================

ax0.plot(
    [line_x, line_x],
    [0.97, 0.77],
    transform=ax0.transAxes,
    color="black",
    lw=1.5,
    clip_on=False,
    alpha=0.4
)


ax0.text(
    label_x,
    0.845,
    "ANTHROPO-\nMETRIC",
    transform=ax0.transAxes,
    ha="right",
    va="center",
    fontsize=8,
    color="black",
    fontweight="semibold",
    alpha=0.4
)

# ==========================
# Blood pressure
# ==========================

ax0.plot(
    [line_x, line_x],
    [0.69, 0.49],
    transform=ax0.transAxes,
    color="black",
    lw=1.5,
    clip_on=False,
    alpha=0.4
)

ax0.text(
    label_x,
    0.59,
    "BLOOD\nPRESSURE",
    transform=ax0.transAxes,
    ha="right",
    va="center",
    fontsize=8,
    color="black",
    fontweight="semibold",
    alpha=0.4
)

# ==========================
# Metabolic
# ==========================

ax0.plot(
    [line_x, line_x],
    [0.41, 0.02],
    transform=ax0.transAxes,
    color="black",
    lw=1.5,
    clip_on=False,
    alpha=0.4
)

ax0.text(
    label_x,
    0.215,
    "METABOLIC",
    transform=ax0.transAxes,
    ha="right",
    va="center",
    fontsize=8,
    color="black",
    fontweight="semibold",
    alpha=0.4
)

# =============================================================================
# LAYOUT
# =============================================================================

fig.subplots_adjust(
    left=0.23,
    right=0.98,
    top=0.74,
    bottom=0.10
)

# =============================================================================
# TOP BURDEN PANELS (FIXED WIDTH)
# =============================================================================

matrix_top = max(ax.get_position().y1 for ax in bottom_axes)

BURDEN_Y = matrix_top + 0.055
BURDEN_W = 0.12
BURDEN_H = 0.14

for disease, matrix_ax in zip(DISEASES, bottom_axes):

    row = burden.loc[burden["Disease"] == disease].iloc[0]

    concordant = row["Concordant"]
    discordant = row["Discordant"]

    bbox = matrix_ax.get_position()

    center_x = (bbox.x0 + bbox.x1) / 2

    burden_ax = fig.add_axes([
        center_x - BURDEN_W / 2,
        BURDEN_Y,
        BURDEN_W,
        BURDEN_H
    ])

    MAX_BURDEN = max(
        burden["Concordant"].max(),
        burden["Discordant"].max()
    )

    burden_ax.set_ylim(0, MAX_BURDEN * 1.08)

    bars = burden_ax.bar(
        [0.40, 0.80],
        [concordant, discordant],
        color=[RED, BLUE],
        width=0.35,
        alpha=0.9
    )

    burden_ax.set_xlim(0.10, 1.10)

    for bar, color in zip(bars, [RED, BLUE]):

        h = bar.get_height()

        burden_ax.text(
            bar.get_x() + bar.get_width()/2,
            h + MAX_BURDEN * 0.015,
            f"{int(h)}",
            ha="center",
            va="bottom",
            fontsize=10,
            color=color,
            fontweight="bold"
        )

    burden_ax.set_xticks([])
    burden_ax.set_yticks([])

    for spine in burden_ax.spines.values():
        spine.set_visible(False)

# =============================================================================
# DISEASE LABELS BETWEEN ROWS
# =============================================================================

LABEL_Y = matrix_top + 0.03

for disease, matrix_ax in zip(DISEASES, bottom_axes):

    bbox = matrix_ax.get_position()

    center_x = (bbox.x0 + bbox.x1) / 2

    fig.text(
        center_x,
        LABEL_Y,
        DISEASE_LABELS[disease],
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold"
    )

# =============================================================================
# EXPORT
# =============================================================================

svg_file = os.path.join(
    OUTPUT_DIR,
    "Figure2_PleioMatrix.svg"
)

png_file = os.path.join(
    OUTPUT_DIR,
    "Figure2_PleioMatrix.png"
)

fig.savefig(
    svg_file,
    bbox_inches="tight"
)

fig.savefig(
    png_file,
    dpi=600,
    bbox_inches="tight"
)

plt.close()

print("Saved:")
print(svg_file)
print(png_file)
