#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgba
import os
from pathlib import Path
import numpy as np

##############################################################################
# INPUT / OUTPUT
##############################################################################

PROJECT_DIR = Path(
    os.environ.get("CVP_PROJECT_DIR", "/path/to/cardiovascular_pleiotropies")
)

TABLE4 = os.environ.get(
    "CVP_TABLE4_NO_HDL_LDL",
    str(PROJECT_DIR / "FinalTables/Table4_PleioLoci_noHDL_LDL.csv")
)

HIGHLIGHT_GENES = os.environ.get(
    "CVP_FIGURE2_HIGHLIGHT_GENES",
    str(PROJECT_DIR / "FinalTables/Figure2_HighlightGenes_noHDL_LDL.csv")
)

OUTDIR = Path(
    os.environ.get("CVP_FIGURES_DIR", str(PROJECT_DIR / "Figures"))
)

OUTDIR.mkdir(parents=True, exist_ok=True)

##############################################################################
# LOAD TABLE
##############################################################################

df = pd.read_csv(TABLE4)
highlight = pd.read_csv(HIGHLIGHT_GENES)

##############################################################################
# SORT LOCI
##############################################################################

def chr_sort(x):
    try:
        return int(x)
    except:
        return 999

df["chr_order"] = (
    df["Chromosome"]
    .astype(str)
    .apply(chr_sort)
)

df = (
    df.sort_values(
        ["chr_order", "Merged_Start"]
    )
    .reset_index(drop=True)
)

##############################################################################
# DISEASES
##############################################################################

DISEASES = [
    "CAD",
    "HT",
    "STR",
    "T2D"
]

DISEASE_LABELS = {
    "CAD": "Coronary\nArtery Disease",
    "HT": "Hypertension",
    "STR": "Stroke",
    "T2D": "Type 2\nDiabetes"
}

##############################################################################
# VARIABLE LOCUS WIDTHS
##############################################################################

WIDTH_MAP = {
    1:0.1,
    2: 1,
    3: 3,
    4: 10
}

CHR_GAP = 5

df["PlotWidth"] = (
    df["Number_of_Diseases"]
    .map(WIDTH_MAP)
    .fillna(1.0)
)

highlight_loci = set()

for _, row in highlight.iterrows():

    for ml in str(
        row["MergedLocusIDs"]
    ).split("|"):

        highlight_loci.add(
            ml.strip()
        )

df.loc[
    df["MergedLocusID"].isin(
        highlight_loci
    ),
    "PlotWidth"
] = 10

x_starts = []
chrom_boundaries = []
chrom_centers = []

current_x = 0

for chrom, sub in df.groupby(
    "Chromosome",
    sort=False
):

    chrom_boundaries.append(current_x)

    chr_start = current_x

    for _, row in sub.iterrows():

        x_starts.append(current_x)

        current_x += row["PlotWidth"]

    chr_end = current_x

    chrom_centers.append(
        (
            chrom,
            (chr_start + chr_end) / 2
        )
    )

    current_x += CHR_GAP

chrom_boundaries.append(current_x)

df["x_start"] = x_starts

df["x_end"] = (
    df["x_start"]
    + df["PlotWidth"]
)

locus_lookup = {}

for _, row in df.iterrows():

    locus_lookup[
        row["MergedLocusID"]
    ] = {

        "x_start":
            row["x_start"],

        "x_end":
            row["x_end"]
    }

total_width = current_x

##############################################################################
# ENLARGE LOCI USED FOR GENE ANNOTATIONS
##############################################################################

print(
    df[
        df["MergedLocusID"].isin(
            ["ML_504","ML_505","ML_506"]
        )
    ][[
        "MergedLocusID",
        "PlotWidth",
        "Chromosome"
    ]]
)

##############################################################################
# COLORS
##############################################################################

RED_BASE = np.array(
    to_rgba("#E41A1C")
)

BLUE_BASE = np.array(
    to_rgba("#0072B2")
)

WHITE = np.array(
    to_rgba("white")
)

def blend(base_color, n_diseases):

    alpha_map = {
        1:0.65,
        2: 0.75,
        3: 0.85,
        4: 0.95
    }

    alpha = alpha_map.get(
        n_diseases,
        0.35
    )

    return (
        WHITE * (1 - alpha)
        + base_color * alpha
    )

##############################################################################
# GENE ANNOTATION SPANS
##############################################################################

gene_spans = []

for _, row in highlight.iterrows():

    gene = row["Label"]

    loci = [
        x.strip()
        for x in str(
            row["MergedLocusIDs"]
        ).split("|")
    ]

    xs = []
    xe = []

    for ml in loci:

        if ml not in locus_lookup:
            continue

        xs.append(
            locus_lookup[ml]["x_start"]
        )

        xe.append(
            locus_lookup[ml]["x_end"]
        )

    if len(xs) == 0:
        continue

    gene_spans.append({

        "gene":
            gene,

        "xmin":
            min(xs),

        "xmax":
            max(xe)
    })

    if gene == "AS3MT/CNNM2/NT5C2/WBP1L/MARCKSL1P1":

        gene_spans[-1]["gene"] = (
            "AS3MT/CNNM2/NT5C2/\n"
            "WBP1L/MARCKSL1P1"
        )
        gene_spans[-1]["label_shift"] = -15

    elif gene == "OVOL1/PCNXL3":

        gene_spans[-1]["label_shift"] = 13

    else:

        gene_spans[-1]["label_shift"] = 0

##############################################################################
# PACK INTO TRACKS
##############################################################################

gene_spans = sorted(
    gene_spans,
    key=lambda x: x["xmin"]
)

tracks = []

for span in gene_spans:

    placed = False

    for track in tracks:

        last = track[-1]

        if span["xmin"] > (
            last["xmax"] + 2
        ):

            track.append(span)
            placed = True
            break

    if not placed:
        tracks.append([span])

n_annotation_tracks = len(tracks)

##############################################################################
# FIGURE SIZE
##############################################################################

fig_width = 18

fig_height = 2.5

fig, ax = plt.subplots(
    figsize=(fig_width, fig_height)
)

##############################################################################
# DRAW LOCI
##############################################################################

for _, row in df.iterrows():

    x0 = row["x_start"]

    width = row["PlotWidth"]

    diseases_present = set(
        str(
            row["Diseases_Involved"]
        ).split("|")
    )

    n_diseases = row["Number_of_Diseases"]

    if row["Pleiotropy_Type"] == "Concordant":

        color = blend(
            RED_BASE,
            n_diseases
        )

    else:

        color = blend(
            BLUE_BASE,
            n_diseases
        )

    for y, disease in enumerate(DISEASES):

        if disease in diseases_present:
            facecolor = color
        else:
            facecolor = "white"

        rect = Rectangle(
            (
                x0,
                len(DISEASES) - 1 - y
            ),
            width,
            1,
            facecolor=facecolor,
            edgecolor="none"
        )

        ax.add_patch(rect)
    
##############################################################################
# GENE ANNOTATIONS
##############################################################################

base_y = -0.15

for track_idx, track in enumerate(tracks):

    y = (
        base_y
        - track_idx * 0.45
    )

    for span in track:

        xmin = span["xmin"]
        xmax = span["xmax"]

        ax.plot(
            [xmin, xmax],
            [y, y],
            color="black",
            linewidth=5,
            solid_capstyle="butt"
        )

        xmid = (
            (xmin + xmax) / 2
            + span.get("label_shift", 0)
        )

        ax.text(
            xmid,
            y - 0.15,
            span["gene"],
            ha="center",
            va="top",
            fontsize=12
        )

##############################################################################
# CHROMOSOME SEPARATORS
##############################################################################

# for boundary in chrom_boundaries:

#     ax.axvline(
#         boundary,
#         color="black",
#         linewidth=0.6
#     )

##############################################################################
# AXES
##############################################################################

ax.set_xlim(
    0,
    total_width
)

bottom_margin = (
    -0.7
    - n_annotation_tracks * 0.45
)

ax.set_ylim(
    bottom_margin,
    len(DISEASES)
)

ax.set_xticks([])

ax.set_yticks(
    np.arange(len(DISEASES))
    + 0.5
)

ax.set_yticklabels(
    [
        DISEASE_LABELS[x]
        for x in DISEASES[::-1]
    ],
    fontsize=12,
    fontweight="bold"
)

##############################################################################
# CHROMOSOME LABELS
##############################################################################

problem_chr = {
    "13","14","18", "19",
    "20","21","22"
}

for chrom, center in chrom_centers:

    rot = 45 if str(chrom) in problem_chr else 0

    ax.text(
        center,
        len(DISEASES)+0.12,
        f"chr{chrom}",
        ha="center",
        va="bottom",
        fontsize=10,
        rotation=rot
    )

##############################################################################
# TITLE
##############################################################################

ax.set_title(
    "Pleiotropic Locus Sharing Across Cardiovascular Diseases",
    fontsize=14,
    fontweight="bold",
    pad=70
)

##############################################################################
# CLEAN LOOK
##############################################################################

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)

##############################################################################
# SAVE
##############################################################################

pdf_file = (
    OUTDIR /
    "Figure2_LocusMatrix_WidthWeighted.pdf"
)

svg_file = (
    OUTDIR /
    "Figure2_LocusMatrix_WidthWeighted.svg"
)

png_file = (
    OUTDIR /
    "Figure2_LocusMatrix_WidthWeighted.png"
)

plt.savefig(
    pdf_file,
    bbox_inches="tight"
)

plt.savefig(
    svg_file,
    bbox_inches="tight"
)

plt.savefig(
    png_file,
    dpi=600,
    bbox_inches="tight"
)

plt.close()

##############################################################################
# FIGURE 2 LEGEND
##############################################################################

fig_leg, ax_leg = plt.subplots(
    figsize=(7, 2.5)
)

ax_leg.axis("off")

##############################################################################
# TITLE
##############################################################################

ax_leg.text(
    0,
    1.85,
    "Pleiotropy degree (number of diseases)",
    fontsize=18,
    fontweight="bold"
)

##############################################################################
# SAME WIDTH RATIOS AS MAIN FIGURE
##############################################################################

legend_widths = {
    1: 0.15,
    2: 0.50,
    3: 0.90,
    4: 1.30
}

##############################################################################
# CONCORDANT
##############################################################################

ax_leg.text(
    0,
    1.15,
    "Concordant pleiotropies",
    fontsize=13,
    va="center"
)

x = 4.0

for nd in [1,2,3,4]:

    w = legend_widths[nd]

    ax_leg.add_patch(
        Rectangle(
            (x, 1.00),
            w,
            0.28,
            facecolor=blend(
                RED_BASE,
                nd
            ),
            edgecolor="none"
        )
    )

    ax_leg.text(
        x + w/2,
        0.83,
        str(nd),
        ha="center",
        fontsize=11
    )

    x += w + 0.40

##############################################################################
# DISCORDANT
##############################################################################

ax_leg.text(
    0,
    0.35,
    "Discordant pleiotropies",
    fontsize=13,
    va="center"
)

x = 4.0

for nd in [1,2,3,4]:

    w = legend_widths[nd]

    ax_leg.add_patch(
        Rectangle(
            (x, 0.20),
            w,
            0.28,
            facecolor=blend(
                BLUE_BASE,
                nd
            ),
            edgecolor="none"
        )
    )

    ax_leg.text(
        x + w/2,
        0.03,
        str(nd),
        ha="center",
        fontsize=11
    )

    x += w + 0.40

##############################################################################
# LIMITS
##############################################################################

ax_leg.set_xlim(
    -0.1,
    9.0
)

ax_leg.set_ylim(
    -0.25,
    2.15
)

##############################################################################
# SAVE
##############################################################################

legend_pdf = (
    OUTDIR /
    "Figure2_Legend.pdf"
)

legend_svg = (
    OUTDIR /
    "Figure2_Legend.svg"
)

legend_png = (
    OUTDIR /
    "Figure2_Legend.png"
)

fig_leg.savefig(
    legend_pdf,
    bbox_inches="tight"
)

fig_leg.savefig(
    legend_svg,
    bbox_inches="tight"
)

fig_leg.savefig(
    legend_png,
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig_leg)

##############################################################################
# QC
##############################################################################

print("\n====================================")
print("FIGURE CREATED")
print("====================================\n")

print(f"Loci plotted: {len(df):,}")
print(f"Total plot width: {total_width:.1f}")

print("\nOutputs:")
print(pdf_file)
print(svg_file)
print(png_file)

print("\nDone.")
