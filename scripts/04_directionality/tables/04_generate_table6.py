#!/usr/bin/env python3

import os
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(
    os.environ.get(
        "CVP_PROJECT_DIR",
        "/path/to/cardiovascular_pleiotropies",
    )
)

INPUT_FILE = Path(
    os.environ.get(
        "CVP_ASSOCIATION_STRENGTH_STATS",
        str(
            PROJECT_DIR
            / "3-Directionality"
            / "3-AssociationStrength"
            / "AssociationStrength_LDPruned_Summary.csv"
        ),
    )
)

OUT_DIR = Path(
    os.environ.get(
        "CVP_FINAL_TABLES_DIR",
        str(PROJECT_DIR / "FinalTables"),
    )
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# =========================================================
# LOAD
# =========================================================

df = pd.read_csv(
    INPUT_FILE
)

expected_analyses = {
    "Primary_LD_pruned",
    "Sensitivity_LD_pruned_no_mixed",
}

observed_analyses = set(
    df["Analysis"]
)

missing = (
    expected_analyses
    - observed_analyses
)

if missing:
    raise ValueError(
        f"Missing expected analyses: {sorted(missing)}"
    )


# =========================================================
# LABEL ANALYSES
# =========================================================

analysis_labels = {
    "Primary_LD_pruned":
        "Primary: LD-pruned",
    "Sensitivity_LD_pruned_no_mixed":
        "Sensitivity: LD-pruned, mixed-context excluded",
}

df["Analysis"] = (
    df["Analysis"]
    .map(analysis_labels)
)

df["Significant_after_FDR"] = (
    df["P_FDR"] < 0.05
).map({
    True: "Yes",
    False: "No",
})


# =========================================================
# SELECT + RENAME
# =========================================================

table6 = df[[
    "Analysis",
    "Disease",
    "N_concordant",
    "N_discordant",
    "Concordant_median_abs_z",
    "Discordant_median_abs_z",
    "Delta_median_abs_z_discordant_minus_concordant",
    "Mann_Whitney_U",
    "P_value",
    "P_FDR",
    "P_Bonferroni",
    "Significant_after_FDR",
]].copy()

table6.columns = [
    "Analysis",
    "Disease",
    "N concordant signals",
    "N discordant signals",
    "Concordant median |Z|",
    "Discordant median |Z|",
    "Difference in median |Z| (discordant - concordant)",
    "Mann-Whitney U",
    "P-value",
    "FDR-adjusted P-value",
    "Bonferroni-adjusted P-value",
    "Significant after FDR",
]


# =========================================================
# ORDER
# =========================================================

analysis_order = [
    "Primary: LD-pruned",
    "Sensitivity: LD-pruned, mixed-context excluded",
]

disease_order = [
    "Coronary Artery Disease",
    "Hypertension",
    "Stroke",
    "Type 2 Diabetes",
]

table6["Analysis"] = pd.Categorical(
    table6["Analysis"],
    categories=analysis_order,
    ordered=True,
)

table6["Disease"] = pd.Categorical(
    table6["Disease"],
    categories=disease_order,
    ordered=True,
)

table6 = (
    table6
    .sort_values(
        ["Analysis", "Disease"]
    )
    .reset_index(drop=True)
)


# =========================================================
# EXPORT
# =========================================================

out_prefix = (
    OUT_DIR
    / "Table6_AssociationStrengthDirectionality_noChol"
)

csv_path = Path(
    f"{out_prefix}.csv"
)

tex_path = Path(
    f"{out_prefix}.tex"
)

table6.to_csv(
    csv_path,
    index=False,
)

latex = table6.to_latex(
    index=False,
    escape=False,
    float_format="%.3g",
)

tex_path.write_text(
    latex
)

print(f"Saved: {csv_path}")
print(f"Saved: {tex_path}")

print("\n=== TABLE 6 ===")
print(
    table6.to_string(
        index=False
    )
)

print("\nDONE")
