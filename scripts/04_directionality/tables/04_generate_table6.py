import os
import sys
import pandas as pd

# =========================================================
# MAIN
# =========================================================

PROJECT_DIR = os.environ.get(
    "CVP_PROJECT_DIR",
    "/path/to/cardiovascular_pleiotropies"
)

input_file = os.environ.get(
    "CVP_DIRECTIONALITY_STATS",
    f"{PROJECT_DIR}/3-Directionality/3-StatisticalComparison/Directionality_StatisticalComparison_noChol.csv"
)
out_dir = os.environ.get(
    "CVP_FINAL_TABLES_DIR",
    f"{PROJECT_DIR}/FinalTables"
)

os.makedirs(out_dir, exist_ok=True)

# =========================================================
# LOAD
# =========================================================

df = pd.read_csv(input_file)

# =========================================================
# SIGNIFICANCE COLUMN
# =========================================================

df["Significant"] = df["P_FDR"] < 0.05

df["Significant"] = df["Significant"].map({
    True: "Yes",
    False: "No"
})

# =========================================================
# SELECT TABLE 6 COLUMNS
# =========================================================

table6 = df[[

    "Disease",

    "N_positive_pleiotropic_SNPs",
    "N_negative_pleiotropic_SNPs",

    "Positive_median_abs_zscore",
    "Negative_median_abs_zscore",

    "Mann_Whitney_U",

    "P_value",
    "P_FDR",
    "P_Bonferroni",

    "Significant"

]].copy()

# =========================================================
# OPTIONAL RENAMING
# =========================================================

table6.columns = [

    "Disease",

    "N positive pleiotropic SNPs",
    "N negative pleiotropic SNPs",

    "Positive median |z|",
    "Negative median |z|",

    "Mann-Whitney U",

    "P-value",
    "FDR-adjusted P-value",
    "Bonferroni-adjusted P-value",

    "Significant"
]

# =========================================================
# EXPORT
# =========================================================

out_prefix = os.path.join(
    out_dir,
    "Table6_EffectSizeDirectionality_noChol"
)

# ---------------------------------------------------------
# CSV
# ---------------------------------------------------------

csv_path = f"{out_prefix}.csv"

table6.to_csv(
    csv_path,
    index=False
)

# ---------------------------------------------------------
# Excel
# ---------------------------------------------------------

xlsx_path = f"{out_prefix}.xlsx"

with pd.ExcelWriter(
    xlsx_path,
    engine='openpyxl'
) as writer:

    table6.to_excel(
        writer,
        index=False,
        sheet_name='Table6'
    )

# ---------------------------------------------------------
# LaTeX
# ---------------------------------------------------------

tex_path = f"{out_prefix}.tex"

latex = table6.to_latex(
    index=False,
    escape=False,
    float_format="%.3g"
)

with open(tex_path, 'w') as f:
    f.write(latex)

# =========================================================
# DONE
# =========================================================

print(f"Saved: {csv_path}")
print(f"Saved: {xlsx_path}")
print(f"Saved: {tex_path}")

print("\nDONE")
