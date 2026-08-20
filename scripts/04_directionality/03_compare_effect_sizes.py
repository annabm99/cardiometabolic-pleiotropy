import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# =========================================================
# CONFIG
# =========================================================

TRANSLATION = {
    "CAD": "Coronary Artery Disease",
    "HT": "Hypertension",
    "STR": "Stroke",
    "T2D": "Type 2 Diabetes"
}

# =========================================================
# FUNCTIONS
# =========================================================

def load_dataset(path):

    return pd.read_csv(
        path,
        sep='\t',
        compression='gzip'
    )


def compute_stats(pos_df, neg_df):

    pos = pos_df["abs_z"].dropna()
    neg = neg_df["abs_z"].dropna()

    if len(pos) == 0 or len(neg) == 0:
        return None

    # -----------------------------------------------------
    # Mann-Whitney U test
    # -----------------------------------------------------

    stat, pval = mannwhitneyu(
        pos,
        neg,
        alternative='two-sided'
    )

    # -----------------------------------------------------
    # Results
    # -----------------------------------------------------

    results = {

        # Counts
        "N_positive_pleiotropic_SNPs":
            len(pos),

        "N_negative_pleiotropic_SNPs":
            len(neg),

        # Medians
        "Positive_median_abs_zscore":
            round(pos.median(), 3),

        "Negative_median_abs_zscore":
            round(neg.median(), 3),

        # Means
        "Positive_mean_abs_zscore":
            round(pos.mean(), 3),

        "Negative_mean_abs_zscore":
            round(neg.mean(), 3),

        # Statistics
        "Mann_Whitney_U":
            round(stat, 3),

        "P_value":
            pval
    }

    return results


# =========================================================
# MAIN
# =========================================================

input_dir = sys.argv[1].rstrip("/")
out_dir = sys.argv[2].rstrip("/")

os.makedirs(out_dir, exist_ok=True)

results = []

# =========================================================
# LOOP DISEASES
# =========================================================

for short in ["CAD", "HT", "STR", "T2D"]:

    print(f"\nProcessing: {short}")

    # -----------------------------------------------------
    # Input paths
    # -----------------------------------------------------

    pos_path = os.path.join(
        input_dir,
        f"{short}_positive_snps_noChol.csv.gz"
    )

    neg_path = os.path.join(
        input_dir,
        f"{short}_negative_snps_noChol.csv.gz"
    )

    # -----------------------------------------------------
    # Check existence
    # -----------------------------------------------------

    if not os.path.exists(pos_path):

        print(f"Missing: {pos_path}")
        continue

    if not os.path.exists(neg_path):

        print(f"Missing: {neg_path}")
        continue

    # -----------------------------------------------------
    # Load datasets
    # -----------------------------------------------------

    pos_df = load_dataset(pos_path)
    neg_df = load_dataset(neg_path)

    # -----------------------------------------------------
    # Compute statistics
    # -----------------------------------------------------

    stats = compute_stats(
        pos_df,
        neg_df
    )

    if stats is None:
        continue

    stats["Disease"] = TRANSLATION[short]

    results.append(stats)

# =========================================================
# GLOBAL ANALYSIS
# =========================================================

print("\nProcessing: GLOBAL")

global_pos_path = os.path.join(
    input_dir,
    "All_positive_snps_noChol.csv.gz"
)

global_neg_path = os.path.join(
    input_dir,
    "All_negative_snps_noChol.csv.gz"
)

if os.path.exists(global_pos_path) and os.path.exists(global_neg_path):

    global_pos_df = load_dataset(global_pos_path)
    global_neg_df = load_dataset(global_neg_path)

    global_stats = compute_stats(
        global_pos_df,
        global_neg_df
    )

    global_stats["Disease"] = "All diseases"

    results.append(global_stats)

# =========================================================
# FINAL DATAFRAME
# =========================================================

stats_df = pd.DataFrame(results)

# ---------------------------------------------------------
# Multiple testing correction
# ---------------------------------------------------------

if len(stats_df) > 0:

    # FDR correction
    stats_df["P_FDR"] = multipletests(
        stats_df["P_value"],
        method='fdr_bh'
    )[1]

    # Bonferroni correction
    stats_df["P_Bonferroni"] = multipletests(
        stats_df["P_value"],
        method='bonferroni'
    )[1]

# ---------------------------------------------------------
# Column order
# ---------------------------------------------------------

stats_df = stats_df[[
    "Disease",

    "N_positive_pleiotropic_SNPs",
    "N_negative_pleiotropic_SNPs",

    "Positive_median_abs_zscore",
    "Negative_median_abs_zscore",

    "Positive_mean_abs_zscore",
    "Negative_mean_abs_zscore",

    "Mann_Whitney_U",

    "P_value",
    "P_FDR",
    "P_Bonferroni"
]]

# =========================================================
# EXPORT
# =========================================================

csv_path = os.path.join(
    out_dir,
    "Directionality_StatisticalComparison_noChol.csv"
)

stats_df.to_csv(
    csv_path,
    index=False
)

print(f"\nSaved: {csv_path}")

print("\nDONE")