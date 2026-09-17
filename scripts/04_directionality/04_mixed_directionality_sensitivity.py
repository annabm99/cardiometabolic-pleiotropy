#!/usr/bin/env python3

"""
Sensitivity analysis for mixed-direction disease-SNP observations.

Primary disease-level classification:
    Concordant:
        no discordant phenotype-pair associations

    Discordant:
        at least one discordant phenotype-pair association

Sensitivity subclassification:
    Pure discordant:
        >=1 discordant pair and 0 concordant pairs

    Mixed discordant:
        >=1 discordant pair and >=1 concordant pair

This script:
    1. quantifies mixed discordant observations;
    2. exports all mixed disease-SNP observations;
    3. compares disease |Z| between concordant observations and:
       a) all discordant observations;
       b) pure-discordant observations only.
"""

import os
import sys

import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


TRANSLATION = {
    "CAD": "Coronary Artery Disease",
    "HT": "Hypertension",
    "STR": "Stroke",
    "T2D": "Type 2 Diabetes",
}


def load_dataset(path):
    return pd.read_csv(
        path,
        sep="\t",
        compression="gzip"
    )


def classify_discordant_subtypes(df):

    pure = df[
        (df["n_negative_pairs"] > 0) &
        (df["n_positive_pairs"] == 0)
    ].copy()

    mixed = df[
        (df["n_negative_pairs"] > 0) &
        (df["n_positive_pairs"] > 0)
    ].copy()

    return pure, mixed


def compare_abs_z(pos_df, neg_df):

    pos = pos_df["abs_z"].dropna()
    neg = neg_df["abs_z"].dropna()

    stat, pval = mannwhitneyu(
        pos,
        neg,
        alternative="two-sided"
    )

    return {
        "N_concordant": len(pos),
        "N_discordant": len(neg),
        "Concordant_median_abs_z": pos.median(),
        "Discordant_median_abs_z": neg.median(),
        "Concordant_mean_abs_z": pos.mean(),
        "Discordant_mean_abs_z": neg.mean(),
        "Mann_Whitney_U": stat,
        "P_value": pval,
    }


def main():

    input_dir = sys.argv[1].rstrip("/")
    out_dir = sys.argv[2].rstrip("/")

    os.makedirs(out_dir, exist_ok=True)

    datasets = [
        ("CAD", TRANSLATION["CAD"]),
        ("HT", TRANSLATION["HT"]),
        ("STR", TRANSLATION["STR"]),
        ("T2D", TRANSLATION["T2D"]),
        ("All", "All diseases"),
    ]

    summary_rows = []
    sensitivity_rows = []
    mixed_cases = []

    for prefix, disease_label in datasets:

        pos = load_dataset(
            os.path.join(
                input_dir,
                f"{prefix}_positive_snps_noChol.csv.gz"
            )
        )

        neg = load_dataset(
            os.path.join(
                input_dir,
                f"{prefix}_negative_snps_noChol.csv.gz"
            )
        )

        pure_neg, mixed = classify_discordant_subtypes(neg)

        # -------------------------------------------------
        # Mixed-direction summary
        # -------------------------------------------------

        summary_rows.append({
            "Disease": disease_label,
            "N_concordant": len(pos),
            "N_discordant_all": len(neg),
            "N_discordant_pure": len(pure_neg),
            "N_mixed_discordant": len(mixed),
            "Mixed_percent_of_discordant":
                100 * len(mixed) / len(neg)
                if len(neg) > 0
                else 0,
        })

        # Avoid duplicating the same observations from
        # disease-specific files and the global file.
        if prefix != "All":

            tmp = mixed.copy()

            tmp["Disease_Label"] = disease_label

            mixed_cases.append(tmp)

        # -------------------------------------------------
        # Original discordant comparison
        # -------------------------------------------------

        original = compare_abs_z(
            pos,
            neg
        )

        # -------------------------------------------------
        # Pure-discordant sensitivity
        # -------------------------------------------------

        pure = compare_abs_z(
            pos,
            pure_neg
        )

        sensitivity_rows.append({
            "Disease": disease_label,

            "N_concordant":
                original["N_concordant"],

            "N_discordant_all":
                original["N_discordant"],

            "N_discordant_pure":
                pure["N_discordant"],

            "N_mixed_removed":
                len(mixed),

            "Concordant_median_abs_z":
                original["Concordant_median_abs_z"],

            "Discordant_all_median_abs_z":
                original["Discordant_median_abs_z"],

            "Discordant_pure_median_abs_z":
                pure["Discordant_median_abs_z"],

            "Concordant_mean_abs_z":
                original["Concordant_mean_abs_z"],

            "Discordant_all_mean_abs_z":
                original["Discordant_mean_abs_z"],

            "Discordant_pure_mean_abs_z":
                pure["Discordant_mean_abs_z"],

            "P_all_discordant":
                original["P_value"],

            "P_pure_discordant":
                pure["P_value"],
        })

    # =====================================================
    # EXPORT SUMMARY
    # =====================================================

    summary = pd.DataFrame(
        summary_rows
    )

    summary.to_csv(
        os.path.join(
            out_dir,
            "MixedDirectionality_Summary.csv"
        ),
        index=False
    )

    # =====================================================
    # EXPORT CASES
    # =====================================================

    cases = pd.concat(
        mixed_cases,
        ignore_index=True
    )

    cases.to_csv(
        os.path.join(
            out_dir,
            "MixedDirectionality_Cases.csv"
        ),
        index=False
    )

    # =====================================================
    # SENSITIVITY STATISTICS
    # =====================================================

    sensitivity = pd.DataFrame(
        sensitivity_rows
    )

    sensitivity["P_all_FDR"] = multipletests(
        sensitivity["P_all_discordant"],
        method="fdr_bh"
    )[1]

    sensitivity["P_pure_FDR"] = multipletests(
        sensitivity["P_pure_discordant"],
        method="fdr_bh"
    )[1]

    sensitivity.to_csv(
        os.path.join(
            out_dir,
            "MixedDirectionality_AbsZ_Sensitivity.csv"
        ),
        index=False
    )

    print("\nMixed-directionality summary:")
    print(
        summary.round(3).to_string(
            index=False
        )
    )

    print("\n|Z| sensitivity:")
    print(
        sensitivity.round(4).to_string(
            index=False
        )
    )

    print("\nDONE")


if __name__ == "__main__":
    main()