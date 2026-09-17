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
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


TRANSLATION = {
    "CAD": "Coronary Artery Disease",
    "HT": "Hypertension",
    "STR": "Stroke",
    "T2D": "Type 2 Diabetes",
}

N_DOWNSAMPLE_ITERATIONS = 5000
RANDOM_SEED = 20260917


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

    # Positive values mean that discordant observations
    # tend to have larger |Z| than concordant observations.
    rank_biserial = (
        1
        - (2 * stat) / (len(pos) * len(neg))
    )

    return {
        "N_concordant": len(pos),
        "N_discordant": len(neg),

        "Concordant_median_abs_z":
            pos.median(),

        "Discordant_median_abs_z":
            neg.median(),

        "Concordant_mean_abs_z":
            pos.mean(),

        "Discordant_mean_abs_z":
            neg.mean(),

        "Mann_Whitney_U":
            stat,

        "P_value":
            pval,

        "Rank_biserial_discordant_gt_concordant":
            rank_biserial,
    }

def downsample_power_diagnostic(
    pos_df,
    neg_df,
    target_n,
    observed_pure_p,
    observed_pure_effect,
    seed,
):

    """
    Downsample the original discordant group to the size of the
    pure-discordant group.

    This estimates how much weakening of the Mann-Whitney P-value
    would be expected from sample-size reduction alone, while
    retaining the original discordant-group composition.
    """

    pos = (
        pos_df["abs_z"]
        .dropna()
        .to_numpy()
    )

    neg = (
        neg_df["abs_z"]
        .dropna()
        .to_numpy()
    )

    if target_n > len(neg):
        raise ValueError(
            "Target downsample size exceeds discordant sample size"
        )

    rng = np.random.default_rng(seed)

    rows = []

    for i in range(N_DOWNSAMPLE_ITERATIONS):

        sampled_neg = rng.choice(
            neg,
            size=target_n,
            replace=False
        )

        stat, pval = mannwhitneyu(
            pos,
            sampled_neg,
            alternative="two-sided"
        )

        rank_biserial = (
            1
            - (2 * stat) /
            (len(pos) * len(sampled_neg))
        )

        rows.append({
            "Iteration": i + 1,
            "P_value": pval,
            "Rank_biserial":
                rank_biserial,
        })

    reps = pd.DataFrame(rows)

    summary = {

        "Downsample_iterations":
            N_DOWNSAMPLE_ITERATIONS,

        "Downsample_target_N":
            target_n,

        "Downsample_median_P":
            reps["P_value"].median(),

        "Downsample_P_2.5pct":
            reps["P_value"].quantile(0.025),

        "Downsample_P_97.5pct":
            reps["P_value"].quantile(0.975),

        "Downsample_percent_P_lt_0.05":
            100 * (
                reps["P_value"] < 0.05
            ).mean(),

        "Downsample_percent_P_ge_observed_pure":
            100 * (
                reps["P_value"] >= observed_pure_p
            ).mean(),

        "Downsample_median_rank_biserial":
            reps["Rank_biserial"].median(),

        "Downsample_rank_biserial_2.5pct":
            reps["Rank_biserial"].quantile(0.025),

        "Downsample_rank_biserial_97.5pct":
            reps["Rank_biserial"].quantile(0.975),

        "Observed_pure_P":
            observed_pure_p,

        "Observed_pure_rank_biserial":
            observed_pure_effect,
    }

    return summary, reps


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

    power_rows = []
    power_replicates = []

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

        # -------------------------------------------------
        # POWER / SAMPLE-SIZE DIAGNOSTIC
        # -------------------------------------------------

        power_summary, power_reps = (
            downsample_power_diagnostic(
                pos_df=pos,
                neg_df=neg,
                target_n=len(pure_neg),
                observed_pure_p=pure["P_value"],
                observed_pure_effect=(
                    pure[
                        "Rank_biserial_discordant_gt_concordant"
                    ]
                ),
                seed=RANDOM_SEED + len(power_rows),
            )
        )

        power_summary["Disease"] = disease_label

        power_summary[
            "Original_rank_biserial"
        ] = original[
            "Rank_biserial_discordant_gt_concordant"
        ]

        power_rows.append(
            power_summary
        )

        power_reps["Disease"] = disease_label

        power_replicates.append(
            power_reps
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

            "Rank_biserial_all_discordant":
                original[
                    "Rank_biserial_discordant_gt_concordant"
                ],

            "Rank_biserial_pure_discordant":
                pure[
                    "Rank_biserial_discordant_gt_concordant"
                ],
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

    # =====================================================
    # POWER / SAMPLE-SIZE DIAGNOSTIC
    # =====================================================

    power = pd.DataFrame(
        power_rows
    )

    power.to_csv(
        os.path.join(
            out_dir,
            "MixedDirectionality_PowerDiagnostic.csv"
        ),
        index=False
    )

    power_reps = pd.concat(
        power_replicates,
        ignore_index=True
    )

    power_reps.to_csv(
        os.path.join(
            out_dir,
            "MixedDirectionality_PowerDiagnostic_Replicates.csv"
        ),
        index=False
    )

    print("\nPower / sample-size diagnostic:")
    print(
        power.round(4).to_string(
            index=False
        )
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