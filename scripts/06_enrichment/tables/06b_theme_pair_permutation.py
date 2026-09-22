#!/usr/bin/env python3

"""
Pair-level permutation test for recurrent GWAS Catalog trait-enrichment
profiles in concordant versus discordant pleiotropy.

Scientific question
-------------------
Among GWAS Catalog trait annotations that recur across at least three
diseases, does thematic composition differ between concordant and
discordant pleiotropy beyond what is expected after phenotype-pair-level
direction-label swapping?

Design
------
- Start from Table 8, the validated primary 24-pair FUMA enrichment table.
- Retain GWAS Catalog enrichments.
- Exclude eight explicitly ambiguous/cross-domain source terms.
- Assign the remaining terms to the predefined trait themes.
- Remove terms that cannot be assigned to one of those themes BEFORE
  observed and permuted recurrence calculations.
- For observed data and every permutation:
    * calculate recurrence separately by direction;
    * retain terms enriched in >=3 diseases;
    * aggregate phenotype-pair contributions by theme.
- Permute by swapping the complete positive/negative enrichment profile
  of each phenotype pair, preserving within-pair correlation.
- Use two-sided empirical permutation P values and BH FDR correction.

Terminology
-----------
Positive = concordant pleiotropy.
Negative = discordant pleiotropy.

The themes represent GWAS Catalog phenotype/trait annotations, not
molecular pathways.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from theme_classification import (
    AMBIGUOUS_THEME_EXCLUSIONS,
    assign_theme,
)


# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_DIR = Path(
    os.environ.get(
        "CVP_PROJECT_DIR",
        "/path/to/cardiovascular_pleiotropies"
    )
)

INPUT_FILE = Path(
    os.environ.get(
        "CVP_TABLE8_PATHWAYS",
        str(
            PROJECT_DIR
            / "FinalTables"
            / "Table8_FUMA_Pathway_Enrichment.csv"
        )
    )
)

OUTPUT_FILE = Path(
    os.environ.get(
        "CVP_THEME_PERMUTATION_OUTPUT",
        str(
            PROJECT_DIR
            / "FinalTables"
            / "Table10_ThemePermutation.csv"
        )
    )
)

N_PERMUTATIONS = int(
    os.environ.get(
        "CVP_THEME_N_PERMUTATIONS",
        "10000"
    )
)

RANDOM_SEED = int(
    os.environ.get(
        "CVP_THEME_RANDOM_SEED",
        "20260917"
    )
)

EXCLUDED_TRAITS = {
    "HDL_t",
    "LDL_t",
}


# =============================================================================
# RECURRENCE AND THEME SUMMARY
# =============================================================================

def calculate_recurrence(data):

    x = (
        data[
            [
                "phenotype_pair",
                "disease",
                "pleiotropy_direction",
                "pathway_name",
                "Theme",
            ]
        ]
        .drop_duplicates()
    )

    recurrence = (
        x.groupby(
            [
                "pleiotropy_direction",
                "pathway_name",
                "Theme",
            ],
            as_index=False,
        )
        .agg(
            N_Phenotype_Pairs=(
                "phenotype_pair",
                "nunique"
            ),
            N_Diseases=(
                "disease",
                "nunique"
            ),
        )
    )

    recurrence = recurrence[
        recurrence["N_Diseases"] >= 3
    ].copy()

    return recurrence


def calculate_theme_summary(data):

    recurrence = calculate_recurrence(
        data
    )

    summary = (
        recurrence.groupby(
            [
                "Theme",
                "pleiotropy_direction",
            ],
            as_index=False,
        )
        .agg(
            Total_Pairs=(
                "N_Phenotype_Pairs",
                "sum"
            )
        )
    )

    if summary.empty:
        return {}

    totals = (
        summary.groupby(
            "pleiotropy_direction"
        )["Total_Pairs"]
        .sum()
        .to_dict()
    )

    output = {}

    for theme in summary["Theme"].unique():

        subset = summary[
            summary["Theme"] == theme
        ]

        values = dict(
            zip(
                subset["pleiotropy_direction"],
                subset["Total_Pairs"]
            )
        )

        positive = values.get(
            "positive",
            0
        )

        negative = values.get(
            "negative",
            0
        )

        positive_total = totals.get(
            "positive",
            0
        )

        negative_total = totals.get(
            "negative",
            0
        )

        positive_percent = (
            100 * positive / positive_total
            if positive_total > 0
            else 0
        )

        negative_percent = (
            100 * negative / negative_total
            if negative_total > 0
            else 0
        )

        fold_change = (
            positive_percent / negative_percent
            if negative_percent > 0
            else np.inf
        )

        output[theme] = {
            "Positive_Percent":
                positive_percent,
            "Negative_Percent":
                negative_percent,
            "Delta_Percent":
                positive_percent
                - negative_percent,
            "Fold_Change":
                fold_change,
        }

    return output


# =============================================================================
# MAIN
# =============================================================================

def main():

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    # -------------------------------------------------------------------------
    # LOAD VALIDATED TABLE 8
    # -------------------------------------------------------------------------

    df = pd.read_csv(
        INPUT_FILE
    )

    required_columns = {
        "Phenotype_Pair",
        "Pleiotropy_Direction",
        "Pathway",
        "Category",
    }

    missing_columns = (
        required_columns
        - set(df.columns)
    )

    if missing_columns:

        raise RuntimeError(
            "Missing required Table 8 columns: "
            + ", ".join(
                sorted(missing_columns)
            )
        )

    df = df.rename(
        columns={
            "Phenotype_Pair":
                "phenotype_pair",
            "Pleiotropy_Direction":
                "pleiotropy_direction",
            "Pathway":
                "pathway_name",
            "Category":
                "category",
        }
    )

    # -------------------------------------------------------------------------
    # PRIMARY 24-PAIR ANALYSIS
    # -------------------------------------------------------------------------

    df = df[
        ~df["phenotype_pair"].apply(
            lambda x: any(
                trait in str(x)
                for trait in EXCLUDED_TRAITS
            )
        )
    ].copy()

    df = df[
        df["category"] == "GWAScatalog"
    ].copy()

    df["disease"] = (
        df["phenotype_pair"]
        .str.split(
            "_d-",
            n=1
        )
        .str[0]
    )

    print(
        "=========================================="
    )
    print(
        "GWAS CATALOG THEME PERMUTATION INPUT"
    )
    print(
        "=========================================="
    )

    print(
        f"Phenotype pairs: "
        f"{df['phenotype_pair'].nunique()}"
    )

    print(
        f"GWAS Catalog rows before "
        f"thematic eligibility: {len(df):,}"
    )

    print(
        f"Unique GWAS Catalog terms: "
        f"{df['pathway_name'].nunique():,}"
    )

    # -------------------------------------------------------------------------
    # EXPLICIT AMBIGUOUS / CROSS-DOMAIN EXCLUSIONS
    # -------------------------------------------------------------------------

    ambiguous_mask = (
        df["pathway_name"]
        .isin(
            AMBIGUOUS_THEME_EXCLUSIONS
        )
    )

    print(
        f"Rows removed by 8 explicit "
        f"term exclusions: "
        f"{ambiguous_mask.sum():,}"
    )

    print(
        f"Unique explicit excluded terms present: "
        f"{df.loc[ambiguous_mask, 'pathway_name'].nunique()}"
    )

    df = df[
        ~ambiguous_mask
    ].copy()

    # -------------------------------------------------------------------------
    # THEME CLASSIFICATION
    # -------------------------------------------------------------------------

    df["Theme"] = (
        df["pathway_name"]
        .apply(
            assign_theme
        )
    )

    unclassified_mask = (
        df["Theme"]
        == "UNCLASSIFIED"
    )

    print(
        f"Unclassified rows removed: "
        f"{unclassified_mask.sum():,}"
    )

    print(
        f"Unique unclassified terms removed: "
        f"{df.loc[unclassified_mask, 'pathway_name'].nunique():,}"
    )

    df = df[
        ~unclassified_mask
    ].copy()

    # -------------------------------------------------------------------------
    # HARD QC FOR THE VALIDATED PRIMARY DATASET
    # -------------------------------------------------------------------------

    n_pairs = (
        df["phenotype_pair"]
        .nunique()
    )

    n_rows = len(df)

    n_terms = (
        df["pathway_name"]
        .nunique()
    )

    directions = set(
        df["pleiotropy_direction"]
        .unique()
    )

    if n_pairs != 24:

        raise RuntimeError(
            f"QC FAIL: expected 24 phenotype pairs, "
            f"found {n_pairs}"
        )

    if n_rows != 2999:

        raise RuntimeError(
            f"QC FAIL: expected 2,999 eligible rows, "
            f"found {n_rows:,}"
        )

    if n_terms != 419:

        raise RuntimeError(
            f"QC FAIL: expected 419 eligible unique terms, "
            f"found {n_terms}"
        )

    if directions != {
        "positive",
        "negative",
    }:

        raise RuntimeError(
            "QC FAIL: unexpected pleiotropy directions: "
            f"{sorted(directions)}"
        )

    print(
        "\nEligible thematic universe:"
    )

    print(
        f"  rows: {n_rows:,}"
    )

    print(
        f"  unique terms: {n_terms:,}"
    )

    print(
        f"  phenotype pairs: {n_pairs}"
    )

    print(
        f"  themes represented before recurrence: "
        f"{df['Theme'].nunique()}"
    )

    # =========================================================================
    # OBSERVED RECURRENCE
    # =========================================================================

    observed_recurrence = (
        calculate_recurrence(
            df
        )
    )

    observed_recurrent_rows = len(
        observed_recurrence
    )

    observed_recurrent_terms = (
        observed_recurrence[
            "pathway_name"
        ]
        .nunique()
    )

    contribution_totals = (
        observed_recurrence
        .groupby(
            "pleiotropy_direction"
        )["N_Phenotype_Pairs"]
        .sum()
        .to_dict()
    )

    if observed_recurrent_rows != 224:

        raise RuntimeError(
            "QC FAIL: expected 224 observed "
            "direction-specific recurrent records, "
            f"found {observed_recurrent_rows}"
        )

    if observed_recurrent_terms != 170:

        raise RuntimeError(
            "QC FAIL: expected 170 observed recurrent "
            f"unique terms, found {observed_recurrent_terms}"
        )

    if contribution_totals.get(
        "positive"
    ) != 1163:

        raise RuntimeError(
            "QC FAIL: expected 1,163 positive "
            "recurrent pair contributions, found "
            f"{contribution_totals.get('positive')}"
        )

    if contribution_totals.get(
        "negative"
    ) != 522:

        raise RuntimeError(
            "QC FAIL: expected 522 negative "
            "recurrent pair contributions, found "
            f"{contribution_totals.get('negative')}"
        )

    print(
        "\nObserved recurrent catalogue "
        "(N_Diseases >= 3):"
    )

    print(
        f"  direction-specific records: "
        f"{observed_recurrent_rows}"
    )

    print(
        f"  unique recurrent terms: "
        f"{observed_recurrent_terms}"
    )

    print(
        f"  positive/concordant pair contributions: "
        f"{contribution_totals['positive']}"
    )

    print(
        f"  negative/discordant pair contributions: "
        f"{contribution_totals['negative']}"
    )

    observed = calculate_theme_summary(
        df
    )

    themes = sorted(
        observed.keys()
    )

    if len(themes) != 15:

        raise RuntimeError(
            f"QC FAIL: expected 15 observed themes, "
            f"found {len(themes)}"
        )

    # =========================================================================
    # PERMUTATIONS
    # =========================================================================

    pairs = sorted(
        df["phenotype_pair"]
        .unique()
    )

    rng = np.random.default_rng(
        RANDOM_SEED
    )

    null_deltas = {
        theme: []
        for theme in themes
    }

    print(
        "\nRunning "
        f"{N_PERMUTATIONS:,} phenotype-pair "
        "label-swap permutations..."
    )

    for iteration in range(
        N_PERMUTATIONS
    ):

        perm = df.copy()

        swap_pairs = {
            pair
            for pair in pairs
            if rng.random() < 0.5
        }

        mask = (
            perm["phenotype_pair"]
            .isin(
                swap_pairs
            )
        )

        perm.loc[
            mask,
            "pleiotropy_direction"
        ] = (
            perm.loc[
                mask,
                "pleiotropy_direction"
            ]
            .map(
                {
                    "positive":
                        "negative",
                    "negative":
                        "positive",
                }
            )
        )

        perm_summary = (
            calculate_theme_summary(
                perm
            )
        )

        for theme in themes:

            delta = (
                perm_summary
                .get(
                    theme,
                    {
                        "Delta_Percent":
                            0
                    }
                )
                ["Delta_Percent"]
            )

            null_deltas[
                theme
            ].append(
                delta
            )

        if (
            (iteration + 1)
            % 1000
            == 0
        ):

            print(
                "Completed "
                f"{iteration + 1:,} / "
                f"{N_PERMUTATIONS:,} "
                "permutations"
            )

    # =========================================================================
    # RESULTS
    # =========================================================================

    rows = []

    for theme in themes:

        obs = observed[
            theme
        ]

        observed_delta = (
            obs[
                "Delta_Percent"
            ]
        )

        null = np.asarray(
            null_deltas[
                theme
            ]
        )

        permutation_p = (
            1
            + np.sum(
                np.abs(null)
                >= abs(
                    observed_delta
                )
            )
        ) / (
            N_PERMUTATIONS
            + 1
        )

        rows.append(
            {
                "Theme":
                    theme,

                "Positive_Percent":
                    obs[
                        "Positive_Percent"
                    ],

                "Negative_Percent":
                    obs[
                        "Negative_Percent"
                    ],

                "Delta_Percent":
                    observed_delta,

                "Fold_Change":
                    obs[
                        "Fold_Change"
                    ],

                "Permutation_P":
                    permutation_p,

                "Null_Delta_Median":
                    np.median(
                        null
                    ),

                "Null_Delta_2.5pct":
                    np.quantile(
                        null,
                        0.025
                    ),

                "Null_Delta_97.5pct":
                    np.quantile(
                        null,
                        0.975
                    ),
            }
        )

    results = pd.DataFrame(
        rows
    )

    results[
        "Permutation_FDR"
    ] = multipletests(
        results[
            "Permutation_P"
        ],
        method="fdr_bh",
    )[1]

    results = (
        results.sort_values(
            [
                "Permutation_FDR",
                "Permutation_P",
            ]
        )
        .reset_index(
            drop=True
        )
    )

    results.to_csv(
        OUTPUT_FILE,
        index=False
    )

    print(
        "\n================================="
    )

    print(
        "PAIR-LEVEL PERMUTATION RESULTS"
    )

    print(
        "================================="
    )

    print(
        results.to_string(
            index=False,
            float_format=lambda x:
                f"{x:.6g}",
        )
    )

    print(
        f"\nSaved:\n{OUTPUT_FILE}"
    )


if __name__ == "__main__":

    main()
