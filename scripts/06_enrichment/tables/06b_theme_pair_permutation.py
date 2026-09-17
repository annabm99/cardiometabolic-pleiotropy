#!/usr/bin/env python3

"""
Pair-level permutation test for biological-theme differences between
concordant and discordant pleiotropy.

Rationale
---------
The original Figure 4B analysis used summed phenotype-pair recurrence
counts in Fisher exact tests. These counts are not independent because
the same phenotype pair can contribute to multiple correlated GWAS
Catalog pathways.

This script uses the phenotype pair as the permutation unit. For each
permutation, the complete positive/negative enrichment profiles of a
phenotype pair are swapped together, preserving within-pair pathway
correlation.

The analysis:
    - starts from significant nonredundant FUMA enrichment results;
    - excludes HDL/LDL pairs;
    - retains GWAS Catalog enrichments;
    - recomputes pathway recurrence separately by direction;
    - retains pathways recurrent across >=3 diseases;
    - summarizes theme composition;
    - tests observed positive-vs-negative percentage differences using
      pair-level permutation;
    - applies Benjamini-Hochberg FDR correction.

Output
------
Table10_ThemePermutation.csv
"""

import ast
import os
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


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
        "CVP_NONREDUNDANT_PATHWAYS",
        str(
            PROJECT_DIR
            / "5-EnrichmentAnalysis/3-Collapse/"
            "all_fuma_pathways_nonredundant.tsv"
        )
    )
)

FINAL_TABLES_DIR = Path(
    os.environ.get(
        "CVP_FINAL_TABLES_DIR",
        str(PROJECT_DIR / "FinalTables")
    )
)

THEME_SCRIPT = Path(__file__).with_name(
    "06_generate_table10_themes.py"
)

OUTPUT_FILE = (
    FINAL_TABLES_DIR
    / "Table10_ThemePermutation.csv"
)

N_PERMUTATIONS = 10000
RANDOM_SEED = 20260917

EXCLUDED_TRAITS = {
    "HDL_t",
    "LDL_t",
}


# =============================================================================
# LOAD CURRENT THEME ASSIGNMENT FUNCTION
# =============================================================================

def load_assign_theme():

    source = THEME_SCRIPT.read_text()

    tree = ast.parse(source)

    func_node = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "assign_theme"
    )

    module = ast.Module(
        body=[func_node],
        type_ignores=[]
    )

    namespace = {}

    exec(
        compile(
            ast.fix_missing_locations(module),
            filename=str(THEME_SCRIPT),
            mode="exec"
        ),
        namespace
    )

    return namespace["assign_theme"]


# =============================================================================
# THEME SUMMARY
# =============================================================================

def calculate_theme_summary(data):

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

    FINAL_TABLES_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    assign_theme = load_assign_theme()

    df = pd.read_csv(
        INPUT_FILE,
        sep="\t"
    )

    df = df[
        ~df["phenotype_pair"].apply(
            lambda x: any(
                trait in x
                for trait in EXCLUDED_TRAITS
            )
        )
    ].copy()

    df = df[
        df["category"] == "GWAScatalog"
    ].copy()

    df["disease"] = (
        df["phenotype_pair"]
        .str.split("_d-")
        .str[0]
    )

    df["Theme"] = (
        df["pathway_name"]
        .apply(assign_theme)
    )

    pairs = sorted(
        df["phenotype_pair"].unique()
    )

    print(
        f"Phenotype pairs: {len(pairs)}"
    )

    print(
        f"GWAS Catalog rows: {len(df):,}"
    )

    print(
        "Unique pathways: "
        f"{df['pathway_name'].nunique():,}"
    )

    # =========================================================================
    # OBSERVED
    # =========================================================================

    observed = calculate_theme_summary(
        df
    )

    themes = sorted(
        observed.keys()
    )

    # =========================================================================
    # PERMUTATIONS
    # =========================================================================

    rng = np.random.default_rng(
        RANDOM_SEED
    )

    null_deltas = {
        theme: []
        for theme in themes
    }

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
            .isin(swap_pairs)
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
            ].append(delta)

        if (
            (iteration + 1) % 1000
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
                    np.median(null),

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
        "\n"
        "================================="
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
                f"{x:.4g}",
        )
    )

    print(
        f"\nSaved:\n{OUTPUT_FILE}"
    )


if __name__ == "__main__":
    main()