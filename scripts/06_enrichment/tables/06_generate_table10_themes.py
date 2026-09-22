#!/usr/bin/env python3

"""
Generate the recurrent GWAS Catalog trait-enrichment catalogue and the
observed thematic summary used for Figure 4B.

Scientific question
-------------------
Which GWAS Catalog trait annotations recur across at least three of the
four diseases, separately for concordant and discordant pleiotropy, and
how are those recurrent annotations distributed across broad phenotype
themes?

Important
---------
These themes summarize GWAS Catalog phenotype/trait annotations. They
must not be interpreted as molecular pathway categories.

This script is descriptive only. It does NOT perform inferential tests.
Pair-level permutation inference is performed separately by:

    06b_theme_pair_permutation.py

Input
-----
Validated primary Table 8 FUMA enrichment results.

Outputs
-------
1. Recurrent GWAS Catalog term catalogue:
   - includes all terms recurrent across >=3 diseases;
   - retains ambiguous terms transparently;
   - records whether each term contributes to Figure 4B.

2. Observed 15-theme summary:
   - recurrent term counts;
   - phenotype-pair contribution counts;
   - percentages;
   - concordant-minus-discordant difference;
   - fold change.

Positive = concordant pleiotropy.
Negative = discordant pleiotropy.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

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

OUTPUT_CATALOGUE = Path(
    os.environ.get(
        "CVP_THEME_CATALOGUE_OUTPUT",
        str(
            PROJECT_DIR
            / "FinalTables"
            / "Table10_Recurrent_GWASCatalog_Terms.csv"
        )
    )
)

OUTPUT_SUMMARY = Path(
    os.environ.get(
        "CVP_THEME_SUMMARY_OUTPUT",
        str(
            PROJECT_DIR
            / "FinalTables"
            / "Table10_ThemeSummary_Observed.csv"
        )
    )
)

EXCLUDED_TRAITS = {
    "HDL_t",
    "LDL_t",
}

DIRECTION_LABELS = {
    "positive": "Concordant",
    "negative": "Discordant",
}


# =============================================================================
# HELPERS
# =============================================================================

def join_sorted(values):

    return "; ".join(
        sorted(
            {
                str(x)
                for x in values
                if pd.notna(x)
            }
        )
    )


def exclusion_reason(pathway, theme):

    if pathway in AMBIGUOUS_THEME_EXCLUSIONS:
        return AMBIGUOUS_THEME_EXCLUSIONS[pathway]

    if theme == "UNCLASSIFIED":
        return (
            "No clear assignment to one of the "
            "predefined GWAS Catalog trait themes."
        )

    return ""


# =============================================================================
# MAIN
# =============================================================================

def main():

    OUTPUT_CATALOGUE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    OUTPUT_SUMMARY.parent.mkdir(
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
        "N_Overlap",
        "P_Value",
        "FDR",
    }

    missing = (
        required_columns
        - set(df.columns)
    )

    if missing:

        raise RuntimeError(
            "Missing required Table 8 columns: "
            + ", ".join(
                sorted(missing)
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
            "N_Overlap":
                "n_overlap",
            "P_Value":
                "p_value",
            "FDR":
                "fdr",
        }
    )

    # -------------------------------------------------------------------------
    # PRIMARY 24-PAIR GWAS CATALOG DATA
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
        "========================================="
    )
    print(
        "RECURRENT GWAS CATALOG TRAIT ENRICHMENT"
    )
    print(
        "========================================="
    )

    print(
        f"GWAS Catalog rows: {len(df):,}"
    )

    print(
        f"Unique terms: "
        f"{df['pathway_name'].nunique():,}"
    )

    print(
        f"Phenotype pairs: "
        f"{df['phenotype_pair'].nunique()}"
    )

    # -------------------------------------------------------------------------
    # RECURRENCE ACROSS DISEASES
    # -------------------------------------------------------------------------

    recurrence = (
        df.groupby(
            [
                "pleiotropy_direction",
                "pathway_name",
            ],
            as_index=False,
        )
        .agg(
            N_Diseases=(
                "disease",
                "nunique"
            ),
            Diseases_Involved=(
                "disease",
                join_sorted
            ),
            N_Phenotype_Pairs=(
                "phenotype_pair",
                "nunique"
            ),
            Mean_N_Overlap=(
                "n_overlap",
                "mean"
            ),
            Best_P_Value=(
                "p_value",
                "min"
            ),
            Best_FDR=(
                "fdr",
                "min"
            ),
        )
    )

    recurrence = recurrence[
        recurrence["N_Diseases"] >= 3
    ].copy()

    # -------------------------------------------------------------------------
    # THEME ASSIGNMENT
    # -------------------------------------------------------------------------

    recurrence[
        "Assigned_Theme"
    ] = (
        recurrence[
            "pathway_name"
        ]
        .apply(
            assign_theme
        )
    )

    recurrence[
        "Explicit_Ambiguous_Exclusion"
    ] = (
        recurrence[
            "pathway_name"
        ]
        .isin(
            AMBIGUOUS_THEME_EXCLUSIONS
        )
    )

    recurrence[
        "Included_in_Figure4B"
    ] = np.where(
        (
            ~recurrence[
                "Explicit_Ambiguous_Exclusion"
            ]
        )
        & (
            recurrence[
                "Assigned_Theme"
            ]
            != "UNCLASSIFIED"
        ),
        "Yes",
        "No",
    )

    recurrence[
        "Exclusion_Reason"
    ] = recurrence.apply(
        lambda row:
            exclusion_reason(
                row["pathway_name"],
                row["Assigned_Theme"],
            ),
        axis=1,
    )

    # For the transparent supplement, ambiguous terms are not
    # presented as if they had a valid unique biological theme.
    recurrence[
        "Biological_Theme"
    ] = np.where(
        recurrence[
            "Included_in_Figure4B"
        ]
        == "Yes",
        recurrence[
            "Assigned_Theme"
        ],
        "Not assigned",
    )

    recurrence[
        "Pleiotropy_Direction_Label"
    ] = (
        recurrence[
            "pleiotropy_direction"
        ]
        .map(
            DIRECTION_LABELS
        )
    )

    # -------------------------------------------------------------------------
    # QC: RECURRENT CATALOGUE
    # -------------------------------------------------------------------------

    n_records = len(
        recurrence
    )

    n_unique_terms = (
        recurrence[
            "pathway_name"
        ]
        .nunique()
    )

    included = recurrence[
        recurrence[
            "Included_in_Figure4B"
        ]
        == "Yes"
    ].copy()

    excluded = recurrence[
        recurrence[
            "Included_in_Figure4B"
        ]
        == "No"
    ].copy()

    if n_records != 234:

        raise RuntimeError(
            "QC FAIL: expected 234 recurrent "
            "direction-specific records, "
            f"found {n_records}"
        )

    if n_unique_terms != 178:

        raise RuntimeError(
            "QC FAIL: expected 178 unique "
            "recurrent GWAS Catalog terms, "
            f"found {n_unique_terms}"
        )

    if len(included) != 224:

        raise RuntimeError(
            "QC FAIL: expected 224 recurrent "
            "records included in Figure 4B, "
            f"found {len(included)}"
        )

    if (
        included[
            "pathway_name"
        ]
        .nunique()
        != 170
    ):

        raise RuntimeError(
            "QC FAIL: expected 170 unique "
            "included recurrent terms"
        )

    if len(excluded) != 10:

        raise RuntimeError(
            "QC FAIL: expected 10 excluded "
            "direction-specific records, "
            f"found {len(excluded)}"
        )

    if (
        excluded[
            "pathway_name"
        ]
        .nunique()
        != 8
    ):

        raise RuntimeError(
            "QC FAIL: expected 8 unique "
            "excluded recurrent terms"
        )

    if (
        (
            excluded[
                "Assigned_Theme"
            ]
            == "UNCLASSIFIED"
        ).sum()
        != 0
    ):

        raise RuntimeError(
            "QC FAIL: observed recurrent catalogue "
            "contains non-ambiguous UNCLASSIFIED terms"
        )

    # -------------------------------------------------------------------------
    # SAVE TRANSPARENT RECURRENT-TERM CATALOGUE
    # -------------------------------------------------------------------------

    catalogue = recurrence[
        [
            "Pleiotropy_Direction_Label",
            "pathway_name",
            "Biological_Theme",
            "Included_in_Figure4B",
            "Exclusion_Reason",
            "N_Diseases",
            "Diseases_Involved",
            "N_Phenotype_Pairs",
            "Mean_N_Overlap",
            "Best_P_Value",
            "Best_FDR",
        ]
    ].copy()

    catalogue = catalogue.rename(
        columns={
            "Pleiotropy_Direction_Label":
                "Pleiotropy Direction",
            "pathway_name":
                "GWAS Catalog Term",
            "Biological_Theme":
                "Biological Theme",
            "Included_in_Figure4B":
                "Included in Figure 4B",
            "Exclusion_Reason":
                "Exclusion Reason",
            "N_Diseases":
                "N Diseases",
            "Diseases_Involved":
                "Diseases Involved",
            "N_Phenotype_Pairs":
                "N Phenotype Pairs",
            "Mean_N_Overlap":
                "Mean N Overlap",
            "Best_P_Value":
                "Best P Value",
            "Best_FDR":
                "Best FDR",
        }
    )

    catalogue["_direction_order"] = (
        catalogue[
            "Pleiotropy Direction"
        ]
        .map(
            {
                "Concordant": 0,
                "Discordant": 1,
            }
        )
    )

    catalogue = (
        catalogue.sort_values(
            [
                "_direction_order",
                "N Diseases",
                "N Phenotype Pairs",
                "GWAS Catalog Term",
            ],
            ascending=[
                True,
                False,
                False,
                True,
            ]
        )
        .drop(
            columns="_direction_order"
        )
        .reset_index(
            drop=True
        )
    )

    catalogue.to_csv(
        OUTPUT_CATALOGUE,
        index=False
    )

    # -------------------------------------------------------------------------
    # OBSERVED THEME SUMMARY
    # -------------------------------------------------------------------------

    theme_long = (
        included.groupby(
            [
                "Assigned_Theme",
                "pleiotropy_direction",
            ],
            as_index=False,
        )
        .agg(
            N_Recurrent_Terms=(
                "pathway_name",
                "nunique"
            ),
            N_Phenotype_Pair_Contributions=(
                "N_Phenotype_Pairs",
                "sum"
            ),
        )
    )

    term_wide = (
        theme_long.pivot(
            index="Assigned_Theme",
            columns="pleiotropy_direction",
            values="N_Recurrent_Terms",
        )
        .fillna(0)
    )

    pair_wide = (
        theme_long.pivot(
            index="Assigned_Theme",
            columns="pleiotropy_direction",
            values="N_Phenotype_Pair_Contributions",
        )
        .fillna(0)
    )

    for direction in [
        "positive",
        "negative",
    ]:

        if direction not in term_wide.columns:
            term_wide[direction] = 0

        if direction not in pair_wide.columns:
            pair_wide[direction] = 0

    positive_total = (
        pair_wide[
            "positive"
        ]
        .sum()
    )

    negative_total = (
        pair_wide[
            "negative"
        ]
        .sum()
    )

    if positive_total != 1163:

        raise RuntimeError(
            "QC FAIL: expected 1,163 concordant "
            "phenotype-pair contributions, "
            f"found {positive_total}"
        )

    if negative_total != 522:

        raise RuntimeError(
            "QC FAIL: expected 522 discordant "
            "phenotype-pair contributions, "
            f"found {negative_total}"
        )

    themes = sorted(
        set(
            term_wide.index
        )
        | set(
            pair_wide.index
        )
    )

    rows = []

    for theme in themes:

        concordant_terms = int(
            term_wide.loc[
                theme,
                "positive"
            ]
        )

        discordant_terms = int(
            term_wide.loc[
                theme,
                "negative"
            ]
        )

        concordant_pairs = int(
            pair_wide.loc[
                theme,
                "positive"
            ]
        )

        discordant_pairs = int(
            pair_wide.loc[
                theme,
                "negative"
            ]
        )

        concordant_percent = (
            100
            * concordant_pairs
            / positive_total
        )

        discordant_percent = (
            100
            * discordant_pairs
            / negative_total
        )

        delta = (
            concordant_percent
            - discordant_percent
        )

        fold = (
            concordant_percent
            / discordant_percent
            if discordant_percent > 0
            else np.inf
        )

        rows.append(
            {
                "Biological Theme":
                    theme,

                "Concordant Recurrent Terms":
                    concordant_terms,

                "Concordant Phenotype-Pair Contributions":
                    concordant_pairs,

                "Concordant (%)":
                    concordant_percent,

                "Discordant Recurrent Terms":
                    discordant_terms,

                "Discordant Phenotype-Pair Contributions":
                    discordant_pairs,

                "Discordant (%)":
                    discordant_percent,

                "Difference (Concordant - Discordant, pp)":
                    delta,

                "Concordant / Discordant Fold Change":
                    fold,
            }
        )

    summary = pd.DataFrame(
        rows
    )

    if len(summary) != 15:

        raise RuntimeError(
            "QC FAIL: expected 15 themes, "
            f"found {len(summary)}"
        )

    summary = (
        summary.sort_values(
            "Difference (Concordant - Discordant, pp)",
            ascending=False,
        )
        .reset_index(
            drop=True
        )
    )

    summary.to_csv(
        OUTPUT_SUMMARY,
        index=False
    )

    # -------------------------------------------------------------------------
    # REPORT
    # -------------------------------------------------------------------------

    print(
        "\nObserved recurrent catalogue:"
    )

    print(
        f"  direction-specific records: "
        f"{len(recurrence)}"
    )

    print(
        f"  unique recurrent terms: "
        f"{n_unique_terms}"
    )

    print(
        f"  included records: "
        f"{len(included)}"
    )

    print(
        f"  included unique terms: "
        f"{included['pathway_name'].nunique()}"
    )

    print(
        f"  excluded records: "
        f"{len(excluded)}"
    )

    print(
        f"  excluded unique terms: "
        f"{excluded['pathway_name'].nunique()}"
    )

    print(
        "\nObserved thematic denominators:"
    )

    print(
        f"  Concordant: "
        f"{int(positive_total)}"
    )

    print(
        f"  Discordant: "
        f"{int(negative_total)}"
    )

    print(
        "\nTheme summary:\n"
    )

    print(
        summary.to_string(
            index=False,
            float_format=lambda x:
                f"{x:.6g}",
        )
    )

    print(
        f"\nSaved recurrent catalogue:\n"
        f"{OUTPUT_CATALOGUE}"
    )

    print(
        f"\nSaved observed theme summary:\n"
        f"{OUTPUT_SUMMARY}"
    )


if __name__ == "__main__":

    main()
