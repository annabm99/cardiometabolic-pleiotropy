#!/usr/bin/env python3

"""
Generate Supplementary Table S7.

S7A
---
Recurrent GWAS Catalog trait enrichments across at least three diseases,
separately for concordant and discordant pleiotropy.

All recurrent terms are retained for transparency, including eight
ambiguous/cross-domain terms that are excluded from Figure 4B thematic
aggregation.

S7B
---
Observed thematic composition of recurrent GWAS Catalog enrichments,
together with pair-level permutation inference.

Important
---------
The thematic categories summarize GWAS Catalog phenotype/trait
annotations. They are not molecular pathway categories.

Positive = concordant pleiotropy.
Negative = discordant pleiotropy.

Inputs
------
FinalTables/Table10_Recurrent_GWASCatalog_Terms.csv
FinalTables/Table10_ThemeSummary_Observed.csv
FinalTables/Table10_ThemePermutation.csv

Outputs
-------
FinalTables/Supplement/
    Supplementary_S7A_RecurrentGWASCatalogTerms.csv
    Supplementary_S7B_GWASCatalogThemeProfiles.csv
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# PATHS
# =============================================================================

PROJECT_DIR = Path(
    os.environ.get(
        "CVP_PROJECT_DIR",
        "/path/to/cardiovascular_pleiotropies"
    )
)

FINAL_TABLES_DIR = Path(
    os.environ.get(
        "CVP_FINAL_TABLES_DIR",
        str(PROJECT_DIR / "FinalTables")
    )
)

OUT_DIR = Path(
    os.environ.get(
        "CVP_SUPPLEMENT_TABLES_DIR",
        str(FINAL_TABLES_DIR / "Supplement")
    )
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

TERMS_FILE = (
    FINAL_TABLES_DIR
    / "Table10_Recurrent_GWASCatalog_Terms.csv"
)

OBSERVED_FILE = (
    FINAL_TABLES_DIR
    / "Table10_ThemeSummary_Observed.csv"
)

PERMUTATION_FILE = (
    FINAL_TABLES_DIR
    / "Table10_ThemePermutation.csv"
)

OUTPUT_S7A = (
    OUT_DIR
    / "Supplementary_S7A_RecurrentGWASCatalogTerms.csv"
)

OUTPUT_S7B = (
    OUT_DIR
    / "Supplementary_S7B_GWASCatalogThemeProfiles.csv"
)


# =============================================================================
# LOAD
# =============================================================================

terms = pd.read_csv(
    TERMS_FILE
)

observed = pd.read_csv(
    OBSERVED_FILE
)

permutation = pd.read_csv(
    PERMUTATION_FILE
)


# =============================================================================
# S7A — RECURRENT GWAS CATALOG TERMS
# =============================================================================

required_terms = [
    "Pleiotropy Direction",
    "GWAS Catalog Term",
    "Biological Theme",
    "Included in Figure 4B",
    "Exclusion Reason",
    "N Diseases",
    "Diseases Involved",
    "N Phenotype Pairs",
    "Mean N Overlap",
    "Best P Value",
    "Best FDR",
]

missing = [
    col
    for col in required_terms
    if col not in terms.columns
]

if missing:

    raise ValueError(
        "Missing required recurrent-term columns: "
        + ", ".join(missing)
    )

s7a = terms[
    required_terms
].copy()

s7a = s7a.rename(
    columns={
        "Biological Theme":
            "GWAS Catalog Trait Theme"
    }
)

# Publication-facing rounding only.
# The canonical production table retains the full-precision mean.
s7a[
    "Mean N Overlap"
] = (
    s7a[
        "Mean N Overlap"
    ]
    .round(2)
)

# -------------------------------------------------------------------------
# S7A QC
# -------------------------------------------------------------------------

assert len(s7a) == 234, (
    f"Expected 234 recurrent direction-specific records, "
    f"found {len(s7a)}"
)

assert (
    s7a["GWAS Catalog Term"]
    .nunique()
    == 178
), (
    "Expected 178 unique recurrent GWAS Catalog terms"
)

included = s7a[
    s7a[
        "Included in Figure 4B"
    ]
    == "Yes"
].copy()

excluded = s7a[
    s7a[
        "Included in Figure 4B"
    ]
    == "No"
].copy()

assert len(included) == 224, (
    f"Expected 224 included records, "
    f"found {len(included)}"
)

assert (
    included[
        "GWAS Catalog Term"
    ]
    .nunique()
    == 170
), (
    "Expected 170 unique included recurrent terms"
)

assert len(excluded) == 10, (
    f"Expected 10 excluded records, "
    f"found {len(excluded)}"
)

assert (
    excluded[
        "GWAS Catalog Term"
    ]
    .nunique()
    == 8
), (
    "Expected 8 unique ambiguous excluded terms"
)

assert (
    excluded[
        "Exclusion Reason"
    ]
    .fillna("")
    .str.strip()
    .ne("")
    .all()
), (
    "Excluded S7A records must have an explicit exclusion reason"
)

assert set(
    s7a[
        "Pleiotropy Direction"
    ]
    .unique()
) == {
    "Concordant",
    "Discordant",
}, (
    "Unexpected pleiotropy direction labels in S7A"
)


# =============================================================================
# S7B — THEMATIC PROFILE + STRICT PERMUTATION INFERENCE
# =============================================================================

required_observed = [
    "Biological Theme",
    "Concordant Recurrent Terms",
    "Concordant Phenotype-Pair Contributions",
    "Concordant (%)",
    "Discordant Recurrent Terms",
    "Discordant Phenotype-Pair Contributions",
    "Discordant (%)",
    "Difference (Concordant - Discordant, pp)",
    "Concordant / Discordant Fold Change",
]

missing = [
    col
    for col in required_observed
    if col not in observed.columns
]

if missing:

    raise ValueError(
        "Missing required observed-theme columns: "
        + ", ".join(missing)
    )

required_permutation = [
    "Theme",
    "Positive_Percent",
    "Negative_Percent",
    "Delta_Percent",
    "Fold_Change",
    "Permutation_P",
    "Permutation_FDR",
]

missing = [
    col
    for col in required_permutation
    if col not in permutation.columns
]

if missing:

    raise ValueError(
        "Missing required permutation columns: "
        + ", ".join(missing)
    )

perm = permutation[
    required_permutation
].copy()

perm = perm.rename(
    columns={
        "Theme":
            "Biological Theme",
        "Permutation_P":
            "Permutation P",
        "Permutation_FDR":
            "Permutation FDR",
    }
)

s7b = observed.merge(
    perm[
        [
            "Biological Theme",
            "Positive_Percent",
            "Negative_Percent",
            "Delta_Percent",
            "Fold_Change",
            "Permutation P",
            "Permutation FDR",
        ]
    ],
    on="Biological Theme",
    how="left",
    validate="one_to_one",
)


# =============================================================================
# CROSS-OUTPUT QC
# =============================================================================

assert len(s7b) == 15, (
    f"Expected 15 thematic profiles, "
    f"found {len(s7b)}"
)

assert s7b[
    "Permutation P"
].notna().all(), (
    "Missing permutation P values in S7B"
)

checks = {

    "Concordant percentage": np.allclose(
        s7b["Concordant (%)"],
        s7b["Positive_Percent"],
        rtol=0,
        atol=1e-12,
    ),

    "Discordant percentage": np.allclose(
        s7b["Discordant (%)"],
        s7b["Negative_Percent"],
        rtol=0,
        atol=1e-12,
    ),

    "Percentage-point difference": np.allclose(
        s7b[
            "Difference (Concordant - Discordant, pp)"
        ],
        s7b["Delta_Percent"],
        rtol=0,
        atol=1e-12,
    ),

    "Fold change": np.allclose(
        s7b[
            "Concordant / Discordant Fold Change"
        ],
        s7b["Fold_Change"],
        rtol=0,
        atol=1e-12,
    ),
}

failed = [
    name
    for name, ok in checks.items()
    if not ok
]

if failed:

    raise AssertionError(
        "Observed/permutation mismatch: "
        + ", ".join(failed)
    )

assert (
    s7b[
        "Concordant Phenotype-Pair Contributions"
    ]
    .sum()
    == 1163
), (
    "Expected 1,163 concordant phenotype-pair contributions"
)

assert (
    s7b[
        "Discordant Phenotype-Pair Contributions"
    ]
    .sum()
    == 522
), (
    "Expected 522 discordant phenotype-pair contributions"
)

assert (
    s7b[
        "Permutation FDR"
    ]
    .lt(0.05)
    .sum()
    == 0
), (
    "Unexpected FDR-significant theme detected"
)


# =============================================================================
# CLEAN PUBLICATION-FACING S7B
# =============================================================================

s7b = s7b[
    [
        "Biological Theme",
        "Concordant Recurrent Terms",
        "Concordant Phenotype-Pair Contributions",
        "Concordant (%)",
        "Discordant Recurrent Terms",
        "Discordant Phenotype-Pair Contributions",
        "Discordant (%)",
        "Difference (Concordant - Discordant, pp)",
        "Concordant / Discordant Fold Change",
        "Permutation P",
        "Permutation FDR",
    ]
].copy()

s7b = s7b.rename(
    columns={
        "Biological Theme":
            "GWAS Catalog Trait Theme"
    }
)

# Order by magnitude of observed concordant-discordant contrast.
s7b["_abs_difference"] = (
    s7b[
        "Difference (Concordant - Discordant, pp)"
    ]
    .abs()
)

s7b = (
    s7b.sort_values(
        "_abs_difference",
        ascending=False
    )
    .drop(
        columns="_abs_difference"
    )
    .reset_index(
        drop=True
    )
)


# =============================================================================
# EXPORT
# =============================================================================

s7a.to_csv(
    OUTPUT_S7A,
    index=False
)

s7b.to_csv(
    OUTPUT_S7B,
    index=False
)


# =============================================================================
# REPORT
# =============================================================================

print("=" * 70)
print("SUPPLEMENTARY TABLE S7 QC")
print("=" * 70)

print("\nS7A — recurrent GWAS Catalog terms")
print(f"Rows: {len(s7a)}")
print(
    "Unique terms:",
    s7a[
        "GWAS Catalog Term"
    ].nunique()
)
print(
    "Included in Figure 4B:",
    len(included)
)
print(
    "Excluded ambiguous records:",
    len(excluded)
)

print("\nS7B — GWAS Catalog trait-enrichment profiles")
print(f"Themes: {len(s7b)}")
print(
    "Concordant contributions:",
    s7b[
        "Concordant Phenotype-Pair Contributions"
    ].sum()
)
print(
    "Discordant contributions:",
    s7b[
        "Discordant Phenotype-Pair Contributions"
    ].sum()
)
print(
    "Permutation FDR < 0.05:",
    int(
        s7b[
            "Permutation FDR"
        ]
        .lt(0.05)
        .sum()
    )
)

print("\nCross-output checks:")
for name, ok in checks.items():
    print(
        f"  {name}: "
        f"{'PASS' if ok else 'FAIL'}"
    )

print("\nQC PASSED")
print(f"\nS7A:\n{OUTPUT_S7A}")
print(f"\nS7B:\n{OUTPUT_S7B}")
