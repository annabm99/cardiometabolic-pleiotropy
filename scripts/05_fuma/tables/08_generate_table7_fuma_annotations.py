#!/usr/bin/env python3

"""
Generate Table 7

Summary of variant- and gene-level FUMA annotations
for concordant and discordant pleiotropic loci.

Outputs:
    - Table7_FUMA_AnnotationSummary.csv
    - Table7_FUMA_AnnotationSummary.xlsx
    - Table7_FUMA_AnnotationSummary.tex
"""

import os
import glob

import numpy as np
import pandas as pd

# ==========================================================
# CONFIGURATION
# ==========================================================

PROJECT_DIR = os.environ.get(
    "CVP_PROJECT_DIR",
    "/path/to/cardiovascular_pleiotropies"
)

# Root directory containing the FUMA outputs
FUMA_DIR = os.environ.get(
    "CVP_FUMA_PLEIOTROPIC_OUT_DIR",
    f"{PROJECT_DIR}/4-FUMA/5-FUMAOut/a-Pleiotropic"
)

# Output directory
OUTPUT_DIR = os.environ.get(
    "CVP_FINAL_TABLES_DIR",
    f"{PROJECT_DIR}/FinalTables"
)

OUTPUT_NAME = "Table7_FUMA_AnnotationSummary"

# Final downstream convergence/FUMA analysis excludes cholesterol traits.
EXCLUDED_TRAITS = {"HDL_t", "LDL_t"}

# Files
LEAD_FILE = "leadSNPs.txt"
GENE_FILE = "genes.txt"

# Columns
FUNC_COL = "func"
CADD_COL = "CADD"
RDB_COL = "RDB"

GENE_ID_COL = "ensg"

PLI_COL = "pLI"
NCRVIS_COL = "ncRVIS"

def get_pair_directories(root_dir):
    """
    Return all phenotype-pair directories.

    Excludes:
        - Dis-Dis
        - txt files
        - hidden folders
    """

    pair_dirs = []

    for entry in sorted(os.listdir(root_dir)):

        full = os.path.join(root_dir, entry)

        if not os.path.isdir(full):
            continue

        if entry == "Dis-Dis":
            continue

        if any(trait in entry for trait in EXCLUDED_TRAITS):
            continue

        pair_dirs.append(full)

    return pair_dirs

# ==========================================================
# LOAD FUMA DATA
# ==========================================================

def load_variant_data(pair_dirs):
    """
    Load annotated lead SNPs from all PosGene/NegGene folders.

    We use:
        leadSNPs.txt -> obtain lead rsIDs
        snps.txt     -> obtain func/CADD/RDB annotations

    Returns
    -------
    positive_snps : DataFrame
    negative_snps : DataFrame
    """

    pos_list = []
    neg_list = []

    for pair_dir in pair_dirs:

        pair = os.path.basename(pair_dir)

        print(f"Processing {pair}")

        for label, container in [("Pos", pos_list), ("Neg", neg_list)]:

            gene_dir = os.path.join(pair_dir, f"{pair}-{label}Gene")

            lead_file = os.path.join(gene_dir, "leadSNPs.txt")
            snp_file = os.path.join(gene_dir, "snps.txt")

            if not (os.path.exists(lead_file) and os.path.exists(snp_file)):
                continue

            lead = pd.read_csv(lead_file, sep="\t")
            snps = pd.read_csv(snp_file, sep="\t")

            # Keep only annotated lead SNPs
            snps = snps[snps["rsID"].isin(lead["rsID"])].copy()

            snps["Pair"] = pair
            snps["Class"] = "Concordant" if label == "Pos" else "Discordant"

            container.append(snps)

    positive_snps = pd.concat(pos_list, ignore_index=True)
    negative_snps = pd.concat(neg_list, ignore_index=True)

    # Global deduplication
    positive_snps = positive_snps.drop_duplicates(subset="rsID")
    negative_snps = negative_snps.drop_duplicates(subset="rsID")

    return positive_snps, negative_snps


# ==========================================================
# LOAD GENE DATA
# ==========================================================

def load_gene_data(pair_dirs):
    """
    Load all mapped genes from FUMA genes.txt.
    """

    pos_list = []
    neg_list = []

    for pair_dir in pair_dirs:

        pair = os.path.basename(pair_dir)

        for label, container in [("Pos", pos_list), ("Neg", neg_list)]:

            gene_file = os.path.join(
                pair_dir,
                f"{pair}-{label}Gene",
                "genes.txt"
            )

            if not os.path.exists(gene_file):
                continue

            genes = pd.read_csv(gene_file, sep="\t")

            genes["Pair"] = pair
            genes["Class"] = "Concordant" if label == "Pos" else "Discordant"

            container.append(genes)

    positive_genes = pd.concat(pos_list, ignore_index=True)
    negative_genes = pd.concat(neg_list, ignore_index=True)

    # Unique genes
    positive_genes = (
        positive_genes
        .sort_values("symbol")
        .drop_duplicates(subset="ensg")
    )

    negative_genes = (
        negative_genes
        .sort_values("symbol")
        .drop_duplicates(subset="ensg")
    )

    return positive_genes, negative_genes

# ==========================================================
# SUMMARY FUNCTIONS
# ==========================================================

def summarize_functional_annotations(pos_snps, neg_snps):
    """
    Summarize functional consequence categories among unique lead SNPs.

    Percentages are calculated relative to the total number of unique
    lead SNPs in each directional class.
    """

    annotations = sorted(
        set(pos_snps[FUNC_COL].dropna()) |
        set(neg_snps[FUNC_COL].dropna())
    )

    total_pos = len(pos_snps)
    total_neg = len(neg_snps)

    rows = []

    for annotation in annotations:
        n_pos = (pos_snps[FUNC_COL] == annotation).sum()
        n_neg = (neg_snps[FUNC_COL] == annotation).sum()

        rows.append({
            "Layer": "Variant-level",
            "Annotation": annotation,
            "Statistic": "Count (% of unique lead SNPs)",
            "Concordant": f"{n_pos} ({100*n_pos/total_pos:.1f}%)",
            "Discordant": f"{n_neg} ({100*n_neg/total_neg:.1f}%)"
        })

    return pd.DataFrame(rows)


def summarize_numeric_metric(
    pos_df,
    neg_df,
    column,
    label,
    layer,
    digits=3
):
    """
    Summarize a numeric annotation separately for concordant and
    discordant entities.

    Reports total N, available N, missingness, mean, median and IQR.
    """

    pos = pd.to_numeric(pos_df[column], errors="coerce")
    neg = pd.to_numeric(neg_df[column], errors="coerce")

    def stats(x):
        n_total = len(x)
        n_available = x.notna().sum()
        n_missing = x.isna().sum()
        missing_pct = 100 * n_missing / n_total if n_total else np.nan

        return {
            "N total": n_total,
            "N available": n_available,
            "Missing, n (%)": (
                f"{n_missing} ({missing_pct:.1f}%)"
            ),
            "Mean": round(x.mean(), digits),
            "Median": round(x.median(), digits),
            "IQR": (
                f"{x.quantile(0.25):.{digits}f}–"
                f"{x.quantile(0.75):.{digits}f}"
            )
        }

    pos_stats = stats(pos)
    neg_stats = stats(neg)

    rows = []

    for statistic in [
        "N total",
        "N available",
        "Missing, n (%)",
        "Mean",
        "Median",
        "IQR"
    ]:
        rows.append({
            "Layer": layer,
            "Annotation": label,
            "Statistic": statistic,
            "Concordant": pos_stats[statistic],
            "Discordant": neg_stats[statistic]
        })

    return pd.DataFrame(rows)


# ==========================================================
# BUILD TABLE 7
# ==========================================================

def build_table7(
    pos_snps,
    neg_snps,
    pos_genes,
    neg_genes
):
    """
    Build a layered annotation summary.

    Variant-level:
        - functional consequence distribution
        - CADD
        - RegulomeDB

    Mapped-gene-level:
        - pLI
        - ncRVIS

    Variant and gene entities are deduplicated within each
    direction before reaching this function.
    """

    # ------------------------------------------------------
    # Variant-level functional consequence composition
    # ------------------------------------------------------

    func_table = summarize_functional_annotations(
        pos_snps,
        neg_snps
    )

    # Sort consequence categories by combined abundance.
    func_table["_Total"] = (
        func_table["Concordant"]
        .str.extract(r"(^\d+)")[0]
        .astype(int)
        +
        func_table["Discordant"]
        .str.extract(r"(^\d+)")[0]
        .astype(int)
    )

    func_table = (
        func_table
        .sort_values("_Total", ascending=False)
        .drop(columns="_Total")
        .reset_index(drop=True)
    )

    # ------------------------------------------------------
    # Variant-level numeric annotations
    # ------------------------------------------------------

    cadd_table = summarize_numeric_metric(
        pos_snps,
        neg_snps,
        CADD_COL,
        "CADD score",
        "Variant-level"
    )

    rdb_table = summarize_numeric_metric(
        pos_snps,
        neg_snps,
        RDB_COL,
        "RegulomeDB score",
        "Variant-level"
    )

    # ------------------------------------------------------
    # Mapped-gene-level constraint annotations
    # ------------------------------------------------------

    pli_table = summarize_numeric_metric(
        pos_genes,
        neg_genes,
        PLI_COL,
        "pLI",
        "Mapped-gene-level"
    )

    ncrvis_table = summarize_numeric_metric(
        pos_genes,
        neg_genes,
        NCRVIS_COL,
        "ncRVIS",
        "Mapped-gene-level"
    )

    final_table = pd.concat(
        [
            func_table,
            cadd_table,
            rdb_table,
            pli_table,
            ncrvis_table
        ],
        ignore_index=True
    )

    return final_table

# ==========================================================
# EXPORT
# ==========================================================

def export_table(table):

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    csv_file = os.path.join(
        OUTPUT_DIR,
        f"{OUTPUT_NAME}.csv"
    )

    xlsx_file = os.path.join(
        OUTPUT_DIR,
        f"{OUTPUT_NAME}.xlsx"
    )

    tex_file = os.path.join(
        OUTPUT_DIR,
        f"{OUTPUT_NAME}.tex"
    )

    table.to_csv(
        csv_file,
        index=False
    )

    try:
        table.to_excel(
            xlsx_file,
            index=False
        )
    except ImportError:
        print(
            "\nWARNING: Excel output was not generated because "
            "openpyxl is not installed in the current Python environment."
        )

    table.to_latex(
        tex_file,
        index=False,
        escape=False
    )

    print("\nTable export completed.")

    print(f"CSV:   {csv_file}")

    if os.path.exists(xlsx_file):
        print(f"Excel: {xlsx_file}")
    else:
        print("Excel: not generated")

    print(f"LaTeX: {tex_file}")

# ==========================================================
# MAIN
# ==========================================================

def main():

    print("\n======================================")
    print("Generating Table 7")
    print("======================================")

    pair_dirs = get_pair_directories(FUMA_DIR)

    print(f"\nDetected {len(pair_dirs)} phenotype pairs.")

    print("\nLoading variant annotations...")

    pos_snps, neg_snps = load_variant_data(pair_dirs)

    print(
        f"  Concordant lead SNPs : {len(pos_snps)}"
    )
    print(
        f"  Discordant lead SNPs : {len(neg_snps)}"
    )

    print("\nLoading gene annotations...")

    pos_genes, neg_genes = load_gene_data(pair_dirs)

    print(
        f"  Concordant genes : {len(pos_genes)}"
    )
    print(
        f"  Discordant genes : {len(neg_genes)}"
    )

    print("\nBuilding summary table...")

    table7 = build_table7(
        pos_snps,
        neg_snps,
        pos_genes,
        neg_genes
    )

    print(table7)

    export_table(table7)

    print("\nDone.")


if __name__ == "__main__":
    main()
