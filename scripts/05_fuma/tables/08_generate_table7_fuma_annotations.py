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

    annotations = sorted(
        set(pos_snps["func"].dropna()) |
        set(neg_snps["func"].dropna())
    )

    total_pos = len(pos_snps)
    total_neg = len(neg_snps)

    rows = []

    for annotation in annotations:

        n_pos = (pos_snps["func"] == annotation).sum()
        n_neg = (neg_snps["func"] == annotation).sum()

        rows.append({
            "Annotation": annotation,
            "Concordant":
                f"{n_pos} ({100*n_pos/total_pos:.1f}%)",
            "Discordant":
                f"{n_neg} ({100*n_neg/total_neg:.1f}%)"
        })

    return pd.DataFrame(rows)


def summarize_numeric_metric(pos_df,
                             neg_df,
                             column,
                             label,
                             digits=2):

    pos = pd.to_numeric(pos_df[column], errors="coerce")
    neg = pd.to_numeric(neg_df[column], errors="coerce")

    return pd.DataFrame([{
        "Annotation": label,
        "Concordant": round(pos.mean(skipna=True), digits),
        "Discordant": round(neg.mean(skipna=True), digits)
    }])

# ==========================================================
# BUILD TABLE 7
# ==========================================================

def build_table7(pos_snps,
                 neg_snps,
                 pos_genes,
                 neg_genes):

    # ------------------------------------------------------
    # Functional annotation summary
    # ------------------------------------------------------

    func_table = summarize_functional_annotations(
        pos_snps,
        neg_snps
    )

    # Sort by total abundance
    func_table["Total"] = (
        func_table["Concordant"].str.extract(r"(^\d+)").astype(int)
        + func_table["Discordant"].str.extract(r"(^\d+)").astype(int)
    )

    func_table = (
        func_table
        .sort_values("Total", ascending=False)
        .drop(columns="Total")
        .reset_index(drop=True)
    )

    # ------------------------------------------------------
    # Numeric summaries
    # ------------------------------------------------------

    metric_tables = [

        summarize_numeric_metric(
            pos_snps,
            neg_snps,
            "CADD",
            "Average CADD score"
        ),

        summarize_numeric_metric(
            pos_snps,
            neg_snps,
            "RDB",
            "Average RegulomeDB score"
        ),

        summarize_numeric_metric(
            pos_genes,
            neg_genes,
            "pLI",
            "Average pLI"
        ),

        summarize_numeric_metric(
            pos_genes,
            neg_genes,
            "ncRVIS",
            "Average ncRVIS"
        )

    ]

    final_table = pd.concat(
        [func_table] + metric_tables,
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

    table.to_excel(
        xlsx_file,
        index=False
    )

    table.to_latex(
        tex_file,
        index=False,
        escape=False
    )

    print("\nTable exported successfully.")
    print(csv_file)
    print(xlsx_file)
    print(tex_file)


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
