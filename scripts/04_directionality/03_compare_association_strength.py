#!/usr/bin/env python3

"""
Compare disease-specific association strength (|Z|) between concordant
and discordant pleiotropic signals after disease-level LD pruning.

Rationale
---------
PleioFDR lead SNPs are LD-pruned within individual phenotype pairs.
After lead SNPs from multiple phenotype pairs are pooled for a disease,
different lead SNPs from different pairs may still be correlated.

This script therefore performs an additional disease-level LD-pruning
step across the pooled concordant and discordant SNPs before statistical
comparison.

Representative SNPs are selected by minimum conjFDR among LD-linked
candidate SNPs. Concordant and discordant SNPs are pruned jointly.

Primary analysis:
    all disease-level observations using the established classification:
        concordant = no discordant pairwise associations
        discordant = >=1 discordant pairwise association

Sensitivity:
    mixed-context discordant observations are removed before LD pruning,
    and pruning/statistical testing is repeated from scratch.

The script requires PLINK 1.9 on PATH.
"""

import glob
import os
import re
import subprocess
import sys
from collections import defaultdict

import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


DISEASES = {
    "CAD_d": {
        "prefix": "CAD",
        "label": "Coronary Artery Disease",
    },
    "HT_d": {
        "prefix": "HT",
        "label": "Hypertension",
    },
    "STR_d": {
        "prefix": "STR",
        "label": "Stroke",
    },
    "T2D_d": {
        "prefix": "T2D",
        "label": "Type 2 Diabetes",
    },
}

EXCLUDED_TRAITS = {"HDL_t", "LDL_t"}

LD_R2_THRESHOLD = 0.1


def usage():
    print(
        "Usage:\n"
        "  python 03_compare_association_strength.py "
        "<aggregate_dir> <pairwise_dir> <plink_reference_dir> <out_dir>\n\n"
        "Arguments:\n"
        "  aggregate_dir        Output of 02_aggregate_no_hdl_ldl.py\n"
        "  pairwise_dir         Pairwise positive/negative a-IndepLead files\n"
        "  plink_reference_dir  Directory containing chr1-22 .bed/.bim/.fam\n"
        "  out_dir              Output directory\n"
    )


def load_baseline(aggregate_dir):
    """
    Load the exact disease-level observations used by the existing
    directionality analysis.

    These files remain authoritative for disease-specific Z-scores and
    concordant/discordant classification.
    """

    frames = []

    for disease, info in DISEASES.items():

        prefix = info["prefix"]

        pos_path = os.path.join(
            aggregate_dir,
            f"{prefix}_positive_snps_noChol.csv.gz",
        )

        neg_path = os.path.join(
            aggregate_dir,
            f"{prefix}_negative_snps_noChol.csv.gz",
        )

        if not os.path.exists(pos_path):
            raise FileNotFoundError(pos_path)

        if not os.path.exists(neg_path):
            raise FileNotFoundError(neg_path)

        pos = pd.read_csv(
            pos_path,
            sep="\t",
            compression="gzip",
        )

        neg = pd.read_csv(
            neg_path,
            sep="\t",
            compression="gzip",
        )

        pos["classification"] = "Positive"
        neg["classification"] = "Negative"

        frames.extend([pos, neg])

    baseline = pd.concat(
        frames,
        ignore_index=True,
    )

    baseline["mixed"] = (
        (baseline["n_positive_pairs"] > 0)
        & (baseline["n_negative_pairs"] > 0)
    )

    # Disease-SNP must already be unique after aggregation.
    duplicates = baseline.duplicated(
        ["Disease", "SNP"],
        keep=False,
    )

    if duplicates.any():
        raise ValueError(
            "Aggregated input contains duplicated disease-SNP observations"
        )

    return baseline


def recover_pairwise_metadata(pairwise_dir):
    """
    Recover chromosome, position and minimum conjFDR for each disease-SNP
    from the pairwise PleioFDR lead-SNP files.
    """

    files = sorted(
        glob.glob(
            os.path.join(
                pairwise_dir,
                "*_positive.csv.gz",
            )
        )
        + glob.glob(
            os.path.join(
                pairwise_dir,
                "*_negative.csv.gz",
            )
        )
    )

    if not files:
        raise FileNotFoundError(
            f"No pairwise positive/negative files found in {pairwise_dir}"
        )

    rows = []

    for path in files:

        df = pd.read_csv(
            path,
            sep="\t",
            compression="gzip",
            low_memory=False,
        )

        z_columns = [
            col
            for col in df.columns
            if col.startswith("Z-")
        ]

        phenotypes = [
            col[2:]
            for col in z_columns
        ]

        if len(phenotypes) != 2:
            raise ValueError(
                f"Expected two Z columns in {path}; found {z_columns}"
            )

        # Match the no-HDL/LDL aggregation exactly.
        if any(
            phenotype in EXCLUDED_TRAITS
            for phenotype in phenotypes
        ):
            continue

        fdr_columns = [
            col
            for col in df.columns
            if col.startswith("conjfdr_")
        ]

        if len(fdr_columns) != 1:
            raise ValueError(
                f"Expected one conjFDR column in {path}; "
                f"found {fdr_columns}"
            )

        fdr_col = fdr_columns[0]

        for disease in DISEASES:

            if disease not in phenotypes:
                continue

            required = [
                "SNP",
                "chrnum",
                "chrpos",
                fdr_col,
            ]

            missing = [
                col
                for col in required
                if col not in df.columns
            ]

            if missing:
                raise ValueError(
                    f"Missing columns in {path}: {missing}"
                )

            tmp = df[required].copy()

            tmp.columns = [
                "SNP",
                "CHR",
                "BP",
                "conjFDR",
            ]

            tmp["Disease"] = disease

            rows.append(tmp)

    if not rows:
        raise ValueError(
            "No disease-containing pairwise rows were recovered"
        )

    meta_raw = pd.concat(
        rows,
        ignore_index=True,
    )

    meta_raw["CHR"] = pd.to_numeric(
        meta_raw["CHR"],
        errors="raise",
    )

    meta_raw["BP"] = pd.to_numeric(
        meta_raw["BP"],
        errors="raise",
    )

    meta_raw["conjFDR"] = pd.to_numeric(
        meta_raw["conjFDR"],
        errors="coerce",
    )

    if meta_raw["conjFDR"].isna().any():
        raise ValueError(
            "Missing/non-numeric conjFDR values found"
        )

    # Position consistency.
    position_qc = (
        meta_raw
        .groupby(["Disease", "SNP"])
        .agg(
            n_chr=("CHR", "nunique"),
            n_bp=("BP", "nunique"),
        )
        .reset_index()
    )

    inconsistent = position_qc[
        (position_qc["n_chr"] != 1)
        | (position_qc["n_bp"] != 1)
    ]

    if not inconsistent.empty:
        raise ValueError(
            "Inconsistent CHR/BP values found for disease-SNP observations:\n"
            + inconsistent.head(20).to_string(index=False)
        )

    metadata = (
        meta_raw
        .sort_values(
            ["conjFDR", "SNP"]
        )
        .groupby(
            ["Disease", "SNP"],
            as_index=False,
        )
        .agg(
            CHR=("CHR", "first"),
            BP=("BP", "first"),
            min_conjFDR=("conjFDR", "min"),
            n_pairwise_rows=("SNP", "size"),
        )
    )

    metadata["CHR"] = metadata["CHR"].astype(int)
    metadata["BP"] = metadata["BP"].astype(int)

    return metadata


def attach_metadata(baseline, metadata):
    """
    Attach LD coordinates and min conjFDR while preserving the existing
    disease-SNP classification and Z-score exactly.
    """

    merged = baseline.merge(
        metadata,
        on=["Disease", "SNP"],
        how="left",
        validate="one_to_one",
        indicator=True,
    )

    missing = merged[
        merged["_merge"] != "both"
    ]

    if not missing.empty:
        raise ValueError(
            "Some aggregated disease-SNP observations lack pairwise "
            "metadata:\n"
            + missing[
                ["Disease", "SNP"]
            ].head(20).to_string(index=False)
        )

    merged = merged.drop(
        columns="_merge"
    )

    required = [
        "CHR",
        "BP",
        "min_conjFDR",
    ]

    if merged[required].isna().any().any():
        raise ValueError(
            "Missing LD-pruning metadata after merge"
        )

    return merged


def run_plink_ld(data, reference_dir, work_dir):
    """
    Calculate LD only among candidate SNPs.

    All same-chromosome candidate pairs are considered; no biological
    distance cutoff is imposed. PLINK reports pairs with R2 >= 0.1.
    """

    os.makedirs(
        work_dir,
        exist_ok=True,
    )

    plink = os.environ.get(
        "PLINK",
        "plink",
    )

    graph = {
        disease: defaultdict(set)
        for disease in DISEASES
    }

    edge_counts = {
        disease: 0
        for disease in DISEASES
    }

    for disease in DISEASES:

        disease_df = data[
            data["Disease"] == disease
        ]

        for chrom in range(1, 23):

            chromosome_df = disease_df[
                disease_df["CHR"] == chrom
            ]

            snps = sorted(
                chromosome_df["SNP"]
                .dropna()
                .unique()
            )

            if not snps:
                continue

            # A single candidate SNP is already independent:
            # there is no pairwise LD to calculate.
            if len(snps) == 1:
                continue

            bfile = os.path.join(
                reference_dir,
                f"chr{chrom}",
            )

            for suffix in [
                ".bed",
                ".bim",
                ".fam",
            ]:
                if not os.path.exists(
                    bfile + suffix
                ):
                    raise FileNotFoundError(
                        bfile + suffix
                    )

            extract_path = os.path.join(
                work_dir,
                f"{disease}.chr{chrom}.snps.txt",
            )

            prefix = os.path.join(
                work_dir,
                f"{disease}.chr{chrom}",
            )

            with open(
                extract_path,
                "w",
            ) as handle:
                for snp in snps:
                    handle.write(
                        f"{snp}\n"
                    )

            command = [
                plink,
                "--bfile",
                bfile,
                "--extract",
                extract_path,
                "--r2",
                "gz",
                "--ld-window",
                "999999",
                "--ld-window-kb",
                "1000000",
                "--ld-window-r2",
                str(LD_R2_THRESHOLD),
                "--out",
                prefix,
            ]

            subprocess.run(
                command,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )

            # Verify all requested SNPs were present.
            log_path = (
                prefix + ".log"
            )

            with open(
                log_path,
                "r",
            ) as handle:
                log_text = handle.read()

            match = re.search(
                r"--extract:\s+(\d+)\s+variants remaining",
                log_text,
            )

            if match is None:
                raise ValueError(
                    f"Could not verify extracted variant count in {log_path}"
                )

            n_extracted = int(
                match.group(1)
            )

            if n_extracted != len(snps):
                raise ValueError(
                    f"{disease} chr{chrom}: requested {len(snps)} SNPs "
                    f"but PLINK retained {n_extracted}"
                )

            ld_path = (
                prefix + ".ld.gz"
            )

            # A chromosome containing a single candidate has no pairwise
            # LD output and requires no pruning.
            if not os.path.exists(
                ld_path
            ):
                if len(snps) == 1:
                    continue

                raise FileNotFoundError(
                    ld_path
                )

            ld = pd.read_csv(
                ld_path,
                sep=r"\s+",
                compression="gzip",
            )

            if ld.empty:
                continue

            ld = ld[
                ld["R2"]
                >= LD_R2_THRESHOLD
            ]

            for snp_a, snp_b in zip(
                ld["SNP_A"],
                ld["SNP_B"],
            ):

                graph[disease][
                    snp_a
                ].add(
                    snp_b
                )

                graph[disease][
                    snp_b
                ].add(
                    snp_a
                )

                edge_counts[
                    disease
                ] += 1

    return graph, edge_counts


def prune_disease(data, graph):
    """
    Greedy LD pruning analogous to the PleioFDR lead-SNP logic:

    1. rank by minimum conjFDR;
    2. retain the strongest remaining SNP;
    3. remove directly LD-linked candidate SNPs;
    4. continue until none remain.

    Concordant and discordant variants are pruned jointly.
    """

    ranked = data.sort_values(
        [
            "min_conjFDR",
            "SNP",
        ],
        ascending=[
            True,
            True,
        ],
    ).copy()

    candidate_snps = set(
        ranked["SNP"]
    )

    removed = set()
    retained = []

    for _, row in ranked.iterrows():

        snp = row["SNP"]

        if snp in removed:
            continue

        retained.append(
            snp
        )

        neighbours = graph.get(
            snp,
            set(),
        )

        removed.update(
            neighbour
            for neighbour in neighbours
            if (
                neighbour in candidate_snps
                and neighbour != snp
            )
        )

    return ranked[
        ranked["SNP"].isin(
            retained
        )
    ].copy()


def compare_abs_z(data):
    """
    Compare disease-specific association strength between concordant and
    discordant LD-independent signals.
    """

    concordant = data.loc[
        data["classification"]
        == "Positive",
        "abs_z",
    ].dropna()

    discordant = data.loc[
        data["classification"]
        == "Negative",
        "abs_z",
    ].dropna()

    if (
        len(concordant) == 0
        or len(discordant) == 0
    ):
        raise ValueError(
            "Both directionality classes are required"
        )

    stat, p_value = (
        mannwhitneyu(
            concordant,
            discordant,
            alternative="two-sided",
        )
    )

    # Positive values mean discordant signals tend to have larger |Z|.
    rank_biserial = (
        1
        - (
            2 * stat
            / (
                len(concordant)
                * len(discordant)
            )
        )
    )

    return {
        "N_concordant":
            len(concordant),
        "N_discordant":
            len(discordant),
        "Concordant_median_abs_z":
            concordant.median(),
        "Discordant_median_abs_z":
            discordant.median(),
        "Delta_median_abs_z_discordant_minus_concordant":
            (
                discordant.median()
                - concordant.median()
            ),
        "Concordant_mean_abs_z":
            concordant.mean(),
        "Discordant_mean_abs_z":
            discordant.mean(),
        "Mann_Whitney_U":
            stat,
        "P_value":
            p_value,
        "Rank_biserial_discordant_gt_concordant":
            rank_biserial,
    }


def run_analysis(data, graph, analysis_name):
    """
    LD-prune independently within each disease and perform the four
    disease-specific tests.
    """

    result_rows = []
    retained_frames = []

    for disease, info in DISEASES.items():

        disease_data = data[
            data["Disease"] == disease
        ].copy()

        retained = prune_disease(
            disease_data,
            graph[disease],
        )

        retained[
            "Analysis"
        ] = analysis_name

        retained_frames.append(
            retained
        )

        stats = compare_abs_z(
            retained
        )

        stats[
            "Analysis"
        ] = analysis_name

        stats[
            "Disease"
        ] = info["label"]

        stats[
            "Disease_ID"
        ] = disease

        result_rows.append(
            stats
        )

    results = pd.DataFrame(
        result_rows
    )

    # Multiple-testing correction across the four disease-specific tests.
    results[
        "P_FDR"
    ] = multipletests(
        results["P_value"],
        method="fdr_bh",
    )[1]

    results[
        "P_Bonferroni"
    ] = multipletests(
        results["P_value"],
        method="bonferroni",
    )[1]

    retained = pd.concat(
        retained_frames,
        ignore_index=True,
    )

    return results, retained


def main():

    if len(sys.argv) != 5:
        usage()
        sys.exit(1)

    aggregate_dir = (
        sys.argv[1]
        .rstrip("/")
    )

    pairwise_dir = (
        sys.argv[2]
        .rstrip("/")
    )

    reference_dir = (
        sys.argv[3]
        .rstrip("/")
    )

    out_dir = (
        sys.argv[4]
        .rstrip("/")
    )

    os.makedirs(
        out_dir,
        exist_ok=True,
    )

    work_dir = os.path.join(
        out_dir,
        "ld_work",
    )

    print(
        "Loading authoritative disease-level observations..."
    )

    baseline = load_baseline(
        aggregate_dir
    )

    print(
        f"Baseline observations: {len(baseline)}"
    )

    print(
        f"Mixed-context observations: {int(baseline['mixed'].sum())}"
    )

    print(
        "\nRecovering pairwise conjFDR/position metadata..."
    )

    metadata = recover_pairwise_metadata(
        pairwise_dir
    )

    prepared = attach_metadata(
        baseline,
        metadata,
    )

    prepared_path = os.path.join(
        out_dir,
        "AssociationStrength_InputWithLDMetadata.csv.gz",
    )

    prepared.to_csv(
        prepared_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    print(
        "\nCalculating disease-level LD..."
    )

    graph, edge_counts = run_plink_ld(
        prepared,
        reference_dir,
        work_dir,
    )

    print(
        "\nLD pairs with R2 >= "
        f"{LD_R2_THRESHOLD}:"
    )

    for disease in DISEASES:
        print(
            f"  {disease}: "
            f"{edge_counts[disease]}"
        )

    # -----------------------------------------------------
    # PRIMARY ANALYSIS
    # -----------------------------------------------------

    primary_stats, primary_snps = (
        run_analysis(
            prepared,
            graph,
            "Primary_LD_pruned",
        )
    )

    # -----------------------------------------------------
    # MIXED-CONTEXT SENSITIVITY
    #
    # Mixed observations are removed BEFORE re-pruning.
    # -----------------------------------------------------

    no_mixed_input = prepared[
        ~prepared["mixed"]
    ].copy()

    sensitivity_stats, sensitivity_snps = (
        run_analysis(
            no_mixed_input,
            graph,
            "Sensitivity_LD_pruned_no_mixed",
        )
    )

    # -----------------------------------------------------
    # OUTPUTS
    # -----------------------------------------------------

    summary = pd.concat(
        [
            primary_stats,
            sensitivity_stats,
        ],
        ignore_index=True,
    )

    summary_path = os.path.join(
        out_dir,
        "AssociationStrength_LDPruned_Summary.csv",
    )

    summary.to_csv(
        summary_path,
        index=False,
    )

    primary_path = os.path.join(
        out_dir,
        "AssociationStrength_LDPruned_Primary_SNPs.csv.gz",
    )

    primary_snps.to_csv(
        primary_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    sensitivity_path = os.path.join(
        out_dir,
        "AssociationStrength_LDPruned_NoMixed_SNPs.csv.gz",
    )

    sensitivity_snps.to_csv(
        sensitivity_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    mixed_cases = prepared[
        prepared["mixed"]
    ].copy()

    mixed_path = os.path.join(
        out_dir,
        "MixedDirectionality_Cases.csv",
    )

    mixed_cases.to_csv(
        mixed_path,
        index=False,
    )

    retention_rows = []

    for disease, info in DISEASES.items():

        original = prepared[
            prepared["Disease"]
            == disease
        ]

        primary = primary_snps[
            primary_snps["Disease"]
            == disease
        ]

        sensitivity = sensitivity_snps[
            sensitivity_snps["Disease"]
            == disease
        ]

        retention_rows.append({
            "Disease":
                info["label"],
            "Original_N":
                len(original),
            "Primary_LD_pruned_N":
                len(primary),
            "No_mixed_LD_pruned_N":
                len(sensitivity),
        })

    retention = pd.DataFrame(
        retention_rows
    )

    retention_path = os.path.join(
        out_dir,
        "AssociationStrength_LDPruning_QC.csv",
    )

    retention.to_csv(
        retention_path,
        index=False,
    )

    print(
        "\n=== ASSOCIATION-STRENGTH RESULTS ==="
    )

    display_columns = [
        "Analysis",
        "Disease",
        "N_concordant",
        "N_discordant",
        "Concordant_median_abs_z",
        "Discordant_median_abs_z",
        "Delta_median_abs_z_discordant_minus_concordant",
        "P_value",
        "P_FDR",
        "P_Bonferroni",
    ]

    print(
        summary[
            display_columns
        ].to_string(
            index=False,
            float_format=lambda value:
                f"{value:.6g}",
        )
    )

    print(
        "\n=== LD-PRUNING RETENTION ==="
    )

    print(
        retention.to_string(
            index=False
        )
    )

    print(
        "\nSaved:"
    )

    for path in [
        summary_path,
        primary_path,
        sensitivity_path,
        mixed_path,
        retention_path,
    ]:
        print(
            f"  {path}"
        )

    print(
        "\nDONE"
    )


if __name__ == "__main__":
    main()
