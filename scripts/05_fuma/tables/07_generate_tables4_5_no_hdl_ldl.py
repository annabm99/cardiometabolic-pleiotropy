#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
from pathlib import Path
import re

##############################################################################
# PATHS
##############################################################################

PROJECT_DIR = Path(
    os.environ.get("CVP_PROJECT_DIR", "/path/to/cardiovascular_pleiotropies")
)

BASE_DIR = Path(
    os.environ.get(
        "CVP_FUMA_PLEIOTROPIC_OUT_DIR",
        str(PROJECT_DIR / "4-FUMA/5-FUMAOut/a-Pleiotropic")
    )
)

OUT_DIR = Path(
    os.environ.get("CVP_FINAL_TABLES_DIR", str(PROJECT_DIR / "FinalTables"))
)

OUT_DIR.mkdir(parents=True, exist_ok=True)

##############################################################################
# HDL / LDL EXCLUSION
##############################################################################

EXCLUDED_TRAITS = {
    "HDL_t",
    "LDL_t"
}

##############################################################################
# UTILITIES
##############################################################################

DISEASES = {"CAD", "HT", "STR", "T2D"}

def extract_diseases(folder_name):

    diseases = set()

    for disease in DISEASES:
        if folder_name.startswith(disease):
            diseases.add(disease)

    return diseases

##############################################################################
# DISCOVER PHENOTYPE FOLDERS
##############################################################################

all_pair_dirs = sorted(
    [x for x in BASE_DIR.iterdir() if x.is_dir()]
)

pair_dirs = []
removed_pairs = []

for pair_dir in all_pair_dirs:

    pair_name = pair_dir.name

    # Dis-Dis is a separate disease-disease analysis and must not
    # enter the disease-trait convergence analysis.
    if pair_name == "Dis-Dis":
        removed_pairs.append(pair_name)
        continue

    if any(trait in pair_name for trait in EXCLUDED_TRAITS):

        removed_pairs.append(pair_name)

    else:

        pair_dirs.append(pair_dir)

print("\n=================================================")
print("HDL / LDL EXCLUSION QC")
print("=================================================\n")

print("Excluded traits:")
for t in sorted(EXCLUDED_TRAITS):
    print(f"  {t}")

print("\nRemoved folders:")
for x in sorted(removed_pairs):
    print(x)

print(f"\nTotal removed folders: {len(removed_pairs)}")
print(f"Remaining folders: {len(pair_dirs)}")

##############################################################################
# LOAD HDL / LDL LOCI FOR LATER ANNOTATION
##############################################################################

chol_loci = []

for pair_dir in all_pair_dirs:

    pair_name = pair_dir.name

    if not any(trait in pair_name for trait in EXCLUDED_TRAITS):
        continue

    for tag, pleio_type in [
        ("PosGene", "Concordant"),
        ("NegGene", "Discordant")
    ]:

        subdirs = list(pair_dir.glob(f"*{tag}"))

        if len(subdirs) == 0:
            continue

        loci_file = subdirs[0] / "GenomicRiskLoci.txt"

        if not loci_file.exists():
            continue

        loci = pd.read_csv(loci_file, sep="\t")

        loci = loci.rename(columns={
            "chr": "Chromosome",
            "start": "Start",
            "end": "End"
        })

        loci["Chromosome"] = loci["Chromosome"].astype(str)

        pretty_name = (
            pair_name
            .replace("_d", "")
            .replace("_t", "")
        )

        loci["Phenotype_Pair"] = pretty_name
        loci["Pleiotropy_Type"] = pleio_type

        chol_loci.append(loci)

chol_loci = pd.concat(chol_loci, ignore_index=True)

print(f"Loaded cholesterol loci: {len(chol_loci):,}")

##############################################################################
# STEP 1
# LOAD FUMA LOCI
##############################################################################

print("\n=================================================")
print("LOADING FUMA LOCI")
print("=================================================\n")

all_loci = []

for pair_dir in pair_dirs:

    pair_name = pair_dir.name

    analyses = [
        ("PosGene", "Concordant"),
        ("NegGene", "Discordant")
    ]

    for tag, pleio_type in analyses:

        subdirs = list(pair_dir.glob(f"*{tag}"))

        if len(subdirs) == 0:
            continue

        fuma_dir = subdirs[0]

        loci_file = fuma_dir / "GenomicRiskLoci.txt"

        if not loci_file.exists():
            continue

        loci = pd.read_csv(
            loci_file,
            sep="\t"
        )

        loci["Phenotype_Pair"] = (
            pair_name
            .replace("_d", "")
            .replace("_t", "")
        )

        loci["Pleiotropy_Type"] = pleio_type

        all_loci.append(loci)

all_loci = pd.concat(
    all_loci,
    ignore_index=True
)

print(f"Input loci: {len(all_loci):,}")

##############################################################################
# STANDARDIZE
##############################################################################

all_loci = all_loci.rename(columns={
    "chr": "Chromosome",
    "start": "Start",
    "end": "End",
    "LeadSNPs": "LeadSNP",
    "p": "P"
})

all_loci["Chromosome"] = all_loci["Chromosome"].astype(str)

all_loci["OriginalLocusID"] = np.arange(
    len(all_loci)
)

##############################################################################
# STEP 2
# MERGE LOCI
##############################################################################

print("\n=================================================")
print("MERGING LOCI")
print("=================================================\n")

merged_rows = []
membership = {}

for pleio_type in ["Concordant", "Discordant"]:

    for chrom in sorted(
        all_loci["Chromosome"].unique()
    ):

        subset = all_loci[
            (all_loci["Chromosome"] == chrom) &
            (all_loci["Pleiotropy_Type"] == pleio_type)
        ].copy()

        if len(subset) == 0:
            continue

        subset = subset.sort_values(
            "Start"
        )

        current = None

        for _, row in subset.iterrows():

            if current is None:

                current = {
                    "start": row["Start"],
                    "end": row["End"],
                    "rows": [row]
                }

                continue

            if row["Start"] <= current["end"]:

                current["end"] = max(
                    current["end"],
                    row["End"]
                )

                current["rows"].append(row)

            else:

                merged_rows.append(
                    (chrom, pleio_type, current)
                )

                current = {
                    "start": row["Start"],
                    "end": row["End"],
                    "rows": [row]
                }

        merged_rows.append(
            (chrom, pleio_type, current)
        )

##############################################################################
# BUILD TABLE 4
##############################################################################

table4 = []

for chrom, pleio_type, block in merged_rows:

    rows = pd.DataFrame(
        block["rows"]
    )

    best_idx = rows["P"].idxmin()
    best_row = rows.loc[best_idx]

    diseases = set()

    for pair in rows["Phenotype_Pair"]:
        diseases |= extract_diseases(pair)

    locus_id = f"ML_{len(table4)+1}"

    for oid in rows["OriginalLocusID"]:
        membership[oid] = locus_id

    table4.append({

        "MergedLocusID":
            locus_id,

        "Chromosome":
            chrom,

        "Merged_Start":
            block["start"],

        "Merged_End":
            block["end"],

        "Representative_Lead_SNP":
            best_row["LeadSNP"],

        "Best_P":
            rows["P"].min(),

        "Pleiotropy_Type":
            pleio_type,

        "Diseases_Involved":
            "|".join(sorted(diseases)),

        "Number_of_Diseases":
            len(diseases),

        "Phenotype_Pairs":
            "|".join(
                sorted(
                    p.replace("_d", "").replace("_t", "")
                    for p in rows["Phenotype_Pair"].unique()
                )
            ),

        "Number_of_Pairs":
            rows["Phenotype_Pair"].nunique(),

        "Num_Original_Loci":
            len(rows)
    })

table4 = pd.DataFrame(table4)

##############################################################################
# ANNOTATE WITH HDL / LDL PAIRS
##############################################################################

chol_concordant = []
chol_discordant = []

for _, locus in table4.iterrows():

    overlaps = chol_loci[
        (chol_loci["Chromosome"] == str(locus["Chromosome"])) &
        (chol_loci["Start"] <= locus["Merged_End"]) &
        (chol_loci["End"] >= locus["Merged_Start"])
    ]

    concordant = sorted(
        overlaps.loc[
            overlaps["Pleiotropy_Type"] == "Concordant",
            "Phenotype_Pair"
        ].unique()
    )

    discordant = sorted(
        overlaps.loc[
            overlaps["Pleiotropy_Type"] == "Discordant",
            "Phenotype_Pair"
        ].unique()
    )

    chol_concordant.append(";".join(concordant))
    chol_discordant.append(";".join(discordant))

table4["Cholesterol_Concordant"] = chol_concordant
table4["Cholesterol_Discordant"] = chol_discordant

##############################################################################
# BUILD LOCUS MEMBERSHIP TABLE
##############################################################################

membership_rows = []

for chrom, pleio_type, block in merged_rows:

    rows = pd.DataFrame(block["rows"])

    locus_id = None

    # recover ML ID assigned earlier
    for oid in rows["OriginalLocusID"]:
        locus_id = membership[oid]
        break

    for _, row in rows.iterrows():

        membership_rows.append({

            "OriginalLocusID":
                row["OriginalLocusID"],

            "Phenotype_Pair":
                row["Phenotype_Pair"],

            "Pleiotropy_Type":
                row["Pleiotropy_Type"],

            "GenomicLocus":
                row["GenomicLocus"],

            "LeadSNP":
                row["LeadSNP"],

            "Chromosome":
                row["Chromosome"],

            "Start":
                row["Start"],

            "End":
                row["End"],

            "MergedLocusID":
                membership[row["OriginalLocusID"]]
        })

table4_membership = pd.DataFrame(
    membership_rows
)

##############################################################################
# QC
##############################################################################

coverage = len(membership)

coverage_pct = (
    coverage /
    len(all_loci)
) * 100

print(f"Input loci: {len(all_loci):,}")
print(f"Covered loci: {coverage:,}")
print(f"Merged loci: {len(table4):,}")
print(f"Coverage: {coverage_pct:.2f}%")

missing = len(all_loci) - coverage

print(f"Missing loci: {missing}")

if missing != 0:
    raise RuntimeError(
        "Coverage not 100%"
    )

print()

print(
    "Concordant merged loci:",
    (table4["Pleiotropy_Type"] == "Concordant").sum()
)

print(
    "Discordant merged loci:",
    (table4["Pleiotropy_Type"] == "Discordant").sum()
)

##############################################################################
# STEP 3
# LOAD GENE EVIDENCE
##############################################################################

print("\n=================================================")
print("LOADING GENE EVIDENCE")
print("=================================================\n")

gene_evidence = []

for pair_dir in pair_dirs:

    pair_name = pair_dir.name

    diseases = extract_diseases(pair_name)

    for tag, pleio_type in [
        ("PosGene", "Concordant"),
        ("NegGene", "Discordant")
    ]:

        subdirs = list(
            pair_dir.glob(f"*{tag}")
        )

        if len(subdirs) == 0:
            continue

        gene_file = (
            subdirs[0] /
            "genes.txt"
        )

        if not gene_file.exists():
            continue

        genes = pd.read_csv(
            gene_file,
            sep="\t"
        )

        keep_cols = [
            "ensg",
            "symbol",
            "chr",
            "start",
            "end",
            "type",
            "GenomicLocus",
            "IndSigSNPs",
            "pLI",
            "ncRVIS"
        ]

        genes = genes[keep_cols].copy()

        genes["Phenotype_Pair"] = (
            pair_name
                .replace("_d", "")
                .replace("_t", "")
        )
        genes["Pleiotropy_Type"] = pleio_type

        genes["Disease"] = "|".join(
            sorted(diseases)
        )

        gene_evidence.append(
            genes
        )

gene_evidence = pd.concat(
    gene_evidence,
    ignore_index=True
)

##############################################################################
# LOAD HDL / LDL GENE EVIDENCE
##############################################################################

chol_gene_evidence = []

for pair_dir in all_pair_dirs:

    pair_name = pair_dir.name

    if not any(trait in pair_name for trait in EXCLUDED_TRAITS):
        continue

    diseases = extract_diseases(pair_name)

    for tag, pleio_type in [
        ("PosGene", "Concordant"),
        ("NegGene", "Discordant")
    ]:

        subdirs = list(pair_dir.glob(f"*{tag}"))

        if len(subdirs) == 0:
            continue

        gene_file = subdirs[0] / "genes.txt"

        if not gene_file.exists():
            continue

        genes = pd.read_csv(gene_file, sep="\t")

        genes = genes[keep_cols].copy()

        genes["Phenotype_Pair"] = (
            pair_name
            .replace("_d", "")
            .replace("_t", "")
        )

        genes["Pleiotropy_Type"] = pleio_type

        genes["Disease"] = "|".join(sorted(diseases))

        chol_gene_evidence.append(genes)

chol_gene_evidence = pd.concat(
    chol_gene_evidence,
    ignore_index=True
)

##############################################################################
# MAP GENE EVIDENCE TO MERGED LOCI
##############################################################################

gene_evidence["GenomicLocus"] = (
    gene_evidence["GenomicLocus"]
    .astype(str)
)

# Some FUMA genes are assigned to more than one genomic risk locus
# (e.g. "63:64"). Expand these into one gene-locus observation per
# constituent FUMA locus before mapping to merged loci.

gene_evidence["GenomicLocus"] = (
    gene_evidence["GenomicLocus"]
    .str.split(":")
)

gene_evidence = gene_evidence.explode(
    "GenomicLocus",
    ignore_index=True
)

gene_evidence["GenomicLocus"] = (
    gene_evidence["GenomicLocus"]
    .str.strip()
)

table4_membership["GenomicLocus"] = (
    table4_membership["GenomicLocus"]
    .astype(str)
)


gene_evidence = gene_evidence.merge(
    table4_membership[
        [
            "Phenotype_Pair",
            "Pleiotropy_Type",
            "GenomicLocus",
            "MergedLocusID"
        ]
    ],
    on=[
        "Phenotype_Pair",
        "Pleiotropy_Type",
        "GenomicLocus"
    ],
    how="left"
)

##############################################################################
# BUILD LOCUS → MAPPED GENES SUMMARY
##############################################################################

locus_gene_map = {}

for ml, g in gene_evidence.groupby("MergedLocusID"):

    if pd.isna(ml):
        continue

    counts = (
        g["symbol"]
        .value_counts()
    )

    ordered_genes = sorted(
        counts.index,
        key=lambda x: (-counts[x], x)
    )

    locus_gene_map[ml] = ";".join(ordered_genes)

table4["Mapped_Genes"] = (
    table4["MergedLocusID"]
    .map(locus_gene_map)
    .fillna("")
)

##############################################################################
# QC
##############################################################################

missing_ml = (
    gene_evidence["MergedLocusID"]
    .isna()
    .sum()
)

print(
    f"Gene observations without ML mapping: "
    f"{missing_ml:,}"
)

##############################################################################
# KEEP ONLY GENES PRESENT IN ALL FOUR DISEASES
##############################################################################

gene_summary = []

for ensg, g in gene_evidence.groupby("ensg"):

    diseases = set()

    for d in g["Disease"]:
        diseases |= {
            x for x in d.split("|")
            if x
        }

    symbol = (
        g["symbol"].dropna().iloc[0]
        if g["symbol"].notna().any()
        else ensg
    )

    gene_summary.append({

        "symbol":
            symbol,

        "Ensembl_ID":
            ensg,

        "Gene_Type":
            g["type"].iloc[0],

        "Chromosome":
            g["chr"].iloc[0],

        "Gene_Start":
            g["start"].min(),

        "Gene_End":
            g["end"].max(),

        "pLI":
            g["pLI"].dropna().iloc[0]
            if g["pLI"].notna().any()
            else np.nan,

        "ncRVIS":
            g["ncRVIS"].dropna().iloc[0]
            if g["ncRVIS"].notna().any()
            else np.nan,

        "Diseases":
            "|".join(
                sorted(diseases)
            ),

        "Number_of_Diseases":
            len(diseases)
    })

gene_summary = pd.DataFrame(
    gene_summary
)

four_disease = gene_summary[
    gene_summary["Number_of_Diseases"] == 4
].copy()

##############################################################################
# BUILD TABLE 5
##############################################################################

table5_rows = []

for _, gene in four_disease.iterrows():

    evidence = gene_evidence[
        gene_evidence["ensg"]
        == gene["Ensembl_ID"]
    ]

    chol = chol_gene_evidence[
        chol_gene_evidence["ensg"]
        == gene["Ensembl_ID"]
    ]

    ml_ids = sorted(
        set(
            evidence[
                "MergedLocusID"
            ]
            .dropna()
            .astype(str)
        )
    )

    
    nonchol_conc = sorted(
        evidence.loc[
            evidence["Pleiotropy_Type"]=="Concordant",
            "Phenotype_Pair"
        ].unique()
    )

    nonchol_disc = sorted(
        evidence.loc[
            evidence["Pleiotropy_Type"]=="Discordant",
            "Phenotype_Pair"
        ].unique()
    )

    chol_conc = sorted(
        chol.loc[
            chol["Pleiotropy_Type"]=="Concordant",
            "Phenotype_Pair"
        ].unique()
    )

    chol_disc = sorted(
        chol.loc[
            chol["Pleiotropy_Type"]=="Discordant",
            "Phenotype_Pair"
        ].unique()
    )

    table5_rows.append({

        "symbol":
            gene["symbol"],

        "Ensembl_ID":
            gene["Ensembl_ID"],

        "Gene_Type":
            gene["Gene_Type"],

        "Chromosome":
            gene["Chromosome"],

        "Gene_Start":
            gene["Gene_Start"],

        "Gene_End":
            gene["Gene_End"],

        "N_Merged_Loci":
            len(ml_ids),

        "MergedLocusIDs":
            "|".join(ml_ids),

        "Diseases":
            gene["Diseases"],

        "pLI":
            gene["pLI"],

        "ncRVIS":
            gene["ncRVIS"],

        "N_Observations":
            len(evidence),

        "N_Phenotype_Pairs":
            evidence["Phenotype_Pair"].nunique(),

        "NonChol_Concordant_Pairs":
            ";".join(nonchol_conc),

        "NonChol_Discordant_Pairs":
            ";".join(nonchol_disc),

        "Cholesterol_Concordant_Pairs":
            ";".join(chol_conc),

        "Cholesterol_Discordant_Pairs":
            ";".join(chol_disc),
    })

table5 = pd.DataFrame(
    table5_rows
)

##############################################################################
# TABLE 5A
# GENE EVIDENCE TABLE
##############################################################################

table5_gene_evidence = gene_evidence.copy()

##############################################################################
# TABLE 5B
# COMPLETE GENE ANNOTATION SUMMARY
##############################################################################

summary_rows = []

for ensg, g in gene_evidence.groupby("ensg"):

    diseases = set()

    for d in g["Disease"]:
        diseases |= {
            x for x in d.split("|")
            if x
        }

    symbol = (
        g["symbol"].dropna().iloc[0]
        if g["symbol"].notna().any()
        else ensg
    )

    merged_loci = sorted(
        set(
            g["MergedLocusID"]
            .dropna()
            .astype(str)
        )
    )

    Mapped_Loci = ";".join(merged_loci)

    phenotype_pairs = sorted(
        g["Phenotype_Pair"]
        .unique()
    )

    ##############################################################################
    # SPLIT NON-CHOLESTEROL PAIRS
    ##############################################################################

    conc_pairs = sorted(
        g.loc[
            g["Pleiotropy_Type"] == "Concordant",
            "Phenotype_Pair"
        ].unique()
    )

    disc_pairs = sorted(
        g.loc[
            g["Pleiotropy_Type"] == "Discordant",
            "Phenotype_Pair"
        ].unique()
    )

    ##############################################################################
    # HDL / LDL SUPPORT
    ##############################################################################

    chol = chol_gene_evidence[
        chol_gene_evidence["ensg"] == ensg
    ]

    chol_conc = sorted(
        chol.loc[
            chol["Pleiotropy_Type"] == "Concordant",
            "Phenotype_Pair"
        ].unique()
    )

    chol_disc = sorted(
        chol.loc[
            chol["Pleiotropy_Type"] == "Discordant",
            "Phenotype_Pair"
        ].unique()
    )

    summary_rows.append({

        "symbol":
            symbol,

        "Ensembl_ID":
            ensg,

        "Gene_Type":
            g["type"].iloc[0],

        "Chromosome":
            g["chr"].iloc[0],

        "Gene_Start":
            g["start"].min(),

        "Gene_End":
            g["end"].max(),

        "pLI":
            g["pLI"].dropna().iloc[0]
            if g["pLI"].notna().any()
            else np.nan,

        "ncRVIS":
            g["ncRVIS"].dropna().iloc[0]
            if g["ncRVIS"].notna().any()
            else np.nan,

        "N_Observations":
            len(g),

        "N_Merged_Loci":
            len(merged_loci),

        "MergedLocusIDs":
            "|".join(merged_loci),

        "Diseases":
            "|".join(sorted(diseases)),

        "N_Diseases":
            len(diseases),

        "N_Phenotype_Pairs":
            len(phenotype_pairs),

        "Concordant_Pairs":
            ";".join(conc_pairs),

        "Discordant_Pairs":
            ";".join(disc_pairs),

        "Cholesterol_Concordant_Pairs":
            ";".join(chol_conc),

        "Cholesterol_Discordant_Pairs":
            ";".join(chol_disc),
    })

table5_gene_summary = (
    pd.DataFrame(summary_rows)
      .sort_values(
          ["N_Diseases",
           "N_Merged_Loci",
           "N_Observations"],
          ascending=False
      )
      .reset_index(drop=True)
)

##############################################################################
# STEP 5
# FIGURE 3 HIGHLIGHT REGIONS
##############################################################################

from itertools import combinations

FOUR_DISEASES = {"CAD", "HT", "STR", "T2D"}

# Current-annotation override for Figure 3 display candidates.
# ENSG00000270316 is currently BORCS7-ASMT
# (formerly C10orf32-ASMT), an ncRNA read-through / NMD candidate.
# Preserve it in the full FUMA evidence, but do not treat it as a
# protein-coding gene for compact Figure 3 labels.
FIGURE3_NONCODING_OVERRIDES = {
    "ENSG00000270316"
}


def split_pipe(value):
    if pd.isna(value):
        return set()

    return {
        x.strip()
        for x in str(value).split("|")
        if x.strip()
    }


def ml_sort_key(value):
    try:
        return int(str(value).replace("ML_", ""))
    except ValueError:
        return str(value)


##############################################################################
# A. GENE-DEFINED FOUR-DISEASE REGIONS
#
# Retain protein-coding genes whose own mapped evidence spans all four
# diseases. Genes belonging to the same connected set of merged loci are
# displayed together as one regional label.
##############################################################################

gene_defined = table5[
    table5["Gene_Type"]
    .astype(str)
    .str.lower()
    .eq("protein_coding")
].copy()

gene_defined = gene_defined[
    ~gene_defined["Ensembl_ID"]
    .astype(str)
    .isin(FIGURE3_NONCODING_OVERRIDES)
].copy()

gene_defined["_DiseaseSet"] = (
    gene_defined["Diseases"]
    .map(split_pipe)
)

gene_defined = gene_defined[
    gene_defined["_DiseaseSet"]
    .map(lambda x: x == FOUR_DISEASES)
].copy()

gene_defined["_LocusSet"] = (
    gene_defined["MergedLocusIDs"]
    .map(split_pipe)
)

##############################################################################
# Group four-disease genes into connected genomic regions according to
# overlap in their merged-locus IDs.
##############################################################################

remaining = set(gene_defined.index)
components = []

while remaining:

    seed = remaining.pop()
    component = {seed}
    component_loci = set(
        gene_defined.at[seed, "_LocusSet"]
    )

    changed = True

    while changed:

        changed = False

        for idx in list(remaining):

            loci = gene_defined.at[
                idx,
                "_LocusSet"
            ]

            if component_loci & loci:

                component.add(idx)
                component_loci |= loci
                remaining.remove(idx)

                changed = True

    components.append(
        (component, component_loci)
    )

highlight_rows = []

for component, locus_ids in components:

    g = gene_defined.loc[
        list(component)
    ].copy()

    # Deterministic display order: genomic position, then symbol.
    g = g.sort_values(
        [
            "Chromosome",
            "Gene_Start",
            "symbol"
        ]
    )

    label_genes = (
        g["symbol"]
        .astype(str)
        .tolist()
    )

    highlight_rows.append({

        "Label":
            "/".join(label_genes),

        "MergedLocusIDs":
            "|".join(
                sorted(
                    locus_ids,
                    key=ml_sort_key
                )
            ),

        "AnnotationType":
            "Gene",

        "LabelBasis":
            "Protein-coding gene(s) with four-disease recurrence"

    })


##############################################################################
# B. LOCUS-DEFINED FOUR-DISEASE REGIONS
#
# Some strict four-disease merged loci contain no individual protein-coding
# gene whose mapped evidence spans all four diseases. For these loci, choose
# the smallest set of protein-coding mapped genes whose combined disease
# evidence covers CAD, HT, STR and T2D.
#
# Among equally small covering sets, prefer the set with the largest summed
# number of unique phenotype-pair observations. Alphabetical ordering is used
# only as the final deterministic tie-breaker.
##############################################################################

annotated_loci = set()

for row in highlight_rows:

    annotated_loci |= split_pipe(
        row["MergedLocusIDs"]
    )


major_loci = table4[
    (table4["Number_of_Diseases"] == 4) &
    (table4["Num_Original_Loci"] >= 5)
].copy()

major_loci = major_loci[
    ~major_loci["MergedLocusID"]
    .isin(annotated_loci)
]


print(
    f"\nAdditional locus-defined annotations: "
    f"{len(major_loci)}"
)


for _, locus in major_loci.iterrows():

    ml = locus["MergedLocusID"]

    genes = gene_evidence[
        gene_evidence["MergedLocusID"]
        .eq(ml)
    ].copy()

    genes = genes[
        genes["type"]
        .astype(str)
        .str.lower()
        .eq("protein_coding")
    ].copy()

    genes = genes[
        ~genes["ensg"]
        .astype(str)
        .isin(FIGURE3_NONCODING_OVERRIDES)
    ].copy()

    genes = genes[
        genes["symbol"].notna()
    ].copy()

    if genes.empty:
        raise ValueError(
            f"{ml}: no protein-coding mapped genes available "
            "for locus-level Figure 3 annotation"
        )

    gene_stats = {}

    for symbol, g in genes.groupby("symbol"):

        disease_set = set()

        for value in g["Disease"]:

            disease_set |= split_pipe(
                value
            )

        gene_stats[symbol] = {

            "Diseases":
                disease_set,

            "N_Phenotype_Pairs":
                g["Phenotype_Pair"]
                .nunique()

        }

    symbols = sorted(
        gene_stats
    )

    best_combo = None
    best_score = None

    for n_genes in range(
        1,
        len(symbols) + 1
    ):

        candidate_combos = []

        for combo in combinations(
            symbols,
            n_genes
        ):

            covered = set()

            for symbol in combo:

                covered |= (
                    gene_stats[symbol]
                    ["Diseases"]
                )

            if FOUR_DISEASES <= covered:

                score = sum(
                    gene_stats[symbol]
                    ["N_Phenotype_Pairs"]

                    for symbol in combo
                )

                candidate_combos.append(
                    (
                        score,
                        combo
                    )
                )

        if candidate_combos:

            # Highest evidence score first;
            # alphabetical tuple as final tie-breaker.
            candidate_combos.sort(
                key=lambda x: (
                    -x[0],
                    x[1]
                )
            )

            best_score, best_combo = (
                candidate_combos[0]
            )

            break

    if best_combo is None:

        raise ValueError(
            f"{ml}: protein-coding mapped genes "
            "do not jointly cover all four diseases"
        )

    highlight_rows.append({

        "Label":
            "/".join(best_combo),

        "MergedLocusIDs":
            ml,

        "AnnotationType":
            "Locus",

        "LabelBasis":
            "Minimal protein-coding gene set covering four-disease locus"

    })


##############################################################################
# FINAL FIGURE 3 ANNOTATION TABLE
##############################################################################

figure3 = pd.DataFrame(
    highlight_rows
)

figure3 = figure3.sort_values(
    [
        "AnnotationType",
        "MergedLocusIDs",
        "Label"
    ]
).reset_index(
    drop=True
)


##############################################################################
# FIGURE 3 QC
##############################################################################

print("\nFigure 3 highlight regions:")
print(
    figure3.to_string(
        index=False
    )
)

assert len(figure3) == 7, (
    f"Expected 7 primary convergence regions, "
    f"found {len(figure3)}"
)

assert (
    figure3["Label"]
    .str.contains(
        r"Y_RNA|snoU13",
        case=False,
        regex=True
    )
    .sum()
    == 0
)

expected_locus_defined = {
    "ML_197",
    "ML_794"
}

observed_locus_defined = set(
    figure3.loc[
        figure3["AnnotationType"]
        .eq("Locus"),
        "MergedLocusIDs"
    ]
)

assert observed_locus_defined == expected_locus_defined, (
    "Unexpected locus-defined Figure 3 regions: "
    f"{sorted(observed_locus_defined)}"
)

##############################################################################
# QC CHOLESTEROL ANNOTATIONS
##############################################################################

print()

print(
    "Loci with cholesterol concordant pairs:",
    (table4["Cholesterol_Concordant"] != "").sum()
)

print(
    "Loci with cholesterol discordant pairs:",
    (table4["Cholesterol_Discordant"] != "").sum()
)

##############################################################################
# EXPORT
##############################################################################

table4.to_csv(
    OUT_DIR /
    "Table4_PleioLoci_noHDL_LDL.csv",
    index=False
)

table5.to_csv(
    OUT_DIR /
    "Table5_GeneLocusMapping_noHDL_LDL.csv",
    index=False
)

figure3.to_csv(
    OUT_DIR /
    "Figure3_HighlightGenes_noHDL_LDL.csv",
    index=False
)

try:
    with pd.ExcelWriter(
        OUT_DIR /
        "Figure3_Framework_noHDL_LDL.xlsx"
    ) as writer:

        table4.to_excel(
            writer,
            sheet_name="Table4_PleioLoci_noHDL_LDL",
            index=False
        )

        table5.to_excel(
            writer,
            sheet_name="Table5_GeneMapping_noHDL_LDL",
            index=False
        )

        figure3.to_excel(
            writer,
            sheet_name="HighlightGenes_noHDL_LDL",
            index=False
        )

        table4_membership.to_excel(
            writer,
            sheet_name="Table4_LocusMembership",
            index=False
        )

        table5_gene_evidence.to_excel(
            writer,
            sheet_name="Table5_GeneEvidence",
            index=False
        )

        table5_gene_summary.to_excel(
            writer,
            sheet_name="Table5_GeneSummary",
            index=False
        )

    print(
        "Excel framework saved:",
        OUT_DIR / "Figure3_Framework_noHDL_LDL.xlsx"
    )

except ImportError:
    print(
        "WARNING: Excel writer dependency not installed; "
        "skipping XLSX export."
    )

table4_membership.to_csv(
    OUT_DIR /
    "Table4_LocusMembership_noHDL_LDL.csv",
    index=False
)

table5_gene_evidence.to_csv(
    OUT_DIR /
    "Table5_GeneEvidence_noHDL_LDL.csv",
    index=False
)

table5_gene_summary.to_csv(
    OUT_DIR /
    "Table5_GeneAnnotationSummary_noHDL_LDL.csv",
    index=False
)

##############################################################################
# FINAL QC
##############################################################################

print("\n=================================================")
print("FINAL QC")
print("=================================================\n")

print("Genes in all four diseases:\n")

for gene in sorted(table5["symbol"]):
    print(gene)

print()
print(f"Table4 rows: {len(table4):,}")
print(f"Table5 rows: {len(table5):,}")
print(f"Highlight regions: {len(figure3):,}")

print("\nDONE")
