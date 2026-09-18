#!/usr/bin/env python3

import os
from pathlib import Path

import pandas as pd


##############################################################################
# PATHS
##############################################################################

PROJECT_DIR = Path(
    os.environ.get(
        "CVP_PROJECT_DIR",
        "/data/samanthafs/scratch/lab_anavarro/anbasquet/"
        "cardiovascular_pleiotropies"
    )
)

INPUT_DIR = PROJECT_DIR / "FinalTables"

OUT_DIR = Path(
    os.environ.get(
        "CVP_SUPPLEMENT_DIR",
        PROJECT_DIR / "FinalTables"
    )
)

OUT_DIR.mkdir(parents=True, exist_ok=True)


##############################################################################
# INPUTS
##############################################################################

TABLE4 = INPUT_DIR / "Table4_PleioLoci_noHDL_LDL.csv"

GENES = (
    INPUT_DIR /
    "Table5_GeneLocusMapping_noHDL_LDL.csv"
)

HIGHLIGHTS = (
    INPUT_DIR /
    "Figure3_HighlightGenes_noHDL_LDL.csv"
)


##############################################################################
# HELPERS
##############################################################################

FOUR_DISEASES = {
    "CAD",
    "HT",
    "STR",
    "T2D"
}

# Current nomenclature/classification overrides.
# Source FUMA/Ensembl annotations are retained separately so that the
# historical analysis remains reproducible.
CURRENT_GENE_ANNOTATION = {
    "ENSG00000270316": {
        "Current_Symbol": "BORCS7-ASMT",
        "Current_Type": "ncRNA read-through (NMD candidate)",
        "Annotation_Note": (
            "Source FUMA symbol C10orf32-ASMT; "
            "current official symbol BORCS7-ASMT"
        ),
    }
}


def split_pipe(value):

    if pd.isna(value):
        return []

    return [
        x.strip()
        for x in str(value).split("|")
        if x.strip()
    ]


def split_semicolon(value):

    if pd.isna(value):
        return []

    return [
        x.strip()
        for x in str(value).split(";")
        if x.strip()
    ]


def join_sorted(values, sep="|"):

    values = {
        str(x).strip()
        for x in values
        if str(x).strip()
    }

    return sep.join(sorted(values))


def ml_number(value):

    return int(
        str(value)
        .replace("ML_", "")
    )


##############################################################################
# LOAD
##############################################################################

table4 = pd.read_csv(TABLE4)
genes = pd.read_csv(GENES)
highlights = pd.read_csv(HIGHLIGHTS)

print("=" * 70)
print("SUPPLEMENTARY TABLE S4")
print("=" * 70)

print("Merged loci:", len(table4))
print("Four-disease recurrent genes:", len(genes))
print("Figure 3 regions:", len(highlights))


##############################################################################
# S4A — RECURRENT FOUR-DISEASE CONVERGENCE REGIONS
##############################################################################

region_rows = []


for _, h in highlights.iterrows():

    label = h["Label"]

    ml_ids = split_pipe(
        h["MergedLocusIDs"]
    )

    ml_set = set(ml_ids)

    loci = table4[
        table4["MergedLocusID"]
        .isin(ml_set)
    ].copy()

    found_ids = set(
        loci["MergedLocusID"]
        .astype(str)
    )

    if found_ids != ml_set:

        raise ValueError(
            f"{label}: missing constituent loci. "
            f"Expected={sorted(ml_set)}, "
            f"found={sorted(found_ids)}"
        )

    chromosomes = set(
        loci["Chromosome"]
        .astype(str)
    )

    if len(chromosomes) != 1:

        raise ValueError(
            f"{label}: region spans multiple chromosomes: "
            f"{sorted(chromosomes)}"
        )

    chromosome = next(
        iter(chromosomes)
    )

    region_start = int(
        loci["Merged_Start"].min()
    )

    region_end = int(
        loci["Merged_End"].max()
    )

    diseases = set()

    for value in loci["Diseases_Involved"]:

        diseases.update(
            split_pipe(value)
        )

    strict_loci = []

    for _, locus in loci.iterrows():

        locus_diseases = set(
            split_pipe(
                locus["Diseases_Involved"]
            )
        )

        if locus_diseases == FOUR_DISEASES:

            strict_loci.append(
                locus["MergedLocusID"]
            )

    strict_loci = sorted(
        strict_loci,
        key=ml_number
    )

    direction_types = set(
        loci["Pleiotropy_Type"]
        .dropna()
        .astype(str)
    )

    if direction_types == {"Concordant"}:

        regional_direction = "Concordant"

    elif direction_types == {"Discordant"}:

        regional_direction = "Discordant"

    elif direction_types == {
        "Concordant",
        "Discordant"
    }:

        regional_direction = "Mixed"

    else:

        raise ValueError(
            f"{label}: unexpected directionality "
            f"{sorted(direction_types)}"
        )

    concordant_pairs = set()
    discordant_pairs = set()

    for _, locus in loci.iterrows():

        pairs = set(
            split_pipe(
                locus["Phenotype_Pairs"]
            )
        )

        if (
            locus["Pleiotropy_Type"]
            == "Concordant"
        ):

            concordant_pairs.update(
                pairs
            )

        elif (
            locus["Pleiotropy_Type"]
            == "Discordant"
        ):

            discordant_pairs.update(
                pairs
            )

    mapped_genes = set()

    for value in loci["Mapped_Genes"]:

        mapped_genes.update(
            split_semicolon(value)
        )

    recurrent_genes = []

    for _, gene in genes.iterrows():

        gene_loci = set(
            split_pipe(
                gene["MergedLocusIDs"]
            )
        )

        if gene_loci & ml_set:

            annotation = CURRENT_GENE_ANNOTATION.get(
                str(gene["Ensembl_ID"]),
                {}
            )

            recurrent_genes.append(
                annotation.get(
                    "Current_Symbol",
                    gene["symbol"]
                )
            )

    constituent_disease_sets = []

    constituent_directions = []

    for _, locus in (
        loci
        .assign(
            _ml=loci[
                "MergedLocusID"
            ].map(ml_number)
        )
        .sort_values("_ml")
        .iterrows()
    ):

        constituent_disease_sets.append(
            f"{locus['MergedLocusID']}:"
            f"{locus['Diseases_Involved']}"
        )

        constituent_directions.append(
            f"{locus['MergedLocusID']}:"
            f"{locus['Pleiotropy_Type']}"
        )

    region_rows.append({

        "Region":
            label,

        "Annotation Type":
            (
                "Gene-defined"
                if h["AnnotationType"] == "Gene"
                else "Locus-defined"
            ),

        "Label Basis":
            h["LabelBasis"],

        "Chr":
            chromosome,

        "Region Start (bp)":
            region_start,

        "Region End (bp)":
            region_end,

        "N Constituent Merged Loci":
            len(ml_ids),

        "Constituent Merged Loci":
            "|".join(
                sorted(
                    ml_ids,
                    key=ml_number
                )
            ),

        "Constituent Locus Directionality":
            "; ".join(
                constituent_directions
            ),

        "Constituent Locus Disease Sets":
            "; ".join(
                constituent_disease_sets
            ),

        "Diseases Across Region":
            "|".join(
                sorted(
                    diseases
                )
            ),

        "N Diseases":
            len(diseases),

        "Contains Strict Four-Disease Merged Locus":
            (
                "Yes"
                if strict_loci
                else "No"
            ),

        "Strict Four-Disease Locus IDs":
            (
                "|".join(
                    strict_loci
                )
                if strict_loci
                else ""
            ),

        "Regional Directionality":
            regional_direction,

        "Concordant Phenotype Pairs":
            join_sorted(
                concordant_pairs
            ),

        "Discordant Phenotype Pairs":
            join_sorted(
                discordant_pairs
            ),

        "Four-Disease Recurrent Genes (All Biotypes)":
            (
                ";".join(
                    sorted(
                        set(
                            recurrent_genes
                        )
                    )
                )
                if recurrent_genes
                else ""
            ),

        "All Mapped Genes":
            (
                ";".join(
                    sorted(
                        mapped_genes
                    )
                )
                if mapped_genes
                else ""
            )
    })


s4a = pd.DataFrame(
    region_rows
)


##############################################################################
# SORT REGIONS
##############################################################################

s4a["_chr_numeric"] = pd.to_numeric(
    s4a["Chr"],
    errors="coerce"
)

s4a = (
    s4a
    .sort_values(
        [
            "_chr_numeric",
            "Region Start (bp)"
        ]
    )
    .drop(
        columns="_chr_numeric"
    )
    .reset_index(
        drop=True
    )
)


##############################################################################
# S4B — TRUE FOUR-DISEASE RECURRENT GENES
##############################################################################

gene_rows = []


for _, gene in genes.iterrows():

    gene_ml = set(
        split_pipe(
            gene["MergedLocusIDs"]
        )
    )

    region_matches = []

    for _, region in s4a.iterrows():

        region_ml = set(
            split_pipe(
                region[
                    "Constituent Merged Loci"
                ]
            )
        )

        if gene_ml & region_ml:

            region_matches.append(
                region["Region"]
            )

    if len(region_matches) != 1:

        raise ValueError(
            f"{gene['symbol']}: expected exactly "
            f"one convergence region, found "
            f"{region_matches}"
        )

    annotation = CURRENT_GENE_ANNOTATION.get(
        str(gene["Ensembl_ID"]),
        {}
    )

    current_symbol = annotation.get(
        "Current_Symbol",
        gene["symbol"]
    )

    current_type = annotation.get(
        "Current_Type",
        gene["Gene_Type"]
    )

    annotation_note = annotation.get(
        "Annotation_Note",
        ""
    )

    gene_rows.append({

        "Convergence Region":
            region_matches[0],

        "Source Gene Symbol":
            gene["symbol"],

        "Current Gene Symbol":
            current_symbol,

        "Ensembl ID":
            gene["Ensembl_ID"],

        "Source Gene Type":
            gene["Gene_Type"],

        "Current Gene Type":
            current_type,

        "Annotation Note":
            annotation_note,

        "Chr":
            gene["Chromosome"],

        "Gene Start (bp)":
            gene["Gene_Start"],

        "Gene End (bp)":
            gene["Gene_End"],

        "N Merged Loci":
            gene["N_Merged_Loci"],

        "Merged Locus IDs":
            gene["MergedLocusIDs"],

        "Diseases":
            gene["Diseases"],

        "N Observations":
            gene["N_Observations"],

        "N Phenotype Pairs":
            gene["N_Phenotype_Pairs"],

        "Non-Cholesterol Concordant Pairs":
            gene["NonChol_Concordant_Pairs"],

        "Non-Cholesterol Discordant Pairs":
            gene["NonChol_Discordant_Pairs"],

        "pLI":
            gene["pLI"],

        "ncRVIS":
            gene["ncRVIS"]
    })


s4b = pd.DataFrame(
    gene_rows
)

s4b["_chr_numeric"] = pd.to_numeric(
    s4b["Chr"],
    errors="coerce"
)

s4b = (
    s4b
    .sort_values(
        [
            "_chr_numeric",
            "Gene Start (bp)",
            "Current Gene Symbol"
        ]
    )
    .drop(
        columns="_chr_numeric"
    )
    .reset_index(
        drop=True
    )
)


##############################################################################
# QC
##############################################################################

print("\n" + "=" * 70)
print("QC")
print("=" * 70)


if len(s4a) != 7:

    raise ValueError(
        f"Expected 7 recurrent regions; "
        f"found {len(s4a)}"
    )


if set(
    s4a["N Diseases"]
) != {4}:

    raise ValueError(
        "Not all convergence regions "
        "span four diseases"
    )


expected_annotation_counts = {
    "Gene-defined": 5,
    "Locus-defined": 2
}

observed_annotation_counts = (
    s4a["Annotation Type"]
    .value_counts()
    .to_dict()
)

if (
    observed_annotation_counts
    != expected_annotation_counts
):

    raise ValueError(
        "Unexpected annotation-type counts: "
        f"{observed_annotation_counts}"
    )


expected_strict = {
    "ML_149",
    "ML_197",
    "ML_794"
}

observed_strict = set()

for value in s4a[
    "Strict Four-Disease Locus IDs"
]:

    observed_strict.update(
        split_pipe(value)
    )


if observed_strict != expected_strict:

    raise ValueError(
        "Unexpected strict four-disease loci: "
        f"{sorted(observed_strict)}"
    )


expected_directionality = {
    "Concordant": 5,
    "Mixed": 2
}

observed_directionality = (
    s4a[
        "Regional Directionality"
    ]
    .value_counts()
    .to_dict()
)

if (
    observed_directionality
    != expected_directionality
):

    raise ValueError(
        "Unexpected regional directionality: "
        f"{observed_directionality}"
    )


if len(s4b) != 12:

    raise ValueError(
        f"Expected 12 four-disease recurrent genes; "
        f"found {len(s4b)}"
    )


if s4b["Ensembl ID"].duplicated().any():

    raise ValueError(
        "Duplicate Ensembl IDs in S4B"
    )


bad_symbols = {
    "Y_RNA",
    "snoU13"
}

if (
    set(
        s4b["Current Gene Symbol"]
    )
    & bad_symbols
):

    raise ValueError(
        "Invalid generic RNA symbols "
        "present among recurrent genes"
    )


##############################################################################
# WRITE
##############################################################################

out_a = (
    OUT_DIR /
    "Supplementary_S4A_RecurrentRegions.csv"
)

out_b = (
    OUT_DIR /
    "Supplementary_S4B_FourDiseaseGenes.csv"
)

s4a.to_csv(
    out_a,
    index=False
)

s4b.to_csv(
    out_b,
    index=False
)


##############################################################################
# REPORT
##############################################################################

print("S4A regions:", len(s4a))
print(
    "Annotation types:",
    observed_annotation_counts
)
print(
    "Regional directionality:",
    observed_directionality
)
print(
    "Strict four-disease loci:",
    sorted(
        observed_strict,
        key=ml_number
    )
)

print("\nS4A:")
print(
    s4a[
        [
            "Region",
            "Annotation Type",
            "Constituent Merged Loci",
            "Contains Strict Four-Disease Merged Locus",
            "Strict Four-Disease Locus IDs",
            "Regional Directionality",
            "Four-Disease Recurrent Genes (All Biotypes)"
        ]
    ].to_string(
        index=False
    )
)

print(
    "\nS4B four-disease recurrent genes:"
)

print(
    s4b[
        [
            "Convergence Region",
            "Source Gene Symbol",
            "Current Gene Symbol",
            "Ensembl ID",
            "Source Gene Type",
            "Current Gene Type"
        ]
    ].to_string(
        index=False
    )
)

print("\nSaved:")
print(out_a)
print(out_b)
print("\nDONE")
