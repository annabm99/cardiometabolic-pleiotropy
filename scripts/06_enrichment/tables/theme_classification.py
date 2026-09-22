#!/usr/bin/env python3

"""
Shared GWAS Catalog theme classification.

Used by the recurrent GWAS Catalog trait-enrichment analysis underlying
Figure 4B and its pair-level permutation sensitivity analysis.

Important
---------
The themes summarize GWAS Catalog phenotype/trait annotations. They are
not molecular pathway categories.

Eight source terms are retained in the recurrent-term catalogue but are
excluded from thematic aggregation because they are composite,
insufficiently specific, or do not have a clear unique assignment to one
of the predefined themes.
"""


# =============================================================================
# EXPLICIT THEMATIC EXCLUSIONS
# =============================================================================

AMBIGUOUS_THEME_EXCLUSIONS = {

    "Body mass index and systole blood pressure (pairwise)":
        "Composite cross-domain phenotype: adiposity and blood pressure.",

    "Alzheimer's disease or HDL levels (pleiotropy)":
        "Composite cross-domain phenotype: neurological and lipid traits.",

    "Intracranial":
        "Insufficiently specific source phenotype label.",

    "Haemorrhoidal disease":
        "No clear unique assignment to the predefined themes.",

    "Hematological and biochemical traits":
        "Broad composite phenotype spanning more than one biological domain.",

    "Endometriosis":
        "No clear unique assignment to the predefined themes.",

    "Arsenic metabolism":
        "No clear unique assignment to the predefined themes.",

    "Spleen volume":
        "No clear unique assignment to the predefined themes.",
}


# =============================================================================
# THEME ASSIGNMENT
# =============================================================================

def assign_theme(pathway):

    p = str(pathway).lower()

    # -------------------------------------------------------------------------
    # Cardiovascular Disease
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "coronary artery disease",
        "coronary heart disease",
        "myocardial infarction",
        "ischemic stroke",
        "cardioembolic stroke",
        "lacunar stroke",
        "intracranial aneurysm",
        "intracerebral hemorrhage",
        "stroke",
        "aortic root",
        "electrocardiogram",
        "intracranial",
        "haemorrhoidal disease"

    ]):

        return "Cardiovascular Disease"

    # -------------------------------------------------------------------------
    # Blood Pressure & Hemodynamics
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "blood pressure",
        "systolic blood pressure",
        "diastolic blood pressure",
        "mean arterial pressure",
        "pulse pressure",
        "hypertension",
        "resting heart rate"

    ]):

        return "Blood Pressure & Hemodynamics"

    # -------------------------------------------------------------------------
    # Adiposity & Body Composition
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "body mass index",
        "waist",
        "hip circumference",
        "hip index",
        "weight",
        "body fat distribution",
        "appendicular lean mass",
        "body shape",
        "adult body size",
        "mineral density"

    ]):

        return "Adiposity & Body Composition"

    # -------------------------------------------------------------------------
    # Lipids & Lipoproteins
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "hdl",
        "ldl",
        "cholesterol",
        "triglyceride",
        "apolipoprotein"

    ]):

        return "Lipids & Lipoproteins"

    # -------------------------------------------------------------------------
    # Glycaemic & Metabolic
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "type 2 diabetes",
        "fasting insulin",
        "glycated hemoglobin",
        "metabolic syndrome",
        "leptin",
        "serum metabolite"

    ]):

        return "Glycaemic & Metabolic"

    # -------------------------------------------------------------------------
    # Hematological
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "platelet",
        "reticulocyte",
        "hematocrit",
        "hemoglobin",
        "white blood cell",
        "myeloid",
        "neutrophil",
        "monocyte",
        "eosinophil",
        "red blood cell",
        "corpuscular",
        "hematological",
        "red cell"

    ]):

        return "Hematological"

    # -------------------------------------------------------------------------
    # Neurological / Psychiatric / Cognitive
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "brain morphology",
        "cortical",
        "subcortical",
        "schizophrenia",
        "autism",
        "anorexia",
        "mood instability",
        "intelligence",
        "cognitive",
        "alzheimer",
        "white matter",
        "hearing difficulty",
        "cbt",
        "sleep duration",
        "multisite chronic pain"

    ]):

        return "Neurological & Psychiatric"

    # -------------------------------------------------------------------------
    # Immune / Inflammatory
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "crohn",
        "ulcerative colitis",
        "inflammatory bowel disease",
        "asthma",
        "vitiligo",
        "c-reactive protein",
        "endometriosis"

    ]):

        return "Immune & Inflammatory"

    # -------------------------------------------------------------------------
    # Renal
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "glomerular filtration",
        "chronic kidney disease",
        "blood urea nitrogen",
        "albumin-to-creatinine",
        "urate",
        "uric acid"

    ]):

        return "Renal Function"

    # -------------------------------------------------------------------------
    # Cancer
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "cancer",
        "glioma",
        "glioblastoma",
        "carcinoma"

    ]):

        return "Cancer"

    # -------------------------------------------------------------------------
    # Ophthalmological
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "glaucoma",
        "intraocular",
        "refractive error"

    ]):

        return "Ophthalmological"

    # -------------------------------------------------------------------------
    # Longevity & Aging
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "parental longevity",
        "parental lifespan",
        "attained age",
        "age at death"

    ]):

        return "Longevity & Aging"

    # -------------------------------------------------------------------------
    # Lifestyle & Behaviour
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "smoking",
        "coffee",
        "fruit consumption",
        "morning person",
        "gym",
        "religious group",
        "sexual intercourse",
        "household income"

    ]):

        return "Lifestyle & Behaviour"

    # -------------------------------------------------------------------------
    # Medication Use
    # -------------------------------------------------------------------------

    if "medication use" in p:

        return "Medication Use"

    # -------------------------------------------------------------------------
    # Liver / Biomarkers
    # -------------------------------------------------------------------------

    if any(x in p for x in [

        "aminotransferase",
        "alkaline phosphatase",
        "non-albumin protein",
        "spleen volume",
        "arsenic metabolism",
        "alanine transaminase",
        "liver enzyme"

    ]):

        return "Liver & Biomarkers"

    return "UNCLASSIFIED"
