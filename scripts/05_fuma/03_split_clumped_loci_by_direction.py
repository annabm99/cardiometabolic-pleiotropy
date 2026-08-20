import os
import sys
import pandas as pd
import logging

#### FUNCTIONS

def setup_logging(outdir, inpair):
    """Configure logging to file and stdout"""
    os.makedirs(outdir, exist_ok=True)
    log_file = os.path.join(outdir, f"{inpair}_log.txt")

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # Clear existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    # File handler
    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.DEBUG)

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    # Formatter
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger


def process_pair(df_pos, df_neg, df_snps, df_loci):
    """
    Filter pleiotropy positive/negative files with SNPs from clumped snps/loci.
    Returns: pos_lead, neg_lead, pos_snps, neg_snps (all DataFrames).
    Logs missing SNPs.
    """

    try:

        pos_snps = set(df_pos["SNP"].astype(str))
        neg_snps  = set(df_neg["SNP"].astype(str))

        pos_lead = df_loci[df_loci["LEAD_SNP"].isin(pos_snps)]
        neg_lead = df_loci[df_loci["LEAD_SNP"].isin(neg_snps)]

        pos_snps = df_snps[df_snps["LEAD_SNP"].isin(pos_snps)]
        neg_snps = df_snps[df_snps["LEAD_SNP"].isin(neg_snps)]

 
        # Log summary stats
        total_leads = len(df_loci)
        total_snps = len(df_snps)

        retained_leads = len(pos_lead) + len(neg_lead)
        retained_snps = len(pos_snps) + len(neg_snps)

        logger.info(f"Lead SNPs retained: {retained_leads}, total initial: {total_leads}, ratio: {retained_leads}/{total_leads} "
                    f"(pos={len(pos_lead)}, neg={len(neg_lead)})")

        logger.info(f"Candidate SNPs retained: {retained_snps}, total initial:{total_snps}, ratio: {retained_snps}/{total_snps} "
                    f"(pos={len(pos_snps)}, neg={len(neg_snps)})")

        return pos_lead, neg_lead, pos_snps, neg_snps

    except Exception as e:
        logger.exception(f"Error in process_pair: {e}")
        raise   # re-raise so main also knows it failed


#### MAIN

if len(sys.argv) != 5:
    print("Usage: python script.py <InPair> <SepDir> <ClumpDir> <OutDir>")
    sys.exit(1)

InPair = sys.argv[1].rstrip()
SepDir = sys.argv[2].rstrip()
ClumpDir = sys.argv[3].rstrip()
OutDir   = os.path.join(sys.argv[4].rstrip(), InPair)

logger = setup_logging(OutDir, InPair)
logger.info(f"Starting analysis for {InPair}")
logger.info(f"SepDir: {SepDir}, ClumpDir: {ClumpDir}, OutDir: {OutDir}")


# Load Files
InPair_NoHyphen = InPair.replace("-", "_")
pos_file = os.path.join(SepDir, f"{InPair_NoHyphen}_positive.csv.gz")
neg_file = os.path.join(SepDir, f"{InPair_NoHyphen}_negative.csv.gz")

snp_file  = os.path.join(ClumpDir, InPair, f"{InPair}-FdrClumped.snps.csv")
loci_file = os.path.join(ClumpDir, InPair, f"{InPair}-FdrClumped.loci.csv")

# Read files
try:
    df_pos = pd.read_csv(pos_file, sep="\t", compression="gzip")
    df_neg = pd.read_csv(neg_file, sep="\t", compression="gzip")
    df_snps = pd.read_csv(snp_file, sep="\t")
    df_loci = pd.read_csv(loci_file, sep="\t")
    logger.info("All input files loaded successfully")

except Exception as e:
    logger.error(f"File reading failed:{e}")
    sys.exit(1)

# Process dfs
pos_lead, neg_lead, pos_snps, neg_snps = process_pair(df_pos, df_neg, df_snps, df_loci)

# Save new dataframes
os.makedirs(OutDir, exist_ok=True)


outfiles = {
    "loci_positive": os.path.join(OutDir, f"{InPair}-FDRClumped.loci_positive.csv"),
    "loci_negative": os.path.join(OutDir, f"{InPair}-FDRClumped.loci_negative.csv"),
    "snps_positive": os.path.join(OutDir, f"{InPair}-FDRClumped.snps_positive.csv"),
    "snps_negative": os.path.join(OutDir, f"{InPair}-FDRClumped.snps_negative.csv"),
}

pos_lead.to_csv(outfiles["loci_positive"], index=False, sep="\t")
neg_lead.to_csv(outfiles["loci_negative"], index=False, sep="\t")
pos_snps.to_csv(outfiles["snps_positive"], index=False, sep="\t")
neg_snps.to_csv(outfiles["snps_negative"], index=False, sep="\t")

for label, path in outfiles.items():
    logger.info(f"Saved {label} file to: {path}")

logger.info("Finished successfully.")
