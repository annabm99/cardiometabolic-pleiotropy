import pandas as pd
import os
import logging
import sys

# Setup logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)

# Input file

InFile=sys.argv[1].rstrip()
logging.info(f"Reading input file: {InFile}")
OutDir=sys.argv[2].rstrip()

# Extract prefix before first '-' in filename
filename = os.path.basename(InFile)
prefix = filename.split('-')[0]
logging.info(f"Filename is {filename}, and prefix is {prefix}")

# Output file
output_file = f"{OutDir}/{prefix}-FUMASumstats.txt.gz"
logging.info(f"Output file will be: {output_file}")

# Columns to keep
keep_cols = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "PVAL", "N", "FRQ"]

# Read gzipped TSV
try:
    df = pd.read_csv(InFile, sep='\t', compression='gzip', dtype=str)
    logging.info(f"Input file loaded successfully with {df.shape[0]} rows and {df.shape[1]} columns")
except Exception as e:
    logging.error(f"Failed to read input file: {e}")
    raise

# Select only the columns that exist
existing_cols = [col for col in keep_cols if col in df.columns]
missing_cols = [col for col in keep_cols if col not in df.columns]

logging.info(f"Columns kept: {existing_cols}")
if missing_cols:
    logging.warning(f"Columns missing from input file and skipped: {missing_cols}")

df_out = df[existing_cols].copy()

# Rename PVAL -> P for FUMA
if "PVAL" in df_out.columns:
    df_out.rename(columns={"PVAL": "P"}, inplace=True)
    logging.info("Renamed 'PVAL' to 'P' for FUMA")

# Convert numeric columns in df_out (not df!)
for col in ["CHR", "BP", "N"]:
    if col in df_out.columns:
        df_out[col] = df_out[col].astype(float).astype(int)

for col in ["BETA", "SE", "P", "FRQ"]:
    if col in df_out.columns:
        df_out[col] = df_out[col].astype(float)

# Save gzipped FUMA-ready file
try:
    df_out.to_csv(output_file, sep='\t', index=False, compression='gzip')
    logging.info(f"FUMA-ready gzipped file saved successfully: {output_file}")
except Exception as e:
    logging.error(f"Failed to write output file: {e}")
    raise