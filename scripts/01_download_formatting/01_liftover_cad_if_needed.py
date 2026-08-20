#!/usr/bin/python

"""
# Liftover Script: Convert genomic coordinates to hg19
# Applies to CAD_d (hg18 → hg19)
"""

# Import required libraries
import numpy as np
import pandas as pd
from pyliftover import LiftOver
import sys
import logging
import os

####### FUNCTIONS #######

def GetInfo(InputFile):
    """Extract filename and phenotype name from input path."""
    filename = os.path.basename(InputFile)
    name = filename.split("-")[0]
    return filename, name


def SetUpLogging(out_dir, name):
    """Setup logging to file in the output directory."""
    log_file_path = f'{out_dir}/{name}-log.txt'
    logging.basicConfig(
        filename=log_file_path,
        filemode='w',  # Overwrites the file each time the script is run
        level=logging.INFO,  # Captures all info-level and higher messages
        format='%(asctime)s - %(levelname)s - %(message)s'  # Log message format
    )
    sys.stdout = open(log_file_path, 'a')  # Append print statements to the log
    sys.stderr = open(log_file_path, 'a')  # Append warnings/errors to the log


def ReadInputFile(input_file):
    """Read gzip-compressed, tab-separated input file into a DataFrame."""

    try:
        return pd.read_csv(input_file, sep="\t", compression="gzip", low_memory=False)
    except Exception as e:
        logging.error(f"Failed to read input file: {e}")
        sys.exit(1)


def ConvertCoordinates(df, chr_column, bp_column, og_build):
    """Convert genomic coordinates to hg19 using pyliftover."""

    # Initialize liftover object
    lo = LiftOver(og_build, 'hg19')
    print(f"Converting from genome build {og_build} to hg19.")

    # Empty list to store converted coordinates
    converted = []

    logging.info(f"There are {len(df)} rows to convert.")

    # Ensure BP column has valid numeric values
    df[bp_column] = pd.to_numeric(df[bp_column], errors='coerce')  # Convert invalid strings to NaN

    for index, row in df.iterrows():
        try:
            bp = int(row[bp_column]) - 1  # Convert to 0-based indexing, has to be integer
            chr = str(row[chr_column])

            converted_coord = lo.convert_coordinate(f'chr{chr}', bp)
            if not converted_coord:
                converted.append((None, None, None, None))
            else:
                converted.append(converted_coord[0])

        except Exception as e:
            logging.warning(f"Failed to convert row {index}: {e}")
            converted.append((None, None, None, None))

    lift = pd.DataFrame(converted, columns=['CHR', 'BP', 'STRAND', 'CONFIDENCE']) # New df with the converted coordinates.
    df.rename(columns={chr_column: 'CHR_ORI', bp_column: 'BP_ORI'}, inplace=True) # Rename the original dataframe.
    df = pd.concat([df, lift[['CHR', 'BP']]], axis=1)  # Concatenate along columns.
    df['CHR'] = df['CHR'].astype(str).str.replace("chr", "", regex=False) # Remove 'chr' before the chromosome number.
    
    # Convert BP back to 1-based indexing
    if 'BP' in df.columns:
        df['BP'] = df['BP'].apply(lambda x: int(x) + 1 if pd.notna(x) else x)
    return df


def FilterNaN(df, columns_to_check):
    """Remove rows with missing or invalid values in specified columns."""

    for col in columns_to_check:
        df[col] = df[col].astype(str)  # Ensure columns are strings for filtering
        df = df[~df[col].str.contains("NaN|NA|Nan", case=False, na=False) & (df[col] != '')]
        df = df[~df[col].isna()]
    return df


def IsChromosome(val):
    """Check if a value is a valid chromosome (1-22, X, Y)."""
  
    try:
        # Return True if the value can be converted to an integer
        int(val)
        return True
    except ValueError:
        # Return True if the value is "X" or "Y" (case-insensitive)
        return str(val).upper() in ["X", "Y"]
    
def IntegerBp(df, bp_column):
    """Ensure base-pair column contains only integers."""
    try:
        # Convert to numeric and coerce errors (non-convertible values become NaN)
        df[bp_column] = pd.to_numeric(df[bp_column], errors='coerce')

        # Drop rows where bp_column is NaN after conversion
        df = df.dropna(subset=[bp_column])

        # Convert column to integer type (ensures all values are integers)
        df[bp_column] = df[bp_column].astype(int)

        print("All bp values are integers.")

        return df
    
    except Exception as e:
        logging.error(f"Error ensuring {bp_column} contains only integers: {e}")
        sys.exit(1)


def SaveDf(df, output_dir, name, termination):
    """Save dataframe as gzip-compressed TSV."""

    outfile = str(output_dir + "/" + name + termination)
    df.to_csv(
        outfile,
        sep='\t',
        index_label=None,
        index=False,
        quoting=None,
        decimal='.',
        compression='gzip'
    )

    return outfile

# -----------------------------
# MAIN SCRIPT
# -----------------------------

# Genome build for each phenotype
BuildDict = {
    "CAD_d" : "hg18",
}

# Get command-line arguments
InputFile = sys.argv[1].rstrip()
OutDir = sys.argv[2].rstrip()

ChrCol = 'CHR' # Name of the chromosome column
BpCol = 'BP' # Name of the base-pair position column

Filename, Name = GetInfo(InputFile)
print (f"Working on phenotype {Name}, in file {Filename}.")

SetUpLogging(OutDir, Name)

df = ReadInputFile(InputFile)
print(f"Input dataframe: {df.shape}")
print(f"Original columns: {df.columns}")

df = ConvertCoordinates(df, ChrCol, BpCol, og_build=BuildDict[Name])
print(f"Dataframe after conversion: {df.shape}, with columns: {df.columns}")
print(df.head(10))

df = FilterNaN(df, [ChrCol, BpCol])
print(f'Dataframe after filtering NaN rows: {df.shape}')


# Remove invalid chromosome rows
df = df[df['CHR'].apply(IsChromosome)]

# Ensure bp column are integers
df = IntegerBp(df, BpCol)
print(f"Dataframe after excluding non-integer bp positions {df.shape}")

print(f'Final dataframe: {df.shape}, with columns: {df.columns}')

# Save the final result as a gzip-compressed file
OutFile = SaveDf(df, OutDir, Name, "-hg19.gz")

print(f"Dataframe saved to {OutFile}")
