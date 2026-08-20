#!/usr/bin/env python3
"""
GWAS Summary Statistics Formatter
--------------------------------

This script standardizes GWAS summary statistics into a consistent format.
It merges input files with reference variant data, harmonizes allele codes,
calculates missing values if needed, and outputs tab-separated,
gzip-compressed file for downstream analysis.

"""

# LIBRARIES
import sys
import os
import pandas as pd
import logging
import warnings
import numpy as np
from scipy.stats import norm
import gzip

# ---------------------------------------------------------------------------
# PHENOTYPE DICTIONARY
# Dictionary to translate from codes to given names and add the sample size (N)
# Note: N values were extracted from the phenotype descriptions.
# Keys: original phenotype codes
# Values (order from left to right)
#      Custom phenotype code
#      Total sample size
#      Number of cases (only in disease binary phenotypes)
#      Number of controls (only in disease binary phenotypes)

phenDict = {
    "GCST000998" : ["CAD_d", 86995, 22233, 64762],
    "T2D-noUKBB" : ["T2D_d", 441894, 18197, 423697],
    "GCST90086088" : ["HT_d", 28391, 28246],
    "MEGASTROKE" : ["STR_d", 446696, 40585, 406111],
    "21001_raw" : ["BMI_t", 359983],
    "48_raw" : ["WC_t", 360564],
    "4080_raw" : ["SBP_t", 340159],
    "4079_raw" : ["DBP_t", 340162],
    "30870_raw" : ["TGL_t", 343991],
    "30780_raw" : ["LDL_t", 343621],
    "30760_raw" : ["HDL_t",315133],
    "30740_raw" : ["FG_t", 314914]
}

# ---------------------------------------------------------------------------

# FUNCTIONS

def custom_warning_handler(message, category, filename, lineno, file=None, line=None):
    """Redirect warnings into logging."""
    logging.warning(f"{filename}:{lineno}: {category.__name__}: {message}")

def name_and_N(filename, phenDict):
    """Retrieve phenotype name and sample size from filename."""
    for key, values in phenDict.items():
        if key in filename:
            # Create the new filename
            new_name = values[0]
            samplesize = values[1]

            # Check if there are additional values (cases and controls in case of binary)
            if len(values) == 4:
                n_cases = values[2]
                n_controls = values[3]
                print(f"File {filename} references phenotype {new_name} (binary), with sample size {samplesize}, {n_cases} cases and {n_controls} controls.")
                return new_name, samplesize, n_cases, n_controls
            else:
                print(f"File {filename} references phenotype {new_name} (continuous), with sample size {samplesize}.")
                return new_name, samplesize, None, None # No cases and controls
    
    # If no match is found, log an error and exit
    logging.error(f"No match found for {filename}. Exiting.")
    sys.exit(1)


# Function to merge with variants file
def variants_merge(df, variants_path, T2D=False, chr_col=None, bp_col=None):
    """Merge GWAS dataframe with reference variants file."""
    # Check input types
    if not isinstance(df, pd.DataFrame):
        logging.error("Input `df` is not a pandas DataFrame. Exiting the program.")
        sys.exit(1)
    if not isinstance(variants_path, str):
        logging.error("`variants_path` must be a string pointing to the file path. Exiting the program.")
        sys.exit(1)

    # Attempt to load the variants file
    try:
        variants = pd.read_csv(variants_path, compression="gzip", sep="\t", low_memory=False)
        print(f"Variants dataframe loaded successfully: {variants.shape}")
        print(variants.head(n=5))

    except FileNotFoundError:
        print(f"Variants file not found: {variants_path}. Exiting the program.")
        sys.exit(1)
    except pd.errors.ParserError as e:
        logging.exception(f"Error while parsing the file: {e}. Exiting the program.")
        sys.exit(1)
    except Exception as e:
        logging.exception(f"Unexpected error while reading the file: {e}. Exiting the program.")
        sys.exit(1)

    # Validate required columns
    required_columns_variants = {'varid', 'chr', 'pos', 'ref', 'alt', 'rsid'}
    if not required_columns_variants.issubset(variants.columns):
        missing = required_columns_variants - set(variants.columns)
        logging.error(f"Variants file is missing required columns: {missing}. Exiting the program.")
        sys.exit(1)
    if not T2D:
        if 'variant' not in df.columns:
            logging.error("Input DataFrame `df` is missing the 'variant' column. Exiting the program.")
            sys.exit(1)

    # Perform the merge
    try:
        if T2D:
            df = pd.merge(df, variants[['chr', 'pos', 'rsid']], left_on=[chr_col, bp_col], right_on=['chr', 'pos'])
        else:
            df = pd.merge(df, variants[['variant', 'chr', 'pos', 'ref', 'alt', 'rsid']], on='variant')
        print("DataFrame merged with variants file successfully.")
    except KeyError as e:
        logging.exception(f"Error during merge: {e}. Exiting the program.")
        sys.exit(1)
    except Exception as e:
        logging.exception(f"Unexpected error during merge: {e}. Exiting the program.")
        sys.exit(1)
    
    # Count the rows with missing rsIDs
    try:
        na_count = df['rsid'].isna().sum()
        print(f"Number of rows with missing RSIDs when merging with Neale Lab reference file: {na_count}")
    except KeyError:
        logging.error("rsid column not found in the merged DataFrame.")
        raise
    except Exception as e:
        logging.error(f"An error occurred while counting NaN rsid values: {str(e)}")
        raise
    
    # Delete variants dataframe
    del variants

    # Return new df
    return df

# Function to calculate missing values in the T2D dataframe
def CalculateMissing(df, confidence_level=0.95):
    """Calculate missing beta, SE, and Z-score from odds ratios and CIs."""
    required_colummns = ['odds_ratio', 'ci_lower', 'ci_upper', 'p_value']
    for col in required_colummns:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}.")
        
    # Initialize a list for results
    print("Missing beta, standard error, and z-score. Calculating them from odds ratio, confidence intervals and p-values.")
    betas, z_scores, standard_errors, inclusions = [], [], [], []
    z_value = norm.ppf(1 - (1 - confidence_level) / 2)
    
    for _,row in df.iterrows():
        
        odds_ratio, ci_lower, ci_upper, p_value = row['odds_ratio'], row['ci_lower'], row['ci_upper'], row['p_value']
        
        # Validate input values
        if odds_ratio <= 0 or ci_lower <= 0 or ci_upper <= 0:
            raise ValueError(f"Invalid value(s) encountered in row: {row}. Odds ratios and CIs must be positive.")

        # Calculate beta from odds ratio
        beta = np.log(odds_ratio)

        # Method 1: calculate z and se from confidenece intervals
        # Log of the confidence intervals
        log_ci_lower, log_ci_upper = np.log(ci_lower), np.log(ci_upper)

        # Standard error and Z
        se_from_ci = (log_ci_upper - log_ci_lower) / (2 * z_value)
        z_from_ci = beta / se_from_ci

        # Method 2: calculate z from p-value
        z_from_p = abs(norm.ppf(p_value / 2))

        # Cross-validation
        tolerance = 0.1
        inclusion = np.isclose(z_from_p, z_from_ci, atol=tolerance)

        # Append all to the list
        betas.append(beta)
        z_scores.append(z_from_ci)
        standard_errors.append(se_from_ci)
        inclusions.append(inclusion)
    
    # Add z-scores and SE to the dataframe
    df['BETA'], df['Z'], df['SE'], df['INCL'] = betas, z_scores, standard_errors, inclusions
    
    print(f"Filtered size: {df['INCL'].sum()} / {len(df)}")
    return df[df['INCL']].drop(columns=['INCL'])


# Function to get interest columns and rename them
def interest_columns(df, name, n_total, n_cases, n_controls):
    """Select and rename columns of interest to the standard schema.
    Standard output columns: [RSID, CHR, BP, A1, A2, BETA, SE, PVAL, FRQ, (N, N_CASES, N_CONTROLS)]
    """
    # Check if input is a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        logging.error("Input is not a pandas DataFrame.")
        raise TypeError("Input must be a pandas DataFrame.")
    
    # Check if 'name' is a string
    if not isinstance(name, str):
        logging.error("Parameter 'name' must be a string.")
        raise TypeError("Parameter 'name' must be a string.")
    
    try:
        # Check if 'name' ends with '_t' or '_d'
        if name.endswith('_t'):
            required_columns = ['rsid', 'chr', 'pos', 'alt', 'ref', 'beta', 'se', 'pval', 'minor_AF'] # in neale lab, the alternative allele is the effect allele (a1)
        elif name.endswith('_d'):
            if 'STR' in name:
                required_columns = ['MarkerName', 'chr', 'pos', 'Allele1', 'Allele2', 'Effect', 'StdErr', 'P-value', 'Freq1']
            elif 'T2D' in name: # Standard errors are all NA, so we have to calculate them.
                required_columns = ['Consolidated_RSID', 'Chr', 'Pos', 'EA', 'NEA', 'Beta', 'SE', 'Pvalue', 'EAF']
            elif 'HT' in name:
                required_columns = ['Consolidated_RSID', 'chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'standard_error', 'p_value', 'freq_effect_allele']
            elif 'CAD' in name:
                required_columns = ['variant_id', 'chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'standard_error', 'p_value', 'effect_allele_frequency']
        else:
            raise ValueError("The name must end with '_t' or '_d'.")
        
        # Ensure all required columns exist in the DataFrame
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            logging.error(f"Variants file is missing required columns. Found columns: {', '.join(df.columns)}")
            logging.error(f"Missing columns: {', '.join(missing_columns)}")
            raise KeyError(f"Missing columns: {', '.join(missing_columns)}")
        
        # Select only the required columns
        df = df[required_columns]

        # Rename columns
        df.columns = ['RSID', 'CHR', 'BP', 'A1', 'A2', 'BETA', 'SE', 'PVAL', 'FRQ']

        # Add a sample size colummn (N)
        df['N'] = n_total

        # Add case and control columns if needed
        if n_cases is not None and n_controls is not None:
            df['N_CASES'] = n_cases
            df['N_CONTROLS'] = n_controls
        
        print(f"Final columns are: {df.columns}")
            
    except KeyError as e:
        logging.exception(f"KeyError: One or more columns are missing from the DataFrame. {e}")
        raise KeyError(f"Missing columns: {e}")
    except ValueError as e:
        logging.exception(f"ValueError: {e}")
        raise ValueError(e)
    except Exception as e:
        logging.exception(f"Unexpected error: {e}")
        raise Exception(f"Unexpected error: {e}")

    return df

def CheckUpper(df, allele_cols):
    """Ensure allele columns are uppercase."""
    for col in allele_cols:
        if col in df.columns:
            # Check if the first row is already uppercase
            first_row_value = df[col].iloc[0]
            if str(first_row_value).isupper():
                print(f"Column '{col}' is already uppercase. No changes made.")
                continue
            
            # Convert entire column to uppercase if not already
            df[col] = df[col].str.upper()
            print(f"Column '{col}' converted to uppercase.")
    return df

def GetRsIDs(df, reference_vcf, chr_col, bp_col): 
    """Retrieve RSIDs by matching CHR/BP against a reference VCF (.gz)."""
    
    ref_data = []  # List to store information of interest

    try:
        print(f"Opening VCF file: {reference_vcf}")
        with gzip.open(reference_vcf, 'rt') as f:  # Open the .gz file in text mode
            print("Reading file line by line...")
            
            # Iterate over each line in the file
            for line in f:
                if line.startswith('#'):
                    continue  # Skip header lines that start with '#'
                
                # Split the line into columns (assuming VCF format with tab-delimited fields)
                columns = line.strip().split('\t')
                
                # Extract relevant data: CHR, BP, and RSID
                ref_data.append({
                    'CHR': columns[0],  # Chromosome
                    'BP': int(columns[1]),  # Position
                    'RSID': columns[2]  # RsID from VCF
                })
        
        print(f"VCF file {reference_vcf} processed successfully.")
        
    except FileNotFoundError:
        logging.error(f"File not found: {reference_vcf}. Please check the file path.")
        raise
    except gzip.BadGzipFile:
        logging.error(f"Failed to open {reference_vcf}. The file may not be a valid .gz file.")
        raise
    except Exception as e:
        logging.error(f"An unexpected error occurred while reading the VCF file: {str(e)}")
        raise


    # Convert the list of dictionaries to a Pandas DataFrame
    try:
        ref_df = pd.DataFrame(ref_data)
        print(f"Reference DataFrame created with {len(ref_df)} rows.")
    except Exception as e:
        logging.error(f"Failed to convert reference data to DataFrame: {str(e)}")
        raise

    print(f"Reference VCF file (converted to pandas df): {ref_df.shape}.")
    print(ref_df.head(n=5))

    # Make sure the types of columns in df are correct
    df[chr_col] = df[chr_col].astype(str)  # Ensure CHR is a string
    df[bp_col] = df[bp_col].astype(int)    # Ensure BP is an integer
    
    # Merge input DataFrame with reference DataFrame to include rsIDs
    try:
        df = pd.merge(df, ref_df,left_on=[chr_col, bp_col], right_on=['CHR', 'BP'], how='left')
        print(f"DataFrame merged successfully to vcf reference file.")
    except KeyError as e:
        logging.error(f"Merge failed due to missing columns: {e}")
        raise
    except Exception as e:
        logging.error(f"An error occurred during the merge: {str(e)}")
        raise
    
    # Count the rows with missing rsIDs
    try:
        na_count = df['RSID'].isna().sum()
        print(f"Number of rows with missing RSIDs when merging with reference VCF: {na_count}")
    except KeyError:
        logging.error("RSID column not found in the merged DataFrame.")
        raise
    except Exception as e:
        logging.error(f"An error occurred while counting NaN RSID values: {str(e)}")
        raise
    
    # Delete vcf dataframe
    del ref_df

    # Return the enriched DataFrame
    return df

def GetChrBp(df, rsid_col, a1_col, a2_col, variants_path): 
    """Retrieve CHR and BP using a variants reference table keyed by (rsid, ref, alt)."""

    # Merge with Neale
    if not isinstance(df, pd.DataFrame):
        logging.error("Input `df` is not a pandas DataFrame. Exiting the program.")
        sys.exit(1)
    if not isinstance(variants_path, str):
        logging.error("`variants_path` must be a string pointing to the file path. Exiting the program.")
        sys.exit(1)

    # Attempt to load the variants file
    try:
        variants = pd.read_csv(variants_path, compression="gzip", sep="\t", low_memory=False)
        print(f"Variants dataframe loaded successfully: {variants.shape}")
        print(variants.head(n=5))

    except FileNotFoundError:
        print(f"Variants file not found: {variants_path}. Exiting the program.")
        sys.exit(1)
    except pd.errors.ParserError as e:
        logging.exception(f"Error while parsing the file: {e}. Exiting the program.")
        sys.exit(1)
    except Exception as e:
        logging.exception(f"Unexpected error while reading the file: {e}. Exiting the program.")
        sys.exit(1)

    # Perform the merge
    try:
        print(f"Variants cols: {variants.columns}")
        df = pd.merge(df, variants[['chr', 'pos', 'rsid', 'ref', 'alt']], left_on=[rsid_col, a1_col, a2_col], right_on=['rsid', 'ref', 'alt'])
        print("DataFrame merged with variants file successfully.")
    except KeyError as e:
        logging.exception(f"Error during merge: {e}. Exiting the program.")
        sys.exit(1)
    except Exception as e:
        logging.exception(f"Unexpected error during merge: {e}. Exiting the program.")
        sys.exit(1)
    
    # Count the rows with missing bp
    try:
        na_count = df['pos'].isna().sum()
        print(f"Number of rows with missing BP when merging with Neale Lab reference file: {na_count}")
    except KeyError:
        logging.error("rsid column not found in the merged DataFrame.")
        raise
    except Exception as e:
        logging.error(f"An error occurred while counting NaN rsid values: {str(e)}")
        raise
    
    # Delete variants dataframe
    del variants

    # Return new df
    return df

def ConsolidateRsids(df):
    """Consolidate RSID/rsid into a single column Consolidated_RSID."""
    try:
        if 'RSID' not in df.columns or 'rsid' not in df.columns:
            missing_cols = [col for col in ['RSID', 'rsid'] if col not in df.columns]
            raise KeyError(f"Missing required columns: {', '.join(missing_cols)}")
        
        mismatched_count = 0  # Counter for mismatched RSIDs

        def resolve_rsids(row):
            nonlocal mismatched_count
            rsid1, rsid2 = row['RSID'], row['rsid']

            if pd.isna(rsid1) and pd.isna(rsid2):
                return np.nan
            elif pd.isna(rsid1):
                return rsid2
            elif pd.isna(rsid2):
                return rsid1
            elif rsid1 == rsid2:
                return rsid1
            else:
                mismatched_count += 1  # Increment counter instead of storing the row
                return np.nan

        df['Consolidated_RSID'] = df.apply(resolve_rsids, axis=1)
        df = df.dropna(subset=['Consolidated_RSID'])

        print(f"Total mismatched RSIDs removed: {mismatched_count}")
        print(f"Number of remaining rows: {len(df)}")


        df = df.drop(columns=['RSID', 'rsid'], errors='ignore')
        return df

    except KeyError as e:
        print(f"Error: {e}")
        logging.error(f"KeyError: {e}")
        raise

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        logging.error(f"Unexpected error: {e}")
        raise

# ---------------------------------------------------------------------------

# MAIN EXECUTION

# Import arguments
FilePath = sys.argv[1].rstrip()
print(f"File path is {FilePath}.")
OutputDir = sys.argv[2].rstrip()
print(f"Output directory is {OutputDir}.")
VariantsFile=sys.argv[3].rstrip()
print(f"Variants file is {VariantsFile}.")
RefVcfFile = sys.argv[4].rstrip()
print(f"Reference VCF file is {RefVcfFile}")

print("ALL ARGUMENTS IMPORTED")

# Get phenotype info
filename = os.path.basename(FilePath)
print(f"Filename is {filename}")
new_filename, n_total, n_cases, n_controls= name_and_N(filename, phenDict)

# Log file
log_file_path=(f'{OutputDir}/{new_filename}-log.txt')
logging.basicConfig(
    filename=log_file_path,
    filemode='w', # overwrites the file each time the script is run
    level=logging.INFO, # captures all INFO-level and higher messages
    format='%(asctime)s - %(levelname)s - %(message)s' # specifies the format of the lof message (time-level-message)
)

# Redirect stdout and stderr to the log file
sys.stdout = open(log_file_path, 'a')  # Append print statements to the log
sys.stderr = open(log_file_path, 'a')  # Append warnings/errors to the log

warnings.showwarning = custom_warning_handler # Custom warning handler, instead of the python default

### STEP 2: OPEN FILE
try:
    df = pd.read_csv(FilePath, sep = " ", low_memory=False)
    print("DataFrame loaded successfully")
    print(df.head(n=5))
except Exception as e:
    logging.exception(f"Error while reading the file: {e}")
    sys.exit(1) # Exit with a non-zero status (indicates failure)

print(f"DataFrame shape after loading: {df.shape}")

### STEP 3: MERGE WITH REFERECE FILE (IF NEEDED)

# Run this if the file comes form UKBB Neale Lab (all trait files)
if new_filename.endswith('_t'):
    df=variants_merge(df, VariantsFile) # Run the merge function

if new_filename.startswith('T2D'):  
    # Retrieve Rsids from the SNPdb VCF file
    df=GetRsIDs(df, RefVcfFile, chr_col='Chr', bp_col='Pos')
    print(f"Dataframe after retrieving RSIDs from reference VCF: {df.shape}")
    print(df.head(n=3))

    # Retrieve Rsids from the Neale Lab variants file (to complete the ones not in the Vcf)    
    df = variants_merge(df, VariantsFile, T2D=True, chr_col='Chr', bp_col='Pos')
    print(f"Dataframe after retrieving rsids from reference variants file: {df.shape}")
    print(df.head(n=3))

    df=ConsolidateRsids(df)

if new_filename.startswith('STR'):
    # Retrieve Chr and Pos from the SNPdb VCF file
    df = CheckUpper(df, ["Allele1", "Allele2"])
    df = GetChrBp(df, rsid_col='MarkerName', a1_col='Allele1', a2_col='Allele2', variants_path=VariantsFile)

if new_filename.startswith('HT'):
    # Retrieve rsids
    df=GetRsIDs(df, RefVcfFile, chr_col='chromosome', bp_col='base_pair_location')
    print(f"Dataframe after retrieving rsids from reference variants file: {df.shape}")
    print(df.head(n=3))

    df = variants_merge(df, VariantsFile, T2D=True, chr_col='chromosome', bp_col='base_pair_location')
    print(f"Dataframe after retrieving rsids from reference variants file: {df.shape}")
    print(df.head(n=3))

    df=ConsolidateRsids(df)

print(f"Final dataframe shape after merging with reference: {df.shape}")
print(df.head(n=5))

### STEP 4: GET NEEDED COLUMNS

df = interest_columns(df, new_filename, n_total, n_cases, n_controls)

print(f"DataFrame shape after getting interest columns: {df.shape}")

# Convert allele columns to upper-case if necessary
df  = CheckUpper(df, ['A1', 'A2'])

print("Final dataframe:")
print(df.head(n=5))

### STEP 5: SAVE NEW FILE

try:
    # Ensure the output directory exists
    os.makedirs(OutputDir, exist_ok=True)
    new_file_path = os.path.join(OutputDir, f"{new_filename}-formatted.csv.gz")

    # Save the DataFrame as a tab-separated, gzip-compressed CSV file
    df.to_csv(new_file_path, sep='\t', index=False, na_rep='NA')

    print(f"File successfully saved as: {OutputDir}/{new_filename}-formatted.csv.gz")

except Exception as e:
    # Log the error if something goes wrong
    logging.exception(f"An error occurred while saving the file: {e}")
