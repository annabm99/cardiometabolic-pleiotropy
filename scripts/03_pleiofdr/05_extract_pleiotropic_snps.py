import sys
import pandas as pd
import logging
import os

def GetNames(FileName):
    parts = os.path.basename(FileName).split("_")
    phen1 = "_".join(parts[:2])  # Join the first two elements
    phen2 = "_".join(parts[2:4])  # Join the next two elements

    return phen1, phen2

def GetPleiotropyInfo (pleiotropies, ref, phen1, phen2):
    # Merge pleiotropic SNPs with reference, to get rsIDs and A1, A2
    interest_columns =  [
    'chrnum',
    'chrpos',
    f'zscore_{phen1}',
    f'zscore_{phen2}',
    f'conjfdr_{phen1}_{phen2}',
    f'pval_{phen1}',
    f'pval_{phen2}'
    ]

    pleiotropies['chrnum'] = pleiotropies['chrnum'].astype(str)
    ref['CHR'] = ref['CHR'].astype(str)

    pleiotropies['chrpos'] = pleiotropies['chrpos'].astype(int)
    ref['BP'] = ref['BP'].astype(int)

    pleiotropies = pd.merge(pleiotropies[interest_columns],
                            ref[['CHR', 'BP', 'SNP', 'A1', 'A2']],
                            left_on=['chrnum', 'chrpos'],
                            right_on=['CHR', 'BP'],
                            how='left')
        
    rename_columns = {
        'CHR': 'Chromosome',
        'BP': 'Position',
        f'zscore_{phen1}': f'Z-{phen1}',
        f'zscore_{phen2}': f'Z-{phen2}',
        f'conjFDR_{phen1}_{phen2}': f'CONJFDR-{phen1}_{phen2}',
        f'pval_{phen1}': f'P-{phen1}',
        f'pval_{phen2}': f'P-{phen2}'
    }

    pleiotropies = pleiotropies.rename(columns=rename_columns)
    return pleiotropies

######

PleiotropiesFile = sys.argv[1].rstrip()
RefFile = sys.argv[2].rstrip()
OutDir = sys.argv[3].rstrip()

print(f"Pleiotropy file is {PleiotropiesFile}, reference file is {RefFile}, output directory is {OutDir}")

# Get phenotype names
phen1, phen2 = GetNames(PleiotropiesFile)
print(f"Phenotype 1 is {phen1}, phentotype 2 is {phen2}")

# Log setup
# Log file
log_file_path=(f'{OutDir}/{phen1}_{phen2}-log.txt')
logging.basicConfig(
    filename=log_file_path,
    filemode='w', # overwrites the file each time the script is run
    level=logging.INFO, # captures all INFO-level and higher messages
    format='%(asctime)s - %(levelname)s - %(message)s' # specifies the format of the lof message (time-level-message)
)

sys.stdout = open(log_file_path, 'a')  # Append print statements to the log
sys.stderr = open(log_file_path, 'a')  # Append warnings/errors to the log

# Open Dfs
pleiotropy_df = pd.read_csv(PleiotropiesFile, sep = ",", low_memory=False)
print(f"Pleiotropies dataframe loaded successfully! {pleiotropy_df.shape}")
print(pleiotropy_df.head(5))

ref_df = pd.read_csv(RefFile, sep="\t", low_memory=False)
print(f"Reference dataframe loaded successfully! {ref_df.shape}")
print(ref_df.head(5))

# Get all info from the reference file
pleiotropy_df = GetPleiotropyInfo(pleiotropy_df, ref_df, phen1, phen2)
print(f"Pleiotropies dataframe was merged with the reference to obtain all info: {pleiotropy_df.shape}")
print(pleiotropy_df.head(5))

# Save Df
pleiotropy_df.to_csv(f"{OutDir}/{phen1}_{phen2}-SNPInfo.txt.gz",sep='\t', compression='gzip')