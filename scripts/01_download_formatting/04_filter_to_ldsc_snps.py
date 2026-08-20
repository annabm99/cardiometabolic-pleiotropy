import pandas as pd # type: ignore
import sys
import os
import logging
import matplotlib.pyplot as plt # type: ignore
import numpy as np # type: ignore

## DICTIONARIES

NamesDict = {
    "CAD_d" : "Coronary Artery Disease",
    "T2D_d" : "Type 2 Diabetes",
    "HT_d" : "Hypertension",
    "STR_d" : "Stroke",
    "BMI_t" : "Body Mass Index (BMI)",
    "WC_t" : "Waist Circumference",
    "SBP_t" : "Systolic Blood Pressure",
    "DBP_t" : "Diastolic Blood Pressure",
    "TGL_t" : "Triglycerides",
    "LDL_t" : "Low-density Lipoprotein (LDL) Cholesterol",
    "HDL_t" : "High-density Lipoprotein (HLD) Cholesterol",
    "FG_t" : "Fasting Glucose"
}


def GetInfo(InputFile):
    filename = os.path.basename(InputFile)
    name = filename.split("-")[0]
    
    return filename, name

def SetUpLogging(out_dir, name):
    """
    Sets up logging configuration.

    Parameters:
        out_dir (str): Output directory for the log file.
        name (str): name used to create the log file name (usually same name as the output file, only the first part).
    """
    log_file_path = f'{out_dir}/{name}-log.txt'
    logging.basicConfig(
        filename=log_file_path,
        filemode='w',  # Overwrites the file each time the script is run
        level=logging.INFO,  # Captures all info-level and higher messages
        format='%(asctime)s - %(levelname)s - %(message)s'  # Log message format
    )
    sys.stdout = open(log_file_path, 'a')  # Append print statements to the log
    sys.stderr = open(log_file_path, 'a')  # Append warnings/errors to the log

# Function to get SNPs in the original GWAS matching the ones after munge filtering

def ReadInputFile(input_file):
    """
    Reads the compressed (gzip) and tab-separated input file as a pandas DataFrame.

    Parameters:
        input_file (str): Path to the input file.

    Returns:
        pd.DataFrame: The input DataFrame.
    """
    try:
        return pd.read_csv(input_file, sep="\t", compression="gzip", low_memory=False)
    except Exception as e:
        logging.error(f"Failed to read input file: {e}")
        sys.exit(1)

def GetFiltered(df_og, df_munge):

    # Drop NA rows from the munged df
    df_munge = df_munge.dropna()
    print(df_munge.head(5))
    # Merge the two DataFrames on the 'SNP' column
    filtered = df_og.merge(df_munge[['SNP', 'Z']], left_on='RSID', right_on='SNP', how='inner')
    filtered = filtered.drop('RSID', axis=1)

    return filtered

def PvalDist(df, name, names_dict):

    p_values=df['PVAL'].dropna()

    # Ensure there is data to plot
    if p_values.empty:
        print("Warning: No valid p-values to plot!")
        return None  # Return nothing if data is empty
    
    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot the histogram
    ax.hist(p_values, bins=50, color = one_color, edgecolor = 'black')

    # Add title and labels
    ax.set_xlabel('P-value')
    ax.set_ylabel('Count')
    ax.set_title(f'Distribution of P-values in {names_dict[name]}', fontsize=20)

    # Return the figure object
    return fig

def Manhattan(df, name, names_dict, threshold):

    colors = two_colors  #  Alternating colors for chromosomes

    # Calculate -log(pval)
    df = df.copy()  # Avoid modifying the original DataFrame
    df['log_pval'] = -np.log10(df['PVAL'])

    # Get unique chromosomes and sort them properly (1, 2, ..., 22, X)
    unique_chromosomes = sorted(df['CHR'].unique(), key=lambda x: (int(x) if str(x).isdigit() else float('inf')))
    chromosome_color_map = {chrom: colors[i % 2] for i, chrom in enumerate(unique_chromosomes)}

    # Map chromosomes to a continuous x-axis, computing genome positions
    x_labels, x_ticks= [], []
    current_pos = 0
    gap = 5e6  # 2 million base pairs gap between chromosomes, adjust as needed

    genome_positions = []
    for chrom in unique_chromosomes:
        chrom_df = df[df['CHR'] == chrom]
        x_labels.append(f"Chr {chrom}")
        x_ticks.append(current_pos + (chrom_df['BP'].max() - chrom_df['BP'].min()) / 2)
        genome_positions.extend(chrom_df['BP'] + current_pos)
        current_pos += chrom_df['BP'].max() - chrom_df['BP'].min() + gap  # Add gap between chromosomes
        
    df['genome_pos'] = genome_positions

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot each chromosome separately to alternate colors
    for chrom in unique_chromosomes:
        chrom_df = df[df['CHR'] == chrom]
        ax.scatter(
            chrom_df['genome_pos'], chrom_df['log_pval'], 
            color=chromosome_color_map[chrom], s=10, label=None  # Remove chromosome labels
        )

    # Add threshold line
    ax.axhline(y=threshold, color='red', linestyle='--', label='Significance threshold (p<5e-08)')

    # Set labels and title
    ax.set_title(f'Manhattan Plot for {names_dict[name]}', fontsize=23)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(p-value)')

    # Adjust x-axis labels
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation=90)

    # Remove and side borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Adjust layout
    plt.legend()
    plt.tight_layout()

    return fig  # Return the figure object

def EffectSizeDist(df, name, names_dict):
    # Set up the figure with 1 row and 2 columns.
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Plot the distribution of effect sizes (β) on the left.
    axes[0].hist(df['BETA'], bins=50, color=one_color, alpha=0.7, edgecolor='black', density=False)
    axes[0].set_title(f"Distribution of Effect Sizes (β)\nin {names_dict[name]}", fontsize=20)
    axes[0].set_xlabel("Effect Size (β)", fontsize=14)
    axes[0].set_ylabel("Density", fontsize=14)

    # Compute mean and standard deviation of BETA.
    beta_mean = df['BETA'].mean()
    beta_std = df['BETA'].std()

    # Generate theoretical normal curve for BETA.
    x_vals_beta = np.linspace(beta_mean - 4 * beta_std, beta_mean + 4 * beta_std, 1000)
    y_vals_beta = (1 / (beta_std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x_vals_beta - beta_mean) / beta_std) ** 2)
    # Scale the density to match the histogram counts.
    bin_width_beta = (df['BETA'].max() - df['BETA'].min()) / 50
    scaling_factor_beta = len(df['BETA']) * bin_width_beta
    axes[0].plot(x_vals_beta, y_vals_beta * scaling_factor_beta, 'k--', label='Normal Fit')
    axes[0].legend(fontsize=12)

    # Plot the distribution of Z scores on the right.
    axes[1].hist(df['Z'], bins=50, color=one_color, alpha=0.7, edgecolor='black', density=False)
    axes[1].set_title(f"Distribution of Z Scores\nin {names_dict[name]}", fontsize=20)
    axes[1].set_xlabel("Z Score", fontsize=14)
    axes[1].set_ylabel("Density", fontsize=14)

    # Generate the standard normal curve (mean=0, std=1) for Z scores.
    x_vals_z = np.linspace(-5, 5, 1000)
    y_vals_z = (1 / np.sqrt(2 * np.pi)) * np.exp(-0.5 * x_vals_z ** 2)
    bin_width_z = (df['Z'].max() - df['Z'].min()) / 50
    scaling_factor_z = len(df['Z']) * bin_width_z
    axes[1].plot(x_vals_z, y_vals_z * scaling_factor_z, 'k--', label='Standard Normal')
    axes[1].legend(fontsize=12)

    plt.tight_layout()

    return fig


# Function to harmonize alleles
def MatchAlleles(row):
    if row['A1'] == row['A1_ref'] and row['A2'] == row['A2_ref']:
        return row
    elif row['A1'] == row['A2_ref'] and row['A2'] == row['A1_ref']:
        row = row.copy()
        row['BETA'] = -row['BETA']
        row['Z'] = -row['Z']
        row['A1'], row['A2'] = row['A1_ref'], row['A2_ref']
        return row
    else:
        return None # if the alleles don't match the reference

##################


# Load the two datasets (munged and formatted summary statistics).
Path_OG = sys.argv[1].rstrip()
print(f"Original File is {Path_OG}.")

Path_Munge = sys.argv[2].rstrip()
print(f"Munged File is {Path_Munge}")

OutDir = sys.argv[3].rstrip()
print(f"Refernce File is: {OutDir}")

filename, name = GetInfo(Path_OG)
print (f"Working on phenotype {name}, in file {filename}")

# Log file
SetUpLogging(OutDir, name)

# Import dataframes
df_og = ReadInputFile(Path_OG)
df_munge = ReadInputFile(Path_Munge)

# Allele reference file
REF_ALLELE_FILE = os.environ.get(
    "CVP_1KG_SNPA1A2_REF",
    "/path/to/reference/1kgPhase3_SNPA1A2.ref.gz"
)
ref = ReadInputFile(REF_ALLELE_FILE)

# Colors
one_color = '#a67664'
two_colors = ['#59433f', '#a38675']

# Filter original dataframe for SNPs in the munged dataframe
filtered = GetFiltered(df_og, df_munge)
print(f"Dataframe after filtering for munged SNPs: {filtered.shape}, with columns: {filtered.columns}")

# Merge on the SNP column with the reference file, and add the A1 and A2 of the reference file
to_compare = pd.merge(filtered, ref, on='SNP', suffixes=('', '_ref'))
print(f"Columns after preparing to compare: {to_compare.columns}")



# Make sure all allleles match the reference (A1, A2)
matched = to_compare.apply(MatchAlleles, axis=1)
print(f"matched A1 and A2 to the reference.")
matched = matched.dropna(how='any')
print(f"Dataframe after matching alleles with reference: {matched.shape}")

# Create p-value distribution figure
pval_fig = PvalDist(matched, name, NamesDict)
print("Pvalue distribution created successfully.")

# Create Manhattan plot figure
manh_fig = Manhattan(matched, name, NamesDict, threshold=-np.log10(5e-8))
print("Manhattan plot created successfully.")

# Create distributions of the effect sizes
sizes_fig = EffectSizeDist(matched, name, NamesDict)
print("Effect size distributions created successfully.")

# Save
matched.to_csv(os.path.join(OutDir, f'{name}-filtered.tsv.gz'), sep='\t', index=False, compression='gzip')

pval_fig.savefig(f"{OutDir}/{name}-Pdist.png", dpi=300) # dpi = dots per inch, 300 is average-good
manh_fig.savefig(f"{OutDir}/{name}-Manhattan.png", dpi=300)
sizes_fig.savefig(f"{OutDir}/{name}-EffectsDist.png", dpi=300)
