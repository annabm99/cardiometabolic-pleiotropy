import glob
import re
import pandas as pd # type: ignore
import numpy as np # type: ignore
import sys
import logging
import matplotlib.pyplot as plt # type:ignore
import seaborn as sns # type: ignore
from sklearn.manifold import MDS # type: ignore
from matplotlib.lines import Line2D # type: ignore

##### FUNCTIONS

def SetUpLogging(out_dir):
    """
    Initializes logging and redirects stdout and stderr to a log file in the specified output directory.

    Parameters:
        out_dir (str): Output directory for the log file.
    
    Output: none
    """
    log_file_path = f'{out_dir}/PairsDf-log.txt'
    logging.basicConfig(
        filename=log_file_path,
        filemode='w',  # Overwrites the file each time the script is run
        level=logging.INFO,  # Captures all info-level and higher messages
        format='%(asctime)s - %(levelname)s - %(message)s'  # Log message format
    )
    sys.stdout = open(log_file_path, 'a')  # Append print statements to the log
    sys.stderr = open(log_file_path, 'a')  # Append warnings/errors to the log


def ParseLogFile(filepath):
    """
    Parse a single LDSC log file to extract:
    - Phenotype pair (from the file name, e.g., "CAD_d_vs_BMI_t")
    - Genetic correlation (rg), its standard error (se) and p-value
    
    Parameters:
        filepath (str): The path to the LDSC log file.
    
    Returns:
        dict: Dictionary with keys 'Phenotype1', 'Phenotype2', 'rg', 'se', and 'p-value'

    """
    # Open file
    with open(filepath, 'r') as f:
        text = f.read()
    
    # Extract the genetic correlation and its standard error.
    rg_match = re.search(r'Genetic Correlation:\s*([\d\.\-eE]+)\s*\(([\d\.\-eE]+)\)', text)
    if not rg_match:
        print(f"Failed to extract genetic correlation from file: {filepath}. Exiting.")
        sys.exit(1)
    rg = float(rg_match.group(1))
    se = float(rg_match.group(2))

    # Extract the p-value.
    p_match = re.search(r'P:\s*([\d\.\-eE]+)', text)

    if not p_match:
        print(f"Failed to extract p-value from file: {filepath}. Exiting.")
        sys.exit(1)
    p_value = float(p_match.group(1))
    
    # Extract phenotype names from the file name.
    # Assuming the file name follows the pattern: "Phenotype1_vs_Phenotype2-GenCorr.log"
    filename = filepath.split('/')[-1]  # get the file name from the full path
    base = filename.replace('-GenCorr.log', '')
    if '_vs_' in base:
        phenotype1, phenotype2 = base.split('_vs_')
    else:
        print(f"Failed to extract phenotype names from filename: {filename}. Only one phenotype detected. Exiting.")
        sys.exit(1)
    
    return {
        'Phenotype1': phenotype1,
        'Phenotype2': phenotype2,
        'rg': rg,
        'se': se,
        'p-value': p_value
    }

def CreateGeneticCorrelationDataframe(directory_path):
    """
    Process all LDSC log files in the specified directory and create a DataFrame with the parsed data:
        - Phenotype 1 code
        - Phenotype 2 code
        - Genetic Correlation (rg)
        - Standard Error of the Genetic Correlation (se)
        - P-value of the Genetic Correlation (p_value).

    Parameters:
        directory_path (str): The path to the directory containing LDSC log files.

    Returns:
        pd.Dataframe: A dataframe containing genetic correlation data between all phenotype pairs.
    """
    # Get all log files in the correlations directory
    log_files = glob.glob(f"{directory_path}/*.log")
    if not log_files:
        print("No log files found in the specified directory. Exiting.")
        sys.exit(1)

    # Get the records for each pair
    records = [ParseLogFile(file) for file in log_files]

    # Return in dataframe format
    return pd.DataFrame(records)

def SaveDf(df, output_dir):
    """
    Saves a dataframe to an output file to a gzipped CSV, and tab-separated.

    Parameters:
        df (pd.Dataframe): The dataframe to be saved.
        output_dir (str): Full path to the output directory.

    Returns:
        str: The path to the output file.
    """

    # Create path to the output file
    outfile = str(output_dir + "/PairsDf.gz")

    # Save in df format
    df.to_csv( outfile, sep='\t', index_label=None, index=False, decimal='.', compression='gzip')

    # Return the path to the output file
    return outfile

def GenerateVisuals(df, output_dir):
    """
    Generates a heatmap and MDS (multi-dimensional scaling) plot from the genetic corerlations matrix.
    
    Input:
        df (pd.Dataframe): DataFrame containing 'Phenotype1', 'Phenotype2', 'rg', 'se', 'p-value'
        output_dir (str): directory where the plot figures will be saved.

    Output:
        str: path to the output files.
    """
    # Dictionary to translate from codes to full names
    translation_dict = {
        'CAD_d': 'Coronary Artery Disease',
        'T2D_d': 'Type 2 Diabetes',
        'HT_d': 'Hypertension',
        'STR_d': 'Stroke',
        'BMI_t': 'Body Mass Index (BMI)',
        'WC_t': 'Waist Circumference',
        'SBP_t': 'Systolic Blood Pressure',
        'DBP_t': 'Diastolic Blood Pressure',
        'TGL_t': 'Triglycerides',
        'LDL_t': 'LDL Cholesterol',
        'HDL_t': 'HDL Cholesterol',
        'FG_t': 'Random Glucose'
    }
    # List of diseases in the data
    diseases = [
        "Coronary Artery Disease",
        "Type 2 Diabetes",
        "Hypertension",
        "Stroke"
    ]

    # Create empty correlation, SE and p-value matrices
    rg_matrix = pd.DataFrame(np.nan, index=translation_dict.keys(), columns=translation_dict.keys())
    se_matrix = pd.DataFrame(np.nan, index=translation_dict.keys(), columns=translation_dict.keys())
    pval_matrix = pd.DataFrame(np.nan, index=translation_dict.keys(), columns=translation_dict.keys())
    
    # Fill the matrices symmetrically
    for _,row in df.iterrows():
        p1, p2 = row['Phenotype1'], row['Phenotype2']
        rg_matrix.loc[p1, p2] = rg_matrix.loc[p2, p1] = row['rg']
        se_matrix.loc[p1, p2] = se_matrix.loc[p2, p1] = row['se']
        pval_matrix.loc[p1, p2] = pval_matrix.loc[p2, p1] = row['p-value']

    # Fill diagonals
    np.fill_diagonal(rg_matrix.values, 1.0)
    np.fill_diagonal(se_matrix.values, 0.0)
    np.fill_diagonal(pval_matrix.values, 0.0)

    # Rename rows and columns to the full names for plotting
    rg_matrix = rg_matrix.rename(index=translation_dict, columns=translation_dict)

    # Ensure row/column order matches to maintain symmetry
    rg_matrix = rg_matrix.loc[rg_matrix.columns]
    
    # Create heatmap
    # Mask to show only the values in the upper diagonal
    mask = np.triu(np.ones_like(rg_matrix, dtype=bool), k=1)
    marker_dict = {pheno: "\u25cf " + pheno if pheno in diseases else "\u25b2" + pheno for pheno in rg_matrix.index}
    heatmap_matrix = rg_matrix.rename(index=marker_dict, columns=marker_dict)

    plt.figure(figsize=(12,10))
    sns.heatmap(
        heatmap_matrix.astype(float),
        mask=mask,
        cmap='RdBu_r',
        center=0,
        annot=True,
        fmt='.2f',
        annot_kws={"size": 14},
        linewidths=0.5,
        cbar_kws = {"label": "Genetic Correlation"}
    )
    cbar = plt.gca().collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label("Genetic Correlation", size=16)

    legend_elements = [
        Line2D([0], [0], marker='o', color='black', markersize=8, linestyle='None', label='Diseases'),
        Line2D([0], [0], marker='^', color='black', markersize=8, linestyle='None', label='Traits')
    ]

    plt.legend(handles=legend_elements, loc='upper right', fontsize=12)
    plt.title("Genetic Correlation Matrix", size = 20)
    plt.xticks(rotation=45, ha='right', fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    heatmap_path = (f"{output_dir}/rg_heatmap.png")
    plt.savefig(heatmap_path)
    plt.close()

    # --- Heatmap: Y-axis diseases only ---
    available_diseases = [d for d in diseases if d in rg_matrix.index]
    rg_y_diseases = rg_matrix.loc[available_diseases, :]

    marker_dict = {
        pheno: "\u25cf " + pheno if pheno in diseases else "\u25b2 " + pheno
        for pheno in rg_matrix.columns
    }

    rg_y_diseases_marked = rg_y_diseases.rename(
        index=marker_dict,
        columns=marker_dict
    )

    plt.figure(figsize=(12, 5))
    sns.heatmap(
        rg_y_diseases_marked.astype(float),
        cmap='RdBu_r',
        center=0,
        vmin=-1,
        vmax=1,
        annot=True,
        fmt='.2f',
        annot_kws={"size": 14},
        linewidths=0.5,
        cbar_kws={"label": "Genetic Correlation"}
    )

    cbar = plt.gca().collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label("Genetic Correlation", size=16)

    plt.xticks(rotation=45, ha='right', fontsize=14)
    plt.yticks(fontsize=14)

    plt.title("Genetic Correlations", fontsize=18)

    plt.tight_layout()
    y_diseases_path = f"{output_dir}/rg_heatmap_y_diseases.png"
    plt.savefig(y_diseases_path)
    plt.close()

    # Create MDS plot
    distance_matrix = 1 - rg_matrix.values.astype(float)
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    coords = mds.fit_transform(distance_matrix)
    df_mds=pd.DataFrame(coords, columns = ['Dim1', 'Dim2'], index=rg_matrix.index)

    marker_map = {pheno: 'o' if pheno in diseases else '^' for pheno in df_mds.index}
    
    fig, ax = plt.subplots(figsize=(8,6))

    for pheno in df_mds.index:
        x, y = df_mds.loc[pheno, 'Dim1'], df_mds.loc[pheno, 'Dim2']
        marker = marker_map.get(pheno, 'o')
        ax.scatter(x, y, s=100, marker=marker, color='black', alpha=0.5)

        if pheno == "Waist Circumference":
            ax.text(x - 0.03, y + 0.05, pheno, fontsize=12, ha='right', va='top')
        else:
            ax.text(x + 0.02, y + 0.01, pheno, fontsize=12)

    ax.legend(handles=legend_elements, loc='upper right', fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title("Genetic Distances Between Phenotypes: MDS Plot", fontsize=16)
    ax.set_xlabel("Dimension 1", fontsize=14)
    ax.set_ylabel("Dimension 2", fontsize=14)
    plt.tight_layout()

    mds_path = f"{output_dir}/mds_plot.png"
    plt.savefig(mds_path)
    plt.close()

    return heatmap_path, y_diseases_path, x_diseases_path, mds_path

##### MAIN SCRIPT

in_directory_path = sys.argv[1].rstrip()
out_directory_path = sys.argv[2].rstrip()
SetUpLogging(out_directory_path)

print(f"Retrieving all log files from {in_directory_path}")
df = CreateGeneticCorrelationDataframe(in_directory_path)

print(f"DataFrame created with shape: {df.shape}")
print(df.head(5))

output_file = SaveDf(df, out_directory_path)
print(f"DataFrame saved to: {output_file}")
heatmap_path, y_diseases_path, x_diseases_path, mds_path = GenerateVisuals(
    df,
    out_directory_path
)

print("Visuals saved to:")
print(f"  Full heatmap: {heatmap_path}")
print(f"  Y-axis diseases heatmap: {y_diseases_path}")
print(f"  X-axis diseases heatmap: {x_diseases_path}")
print(f"  MDS plot: {mds_path}")
