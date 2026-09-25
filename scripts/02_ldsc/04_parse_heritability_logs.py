import glob
import re
import pandas as pd
import sys
import logging
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import scipy.stats as stats


##### FUNCTIONS


def SetUpLogging(out_dir):
    """
    Initializes logging and redirects stdout and stderr.
    """

    log_file_path = f'{out_dir}/H2_Df-log.txt'

    logging.basicConfig(
        filename=log_file_path,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    sys.stdout = open(log_file_path, 'a')
    sys.stderr = open(log_file_path, 'a')

def ParseLogFile(filepath):
    """
    Parse LDSC heritability log file.
    """

    with open(filepath, 'r') as f:
        text = f.read()

    # Heritability
    h2_match = re.search(
        r'Total Observed scale h2:\s*([\d\.\-eE]+)\s*\(([\d\.\-eE]+)\)',
        text
    )

    if not h2_match:
        raise ValueError(f'Could not extract h2 from {filepath}')

    h2 = float(h2_match.group(1))
    se = float(h2_match.group(2))

    # Z-score
    z_score = h2 / se

    # P-value
    p_value = 2 * stats.norm.sf(abs(z_score))

    # Mean chi-square
    chi_match = re.search(r'Mean Chi\^2:\s*([\d\.\-eE]+)', text)

    if not chi_match:
        raise ValueError(f'Could not extract Mean Chi2 from {filepath}')

    mean_chi2 = float(chi_match.group(1))

    # Lambda GC
    lambda_match = re.search(r'Lambda GC:\s*([\d\.\-eE]+)', text)

    if not lambda_match:
        raise ValueError(f'Could not extract Lambda GC from {filepath}')

    lambda_gc = float(lambda_match.group(1))

    # LDSC intercept
    intercept_match = re.search(r'Intercept:\s*([\d\.\-eE]+)', text)

    if not intercept_match:
        raise ValueError(f'Could not extract intercept from {filepath}')

    intercept = float(intercept_match.group(1))

    # Phenotype name
    filename = filepath.split('/')[-1]
    phenotype = filename.replace('.log', '')

    return {
        'Phenotype': phenotype,
        'h2': round(h2, 4),
        'SE': round(se, 4),
        'Z_score': round(z_score, 2),
        'P_value': p_value,
        'Lambda_GC': round(lambda_gc, 3),
        'Mean_Chi2': round(mean_chi2, 3),
        'LDSC_Intercept': round(intercept, 3)
    }

def CreateGeneticHeritabilitiesDataframe(directory_path):
    """
    Process all LDSC heritability log files.
    """

    log_files = glob.glob(f"{directory_path}/*.log")

    if not log_files:
        print("No log files found.")
        sys.exit(1)

    records = [ParseLogFile(file) for file in log_files]

    return pd.DataFrame(records)



def SaveDf(df, output_dir):
    """
    Save canonical dataframe.
    """

    outfile = f"{output_dir}/HeritabilitiesDf.tsv"

    df.to_csv(
        outfile,
        sep='\t',
        index=False,
    )

    return outfile

def GenerateVisuals(df, sample_sizes, output_dir):
    """
    Generate heritability barplot.
    """

    fullnames_dict = {
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

    diseases = [
        'Coronary Artery Disease',
        'Type 2 Diabetes',
        'Hypertension',
        'Stroke'
    ]

    df = df.sort_values('h2', ascending=True).reset_index(drop=True)

    sample_sizes_series = df['Phenotype'].map(sample_sizes)
    df['SampleSize'] = sample_sizes_series

    df['Phenotype'] = df['Phenotype'].replace(fullnames_dict)

    new_labels = []

    for phenotype in df['Phenotype']:

        if phenotype in diseases:
            new_labels.append(f'●  {phenotype}')
        else:
            new_labels.append(f'▲  {phenotype}')

    norm = mcolors.Normalize(
        vmin=df['SampleSize'].min(),
        vmax=df['SampleSize'].max()
    )

    cmap = cm.get_cmap('Purples')

    colors = [cmap(norm(size)) for size in df['SampleSize']]

    plt.figure(figsize=(6, 7), constrained_layout=True)

    plt.barh(
        df['Phenotype'],
        df['h2'],
        xerr=df['SE'],
        color=colors,
        edgecolor='black',
        linewidth=1,
        capsize=7
    )

    plt.xticks(fontsize=12)

    plt.yticks(
        ticks=range(len(df)),
        labels=new_labels,
        fontsize=12
    )

    plt.ylabel('Phenotype', fontsize=14, labelpad=20)
    plt.xlabel('Heritability (h^2)', fontsize=14, labelpad=20)

    plt.title(
        'Heritability Estimates with\nStandard Errors and Sample Sizes',
        fontsize=18,
        pad=10
    )

    plt.grid(axis='x', linestyle='--', alpha=0.7)

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=plt.gca())
    cbar.set_label('Sample Size', labelpad=20, fontsize=14)

    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='black', markersize=8,
                   linestyle='', label='Disease'),
        plt.Line2D([0], [0], marker='^', color='black', markersize=8,
                   linestyle='', label='Trait')
    ]

    plt.legend(handles=legend_elements, loc='lower right', fontsize=14)

    fig_path = f'{output_dir}/HeritabilitiesBarplot.png'

    plt.savefig(fig_path, bbox_inches='tight', dpi=300)
    plt.close()

    return fig_path

##### MAIN SCRIPT

in_directory_path = sys.argv[1].rstrip()
out_directory_path = sys.argv[2].rstrip()

sample_sizes = {
    'CAD_d': 86995,
    'T2D_d': 441894,
    'HT_d': 56637,
    'STR_d': 446696,
    'BMI_t': 359983,
    'WC_t': 360564,
    'SBP_t': 340159,
    'DBP_t': 340162,
    'TGL_t': 343991,
    'LDL_t': 343621,
    'HDL_t': 315133,
    'FG_t': 314914
}

SetUpLogging(out_directory_path)

print(f'Retrieving log files from {in_directory_path}')

# Create dataframe

df = CreateGeneticHeritabilitiesDataframe(in_directory_path)

print(f'Dataframe created with shape {df.shape}')
print(df.head())

# Save canonical dataframe

output_file = SaveDf(df, out_directory_path)

print(f'Saved dataframe to: {output_file}')

# Generate figure

output_visuals = GenerateVisuals(df, sample_sizes, out_directory_path)

print(f'Visuals saved to: {output_visuals}')