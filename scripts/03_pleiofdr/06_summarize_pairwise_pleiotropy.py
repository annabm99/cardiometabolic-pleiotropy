import sys
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.patches import Patch
import seaborn as sns
from scipy.stats import mannwhitneyu

#######
# UTILITY FUNCTIONS
#######

def GetPhenotypeNames(filepath):
  base = os.path.basename(filepath)
  names = base.split('-')[0]
  parts = names.split('_')
  phen1 = "_".join(parts[:2])
  phen2 = "_".join(parts[2:4])
  return phen1, phen2

def ProcessFile(file_path):
  df = pd.read_csv(file_path, compression='gzip', sep='\t', low_memory=False)
  phen1, phen2 = GetPhenotypeNames(file_path)

  z1_col = f'Z-{phen1}'
  z2_col = f'Z-{phen2}'

  # Check columns exist
  if z1_col not in df.columns or z2_col not in df.columns:
    raise ValueError(f"Expected columns {z1_col} and {z2_col} not found in {file_path}.")

  # Incorporate column with the z product
  df['Z_product'] = df[z1_col] * df[z2_col]
  
  # Classify pleiotropy direction
  df["Pleiotropy"] = np.where(
        df["Z_product"] > 0, "Positive",
        np.where(df["Z_product"] < 0, "Negative", "Zero")
    )
  
  # Per-SNP dataframe
  per_snp = df[['SNP' ,z1_col, z2_col, "Z_product", "Pleiotropy"]].copy()
  per_snp = per_snp.rename(columns={
        z1_col: 'Z_Phenotype1',
        z2_col: 'Z_Phenotype2'
    })
  per_snp["Phenotype1"] = phen1
  per_snp["Phenotype2"] = phen2
  
  # Summary
  summary = {
        "Phenotype1": phen1,
        "Phenotype2": phen2,
        "Count_Positive": (df["Pleiotropy"] == "Positive").sum(),
        "Count_Negative": (df["Pleiotropy"] == "Negative").sum(),
    }
  return summary, per_snp

def CreatePleiotropySummaryDf(file_list):
  records = []
  snp_records = []
  for file_path in file_list:
    summary, per_snp = ProcessFile(file_path)
    if summary is not None:
      records.append(summary)
      snp_records.append(per_snp)

    df_summary = pd.DataFrame(records)
    df_snps = pd.concat(snp_records, ignore_index=True)
  
  return df_summary, df_snps

########
# PLOTTING FUNCTIONS
########

def TraitDiseaseBarplot(dfs, pos_color, neg_color, translation_dict, traits_order):

    legend_elements = [
       Patch(facecolor=neg_color, label = 'Negative Pleiotropies (opposite direction)'),
       Patch (facecolor=pos_color, label = 'Positive pleiotropies (same direction)')
    ]

    # Compute max values for each DataFrame
    max_values = np.array([df[['Count_Negative', 'Count_Positive']].sum(axis=1).max() for df in dfs])

    # Get sorting indices based on max_values (sorted from smallest to largest)
    sorted_indices = np.argsort(max_values)

    # Sort max_values using the indices
    sorted_max_values = max_values[sorted_indices]

    # Sort dataframes using the same indices
    sorted_dataframes = [dfs[i] for i in sorted_indices]

    y_pos = range(len(traits_order))

    def format_trait_name(trait_name):
        words = trait_name.split()
        if len(words) > 2:
            # Split into two lines
            return '\n'.join([' '.join(words[:len(words)//2]), ' '.join(words[len(words)//2:])])
        else:
            # Return as it is if only one or two words
            return trait_name

    fig, axes = plt.subplots(
    1, len(sorted_dataframes),
    figsize=(3*len(sorted_dataframes), 6),
    sharey=True,
    gridspec_kw={'width_ratios': sorted_max_values/5}
    )
    # Adjust space between subplots
    fig.subplots_adjust(wspace=0.5)

    # Iterate over the sorted dataframes, axes, and max values
    for i, (df_subset, ax, max_x) in enumerate(zip(sorted_dataframes, axes.flatten(), sorted_max_values)):
        # Extract the disease code from the dataframe
        disease_code = df_subset.iloc[0, 1]  # Assuming disease code is in column 2 (index 1)

        # Access the full name from the names dictionary
        disease_name = translation_dict.get(disease_code, disease_code)  # Default to disease_code if not found in dictionary

        # Plot negative pleiotropy counts
        ax.barh(y_pos, df_subset['Count_Negative'], color=neg_color, label='Negative pleiotropies')

        # Plot positive pleiotropy counts (stacked on top)
        ax.barh(y_pos, df_subset['Count_Positive'], left=df_subset['Count_Negative'], color=pos_color, label='Positive pleiotropies')

        # Set y-axis labels
        labels = [format_trait_name(translation_dict.get(trait, trait)) for trait in traits_order]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=13)

        # Set title
        if disease_code == 'CAD_d':
            ax.set_title("Coronary\nArtery\nDisease", fontsize=13, pad=10)
        else:
            ax.set_title(disease_name, fontsize=13, pad=10)  # Use full name from the dictionary

        # Draw a vertical reference line at x=0
        ax.axvline(0, color='black', linewidth=1)

        # Add gridlines for readability
        ax.xaxis.grid(True, linestyle='--', alpha=0.6)

        # Remove top and right borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.margins(x=0.05)

        # Make axis the same as in the df
        ax.invert_yaxis()

    # Set the x-axis label once for the entire figure
    fig.text(0.5, -0.01, "Pleiotropy Count", ha='center', fontsize=13)

    # Improve legend placement
    fig.legend(handles=legend_elements, fontsize=13, bbox_to_anchor=(0.8, -0.02), ncol=2)

    # Add a main title
    fig.suptitle('Pleiotropies between Diseases and Traits', fontsize=15, fontweight='bold')

    # Adjust layout and display the plot
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    return fig

def DiseaseDiseaseBarplot(df, pos_color, neg_color, translation_dict, diseases):

    def format_dis_name(dis_name):
        words = dis_name.split()
        if len(words) > 2:
            return '\n'.join([' '.join(words[:len(words)//2]), ' '.join(words[len(words)//2:])])
        else:
            return dis_name
    
    # Create a figure with multiple subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    for i, (disease, ax) in enumerate(zip(diseases, axes.flatten())):
        # Filter for the disease
        df_disease = df[(df['Phenotype1'] == disease) | (df['Phenotype2'] == disease)].copy()

        # Swap axes to ensure 'disease' is always in Phenotype1
        for index, row in df_disease.iterrows():
            if row['Phenotype1'] != disease and row['Phenotype2'] == disease:
                df_disease.loc[index, ['Phenotype1', 'Phenotype2']] = row[['Phenotype2', 'Phenotype1']].values

        # Get sorted traits
        other_traits = df_disease['Phenotype2'].unique().tolist()
        traits_order = sorted(other_traits)
        df_disease = df_disease.set_index('Phenotype2').reindex(traits_order).reset_index()

        # Define y positions for horizontal bars
        y_pos = np.arange(len(df_disease))

        # Plot the bars as horizontal bars
        ax.barh(y_pos, df_disease['Count_Negative'], color=neg_color, label='Negative Pleiotropy')
        ax.barh(y_pos, df_disease['Count_Positive'], left=df_disease['Count_Negative'], color=pos_color, label='Positive Pleiotropy')

        # Set y labels
        labels = [format_dis_name(translation_dict.get(trait, trait)) for trait in traits_order]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=15)

        ax.set_xlabel("Pleiotropy Count", fontsize=15)
        ax.set_title(f"{translation_dict.get(disease, disease)}", fontsize=15, pad=20)

        # Hide legend for all but the first plot
        if i != 0:
            ax.legend().set_visible(False)

    # Add a single legend outside the loop
    fig.legend(["Negative Pleiotropy", "Positive Pleiotropy"], loc='upper center', bbox_to_anchor=(0.5, -0.01), ncol=2, fontsize=15)

    # Add a main title
    fig.suptitle('Pleiotropies between Diseases', fontsize=17, fontweight='bold')

    # Final layout adjustments
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    return fig

def PlotAndTestEffectSizes(df_snps, outdir, disease_list, positive_color, negative_color, translation_dict):
    results = []

    for disease in disease_list:
            # Select SNPs where this disease is Phenotype1 or Phenotype2
            subset = df_snps[(df_snps["Phenotype1"] == disease) | (df_snps["Phenotype2"] == disease)].copy()
            if subset.empty:
                continue

            # Extract the Z-score corresponding to the disease
            subset["Z_disease"] = subset.apply(
                lambda row: row["Z_Phenotype1"] if row["Phenotype1"] == disease else row["Z_Phenotype2"], axis=1
            )

            # Separate positive vs negative pleiotropies
            #pos = subset.loc[subset["Pleiotropy"] == "Positive", "Z_disease"]
            #neg = subset.loc[subset["Pleiotropy"] == "Negative", "Z_disease"]

            subset["Z_disease_abs"] = subset["Z_disease"].abs()
            pos = subset.loc[subset["Pleiotropy"] == "Positive", "Z_disease_abs"]
            neg = subset.loc[subset["Pleiotropy"] == "Negative", "Z_disease_abs"]


            if len(pos) > 0 and len(neg) > 0:
                # Mann–Whitney test
                stat, pval = mannwhitneyu(pos, neg, alternative="two-sided")
                
                # Determine significance
                alpha = 0.05
                sig = "Yes" if pval < alpha else "No"

                # Determine which group has higher median
                median_pos = pos.median()
                median_neg = neg.median()
                direction = "Positive > Negative" if median_pos > median_neg else "Negative > Positive"

                # Get sample sizes
                n_pos = len(pos)
                n_neg = len(neg)

                results.append({
                    "Disease": disease,
                    "U_stat": stat,
                    "pval": pval,
                    "Significant": sig,
                    "N_Positive": n_pos,
                    "N_Negative": n_neg,
                    "Positive median": median_pos,
                    "Negative median": median_neg,
                    "Higher": direction
                })

             # Violin plot with dots
            plt.figure(figsize=(8, 6))
            sns.violinplot(
                data=subset,
                x="Pleiotropy",
                y="Z_disease_abs",
                order=["Negative", "Positive"],
                inner=None,
                palette={"Negative": negative_color, "Positive": positive_color}
            )
            sns.stripplot(
                data=subset,
                x="Pleiotropy",
                y="Z_disease_abs",
                order=["Negative", "Positive"],
                hue="Pleiotropy",
                dodge=False,
                palette={"Negative": negative_color, "Positive": positive_color},
                size=4,
                jitter=True,
                alpha=0.7,
                marker="o"
            )

            # Add median lines
            medians = subset.groupby("Pleiotropy")["Z_disease_abs"].median()
            for i, pleio in enumerate(["Negative", "Positive"]):
                plt.hlines(
                    y=medians[pleio],
                    xmin=i-0.3,
                    xmax=i+0.3,
                    colors='black',
                    linestyles='dashed',
                    linewidth=2
                )

            # Annotate sample sizes above violins
            counts = subset["Pleiotropy"].value_counts()
            for i, pleio in enumerate(["Negative", "Positive"]):
                n = counts.get(pleio, 0)
                plt.text(i, subset["Z_disease_abs"].max() + 0.1, f"n={n}", ha='center', fontsize=12)

            plt.title(f"{translation_dict.get(disease, disease)} Absolute Z-scores by Pleiotropy Type\n(p={pval:.2e})")
            plt.ylabel("Absolute Z-score")
            plt.xlabel("")
            plt.legend([],[], frameon=False)  # hide duplicate legend
            plt.tight_layout()
            plt.savefig(f"{outdir}/{disease}_EffectSizeViolin.png")
            plt.close()

    return pd.DataFrame(results)

######
# DICTIONARIES AND CONSTANTS
######

# Dictionary to get the full phenotype names
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
    'FG_t': 'Fasting Glucose'
}

traits_order = ['BMI_t', 'WC_t', 'SBP_t', 'DBP_t', 'TGL_t', 'LDL_t', 'HDL_t', 'FG_t']
diseases = ['CAD_d', 'T2D_d', 'HT_d', 'STR_d']

positive_color = '#c13639'
negative_color = '#377eb8'

########
# MAIN
########

FileListPath = sys.argv[1].rstrip()

with open(FileListPath, "r") as f:
    FileList = [line.strip() for line in f if line.strip()]  # remove empty lines

print(f"Found {len(FileList)} input files")

OutDir = sys.argv[2].rstrip()
print(f"Output directory is: {OutDir}")

# Generate datasets
df_summary, df_snps = CreatePleiotropySummaryDf(FileList)

# Save outputs
df_summary.to_csv(f"{OutDir}/PleiotropiesSummary.csv.gz", sep='\t', index=False, compression='gzip')
df_snps.to_csv(f"{OutDir}/PleiotropyPerSNP.csv.gz", sep='\t', index=False, compression='gzip')

# Disease-trait plots
dis_trait_dfs = []
for disease in diseases:
    # STEP 1. Filter for the current disease
    df_disease = df_summary[(df_summary['Phenotype1'] == disease) | (df_summary['Phenotype2'] == disease)].copy()

    # STEP 2. Swap Phenotype1 and Phenotype2 so that disease is always in Phenotype1, vectorized swapping
    mask = df_disease['Phenotype2'] == disease
    df_disease.loc[mask, ['Phenotype1', 'Phenotype2']] = df_disease.loc[mask, ['Phenotype2', 'Phenotype1']].values

    # STEP 3. Reindex by traits_order
    df_disease = df_disease.set_index('Phenotype2').reindex(traits_order).reset_index()

    # Append to list
    dis_trait_dfs.append(df_disease)

dis_trait_plot = TraitDiseaseBarplot(dis_trait_dfs, positive_color, negative_color, translation_dict, traits_order)
# Save plot
dis_trait_plot.savefig(f"{OutDir}/trait_disease_pleiotropies.png", dpi=300, bbox_inches="tight")

# Disease-disease plots
disease_df = df_summary[(df_summary['Phenotype1'].isin(diseases)) & (df_summary['Phenotype2'].isin(diseases))]
dis_dis_plot = DiseaseDiseaseBarplot(disease_df, positive_color, negative_color, translation_dict, diseases)
dis_dis_plot.savefig(f"{OutDir}/disease_disease_pleiotropies.png", dpi=300, bbox_inches="tight")

# Run tests + plots
df_tests = PlotAndTestEffectSizes(df_snps, OutDir, diseases, positive_color, negative_color, translation_dict)

# Save test results