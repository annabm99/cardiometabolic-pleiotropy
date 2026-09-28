import pandas as pd
import sys
import os


PROJECT_DIR = os.environ.get(
    "CVP_PROJECT_DIR",
    "/path/to/cardiovascular_pleiotropies"
)

input_file = os.environ.get(
    "CVP_HERITABILITY_DF",
    f"{PROJECT_DIR}/1-LDSC/1-Heritability/JointData/HeritabilitiesDf.tsv"
)
output_dir = os.environ.get(
    "CVP_FINAL_TABLES_DIR",
    f"{PROJECT_DIR}/FinalTables"
)


# Load dataframe

df = pd.read_csv(input_file, sep='\t')


# =========================
# CSV export
# =========================

csv_path = f'{output_dir}/Table2_HeritabilityEstimates.csv'


df.to_csv(csv_path, index=False)

print(f'Saved CSV table: {csv_path}')

# =========================
# LaTeX export
# =========================

latex_path = f'{output_dir}/Table2_HeritabilityEstimates.tex'

latex_table = df.to_latex(
    index=False,
    escape=False,
    float_format='%.4f',
    caption='SNP-based heritability estimates across analyzed phenotypes.',
    label='tab:heritability_estimates'
)

with open(latex_path, 'w') as f:
    f.write(latex_table)

print(f'Saved LaTeX table: {latex_path}')
