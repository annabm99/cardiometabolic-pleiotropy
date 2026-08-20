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
# Excel export
# =========================

excel_path = f'{output_dir}/Table2_HeritabilityEstimates.xlsx'

with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:

    df.to_excel(
        writer,
        sheet_name='Table2_Heritability',
        index=False
    )

    worksheet = writer.sheets['Table2_Heritability']

    # Auto-adjust columns
    for column_cells in worksheet.columns:

        length = max(
            len(str(cell.value)) if cell.value is not None else 0
            for cell in column_cells
        )

        worksheet.column_dimensions[
            column_cells[0].column_letter
        ].width = length + 3

print(f'Saved Excel table: {excel_path}')


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
