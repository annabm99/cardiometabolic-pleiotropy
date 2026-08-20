import sys
import os
import pandas as pd
import numpy as np

#### FUNCTIONS

def GetPhenotypeNames(filepath):
    base = os.path.basename(filepath)
    names = base.split('-')[0]
    parts = names.split('_')
    phen1 = "_".join(parts[:2])
    phen2 = "_".join(parts[2:4])
    return phen1, phen2


def SplitFile(pleio_df, phen1, phen2):
    z1_col = f'Z-{phen1}'
    z2_col = f'Z-{phen2}'
    if z1_col not in pleio_df.columns or z2_col not in pleio_df.columns:
        raise ValueError(f"Expected columns {z1_col} and {z2_col} not found.")

    pleio_df['Z_product'] = pleio_df[z1_col] * pleio_df[z2_col]
    pleio_df["Pleiotropy"] = np.where(pleio_df["Z_product"] > 0, "Positive",
                                 np.where(pleio_df["Z_product"] < 0, "Negative", "Zero"))

    # Positive pleiotropies
    pos = pleio_df[pleio_df["Pleiotropy"] == "Positive"]
    neg = pleio_df[pleio_df["Pleiotropy"] == "Negative"]

    return pos, neg

#### MAIN SCRIPT
PleioFile = sys.argv[1].rstrip()
OutDir = sys.argv[2].rstrip()

phen1, phen2 = GetPhenotypeNames(PleioFile)


pleio_df = pd.read_csv(PleioFile, sep = "\t", low_memory=False)

print(f"Dataframe looks like this: {pleio_df.head(n=3)}")

pos_pleio, neg_pleio = SplitFile(pleio_df, phen1, phen2)

basename = os.path.basename(PleioFile).replace("-SNPInfo.txt.gz", "").replace("-LociInfo.txt.gz", "")
pos_path = os.path.join(OutDir, f"{basename}_positive.csv.gz")
neg_path = os.path.join(OutDir, f"{basename}_negative.csv.gz")
pos_pleio.to_csv(pos_path, sep='\t', index=False, compression='gzip')
neg_pleio.to_csv(neg_path, sep='\t', index=False, compression='gzip')

print(f"Saved positive pleiotropies to: {pos_path}")
print(f"Saved negative pleiotropies to: {neg_path}")