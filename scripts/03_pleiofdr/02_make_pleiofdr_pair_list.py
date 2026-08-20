# File paths
# Genome-wide PleioFDR analyses retain HDL and LDL, but trait-trait
# comparisons are excluded.
import sys

list_file = sys.argv[1] if len(sys.argv) > 1 else "/path/to/2-PleioFDR/2-PairedAnalysis/Input.txt"
out_file = sys.argv[2] if len(sys.argv) > 2 else "/path/to/2-PleioFDR/2-PairedAnalysis/Pairs.txt"

# Read phenotypes
with open(list_file) as f:
    phenotypes = [line.strip() for line in f if line.strip()]

pairs = []

# Generate pairs according to rules
for i, pheno1 in enumerate(phenotypes):
    for pheno2 in phenotypes[i+1:]:
        if pheno1.endswith('_t.mat') and pheno2.endswith('_t.mat'):
            continue  # skip t-t pairs
        pairs.append((pheno1, pheno2))

# Write pairs to file
with open(out_file, 'w') as f:
    for p1, p2 in pairs:
        f.write(f"{p1},{p2}\n")

print(f"Generated {len(pairs)} pairs.")
