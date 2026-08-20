# Libraries
import sys
import itertools

# Read inputs
InputFile = sys.argv[1] if len(sys.argv) > 1 else "/path/to/1-LDSC/2-Correlation/Input.txt"
OutputFile = sys.argv[2] if len(sys.argv) > 2 else "/path/to/1-LDSC/2-Correlation/Pairs.txt"

### IN CASE WE JUST WANT TO COMBINE 1 TRAIT TO ALL OTHERS (InputFile should contain the paths of all the traits to combine)
#CombTraitPath = ["/path/to/0-Download/4-Munge/HT_d.sumstats.gz"]

with open(InputFile, 'r') as fh:
    PathList = [line.strip() for line in fh if line.strip()]


# Make pairs of all elements on the list
pairs = itertools.combinations(PathList, 2) # For all traits with all traits
#pairs = itertools.product(CombTraitPath, PathList) # For 1 trait with all traits

# Open filehandle to save output
out = open(OutputFile, "w")

# Iterate through pairs and write to file in the adequate format for LDSC
for p in pairs:
    out.write(p[0]+","+p[1]+"\n")

out.close() # close filehandle
