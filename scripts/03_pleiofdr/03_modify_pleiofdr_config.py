import sys
import gc  # Import garbage collection module
import os

### FUNCTION ####

def ModifyTemplate(ref_file, file1, file2, phen1, phen2, pair_dir, stat_type, fdr_threshold, randprune_n, exclude_chr_pos):
    try:
        # Read the template file
        pleiofdr_dir = os.environ.get("CVP_PLEIOFDR_DIR", "/path/to/pleiofdr")
        default_config = os.environ.get(
            "CVP_PLEIOFDR_CONFIG_TEMPLATE",
            os.path.join(pleiofdr_dir, "config_default.txt")
        )
        with open(default_config, "r") as f:
            template = f.readlines()

        # Create a list to hold modified lines
        modified_lines = []

        # Replace variables with the specific values
        for line in template:
            if line.startswith("reffile="):
                line = "reffile={}\n".format(ref_file)
            elif line.startswith("traitfolder="):
                line = "traitfolder=\n"
            elif line.startswith("traitfile1="):
                line = "traitfile1={}\n".format(file1)
            elif line.startswith("traitname1="):
                line = "traitname1={}\n".format(phen1)
            elif line.startswith("traitfiles="):
                line = "traitfiles={{'{}'}}\n".format(file2)
            elif line.startswith("traitnames="):
                line = "traitnames={{'{}'}}\n".format(phen2)
            elif line.startswith("outputdir="):
                line = "outputdir={}\n".format(pair_dir)
            elif line.startswith("stattype="):
                line = "stattype={}\n".format(stat_type)
            elif line.startswith("fdrthresh="):
                line = "fdrthresh={}\n".format(fdr_threshold)
            elif line.startswith("randprune_n="):
                line = "randprune_n={}\n".format(randprune_n)
            elif line.startswith("# Exclusion regions"):
                if exclude_chr_pos:
                    line += "{}\n".format(exclude_chr_pos)
            modified_lines.append(line)

        # Write the modified content to a new configuration file
        outfile = os.path.join(
            pleiofdr_dir,
            "config-{}_{}.txt".format(phen1, phen2)
        )
        with open(outfile, "w") as f:
            f.writelines(modified_lines)

        print("New config file written: {}".format(outfile))

    except Exception as e:
        print("Error:", e)
        sys.exit(1)

    finally:
        # Clean-up memory by deleting large variables
        del template
        del modified_lines
        # Explicitly run garbage collection
        gc.collect()

#### SCRIPT ####

# Check if the correct number of arguments is provided
if len(sys.argv) != 6:
    print("Usage: python 03_modify_pleiofdr_config.py path1 path2 phen1 phen2 pair_dir")
    sys.exit(1)

# Import variables from input arguments
File1, File2, Phen1, Phen2, Dir = [arg.strip() for arg in sys.argv[1:]]

# Execute modifyer function
try:
    ModifyTemplate(
        ref_file=os.environ.get(
            "CVP_PLEIOFDR_REF_MAT",
            "/path/to/reference/ref9545380_1kgPhase3eur_LDr2p1.mat"
        ),
        file1=File1,
        file2=File2,
        phen1=Phen1,
        phen2=Phen2,
        pair_dir=Dir,
        stat_type="conjfdr",
        fdr_threshold="0.05",
        randprune_n="500",
        exclude_chr_pos=""
    )
except Exception as e:
    print("Error:", e)
    sys.exit(1)

# Clean up after execution
gc.collect()  # Run garbage collection manually, if necessary

print("ConfigFile generated for the pair {} and {}".format(Phen1, Phen2))
