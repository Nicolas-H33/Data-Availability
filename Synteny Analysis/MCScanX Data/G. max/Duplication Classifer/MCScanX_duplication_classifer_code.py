import subprocess
import os

#Working directory & paths
working_dir = "/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/dNdS analysis/MCScanX_temp"
os.chdir(working_dir)

#Run MCScanX
subprocess.run(["duplicate_gene_classifier", "Glycine_max"], check=True)