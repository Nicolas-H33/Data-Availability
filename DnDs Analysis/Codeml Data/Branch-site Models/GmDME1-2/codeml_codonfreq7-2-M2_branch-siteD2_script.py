import subprocess
import os

codeml_path = "/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/dNdS analysis/paml/src/codeml"
working_dir = "/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/dNdS analysis/Gm 5 seqs_notrim/Control file8 - branch-siteDME2"
control_file = "codeml_codonfreq7-2-M2_branch-siteD2.ctl"

os.chdir(working_dir)

#Run codeml
try:
    subprocess.run([codeml_path, control_file], check=True)
    print("codeml run complete.")
except subprocess.CalledProcessError as e:
    print(f"codeml failed with error: {e}")