import subprocess
import os
import shutil

#Working directory & paths
working_dir = "/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/dNdS analysis/MCScanX_temp"
os.makedirs(working_dir, exist_ok=True)
os.chdir(working_dir)

protein = os.path.join(working_dir, "Glycine_max_filtered_longest.fa")
gff = os.path.join(working_dir, "Glycine_max_filtered.gff")

for f in ["Glycine_max.blast", "Glycine_max.collinearity"]:
    if os.path.exists(os.path.join(working_dir, f)):
        os.remove(os.path.join(working_dir, f))

#BLASTP database
subprocess.run(["makeblastdb", "-in", protein, "-dbtype", "prot", "-out", "Glycine_max_db"], check=True)

subprocess.run(["blastp", "-query", protein, "-db", "Glycine_max_db", "-out", "Glycine_max.blast", "-evalue", "1e-4", "-outfmt", "6", "-num_threads", "7"], check=True)

#Run MCScanX
if os.path.basename(gff) != "Glycine_max.gff":
    shutil.copy(gff, os.path.join(working_dir, "Glycine_max.gff"))

subprocess.run(["/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/dNdS analysis/MCScanX_temp/MCScanX-master/MCScanX", "Glycine_max"], check=True)