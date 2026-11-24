import subprocess
import os
import shutil

#Working directory & paths
working_dir = "/Users/nickh245/Commonbean_Glycinemax_synteny"
os.makedirs(working_dir, exist_ok=True)
os.chdir(working_dir)

protein = os.path.join(working_dir, "Pv_Gm.fasta")
gff = os.path.join(working_dir, "Pv_Gm.gff")

makeblastdb_path = "/Users/nickh/anaconda3/envs/SoyPhyloTree/bin/makeblastdb"
blastp_path = "/Users/nickh/anaconda3/envs/SoyPhyloTree/bin/blastp"
mcscanx_path = "/Users/nickh/Library/CloudStorage/OneDrive-Personal/Nick_MS_Thesis copy/Objective 1/dNdS analysis/MCScanX_temp/MCScanX-master/MCScanX"

#Checking & clearing
if not os.path.exists(protein):
    raise FileNotFoundError(f"Protein FASTA not found: {protein}")
if not os.path.exists(gff):
    raise FileNotFoundError(f"GFF file not found: {gff}")

old_files = ["Pv_Gm.blast", "Pv_Gm.collinearity", "Pv_Gm.tandem"]
for f in old_files:
    file_path = os.path.join(working_dir, f)
    if os.path.exists(file_path):
        os.remove(file_path)

#BLASTP database
subprocess.run([makeblastdb_path, "-in", "Pv_Gm.fasta", "-dbtype", "prot", "-blastdb_version", "4", "-out", "Pv_Gm_db"], check=True)

subprocess.run([blastp_path, "-query", "Pv_Gm.fasta", "-db", "Pv_Gm_db", "-out", "Pv_Gm.blast", "-evalue", "1e-4", "-outfmt", "6", "-num_threads", "7"], check=True)

if os.path.basename(gff) != "Pv_Gm.gff":
    shutil.copy(gff, os.path.join(working_dir, "Pv_Gm.gff"))

#Run MCScanX
subprocess.run([mcscanx_path, "Pv_Gm"], check=True)
