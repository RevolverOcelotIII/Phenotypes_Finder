import sys
sys.path.append('/home/domdeny/src/bioinfo/genome-api/pipeline/Phenotypes_Finder/implementation')
import os
import subprocess

# Alinhar os consensos com a sequência de referência

RESULT_FOLDER = "/home/domdeny/src/bioinfo/genome-api/pipeline/Phenotypes_Finder/result/Trimmomatic/"

def get_ref_seq(file_name):
    if "RHD" in file_name:
        return f"/home/domdeny/src/bioinfo/genome-api/pipeline/Phenotypes_Finder/implementation/RefSeq/RefSeq_RHD.fasta"
    elif "RHCE" in file_name:
        return f"/home/domdeny/src/bioinfo/genome-api/pipeline/Phenotypes_Finder/implementation/RefSeq/RefSeq_RHCE.fasta"
    
def single_alignment():
    consensus_list = os.listdir(RESULT_FOLDER)
    for consensus in consensus_list:
        consensus_file = [f for f in os.listdir(f"{RESULT_FOLDER}{consensus}") if f.endswith("_new_consensus.fasta")]
        ref_seq_file = get_ref_seq(consensus)
        print(f"cat {ref_seq_file} \'\\n\' {RESULT_FOLDER}{consensus}/{consensus_file[0]} > {RESULT_FOLDER}{consensus}/{consensus}_merged.fasta")
        subprocess.getoutput(f"cat {ref_seq_file} \'\n\' {RESULT_FOLDER}{consensus}/{consensus_file[0]} > {RESULT_FOLDER}{consensus}/{consensus}_merged.fasta")
        command = f"mafft --auto {RESULT_FOLDER}{consensus}/{consensus}_merged.fasta > {RESULT_FOLDER}{consensus}/{consensus}_mafft.fasta"
        print("Aligning with MAFFT: " + command)
        subprocess.call(command, shell=True)

single_alignment()