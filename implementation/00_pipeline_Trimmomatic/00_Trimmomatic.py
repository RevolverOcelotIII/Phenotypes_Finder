import sys
sys.path.append('/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/implementation')

import subprocess
import os
# import Gene as Gene

import time

def trimmomatic():
    fastq_gz_folder = "/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/dataset/FASTQ/"

    print("Running Trimmomatic pipeline...\n")

    r1_files = [f for f in os.listdir(fastq_gz_folder) if (f.endswith(".fastq.gz") and "_R1" in f.upper())]
    print(r1_files)
    
    for r1_file in r1_files:
        file_name = r1_file.strip(".fastq.gz")
        ref_seq = get_ref_seq(file_name)
        r2_file = get_r2_file(fastq_gz_folder, file_name)
        r2_file.join(".fastq.gz")
        print(r2_file)
        #trimmed_folder = f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/result/Trimmomatic/trimmed/{file_name}/"
        output_folder = f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/result/Trimmomatic/{file_name}/"
        #if not os.path.exists(trimmed_folder):
        #    os.makedirs(trimmed_folder)
        #if not os.path.exists(output_folder):
        #    os.makedirs(output_folder)

        if len(ref_seq)> 0:
            singleGeneRun(fastq_gz_folder, file_name, r1_file, r2_file, ref_seq)
        else:
            bothGenesRun(fastq_gz_folder, file_name, r1_file, r2_file)

        '''
        output_folder = f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/result/Trimmomatic/{file_name}/"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        step_01(fastq_gz_folder, r1_file, r2_file, output_folder)
        step_02_minimap2(file_name, output_folder, r1_file, r2_file, ref_seq)
        #step_03(ref_seq, output_folder, file_name)
        step_04(file_name, output_folder)
        step_05(ref_seq, output_folder, file_name)
        '''
        

def singleGeneRun(fastq_gz_folder, file_name, r1_file, r2_file, ref_seq):
    trimmed_folder = f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/result/Trimmomatic/trimmed/{file_name}/"
    output_folder = f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/result/Trimmomatic/{file_name}/"
    if not os.path.exists(trimmed_folder):
        os.makedirs(trimmed_folder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    step_01(fastq_gz_folder, r1_file, r2_file, trimmed_folder)
    step_02_minimap2(file_name, output_folder, trimmed_folder, r1_file, r2_file, ref_seq[0])
    #step_03(ref_seq, output_folder, file_name)
    step_04(file_name, output_folder)
    step_05(ref_seq[1], output_folder, file_name)

def bothGenesRun(fastq_gz_folder, file_name, r1_file, r2_file):
    #output_folder = f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/result/Trimmomatic/{file_name}/"
    trimmed_folder = f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/result/Trimmomatic/trimmed/{file_name}/"
    if not os.path.exists(trimmed_folder):
        os.makedirs(trimmed_folder)
    # for both
    step_01(fastq_gz_folder, r1_file, r2_file, trimmed_folder)

    # FOR RHD
    rhd_output_folder = f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/result/Trimmomatic/{file_name}_RHD/"
    if not os.path.exists(rhd_output_folder):
        os.makedirs(rhd_output_folder)
    step_02_minimap2(file_name, rhd_output_folder, trimmed_folder, r1_file, r2_file, get_ref_seq('RHD')[0])
    #step_03(get_ref_seq('RHD'), rhd_output_folder, file_name)
    step_04(file_name, rhd_output_folder)
    step_05(get_ref_seq('RHD')[1], rhd_output_folder, file_name)

    # FOR RHCE
    rhce_output_folder = f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/result/Trimmomatic/{file_name}_RHCE/"
    if not os.path.exists(rhce_output_folder):
        os.makedirs(rhce_output_folder)
    step_02_minimap2(file_name, rhce_output_folder, trimmed_folder, r1_file, r2_file, get_ref_seq('RHCE')[1])
    #step_03(get_ref_seq('RHCE'), rhce_output_folder, file_name)
    step_04(file_name, rhce_output_folder)
    step_05(get_ref_seq('RHCE')[1], rhce_output_folder, file_name)

def get_ref_seq(file_name):
    # Get RefSeq
    if 'RHD' in file_name:
        return [
            f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/implementation/RefSeq/RefSeq_RHD.fasta",
            f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/implementation/RefSeq/RefSeq_RHD_Exome_new.fasta"
        ]
    elif 'RHCE' in file_name:
        return [
            f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/implementation/RefSeq/RefSeq_RHCE.fasta",
            f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/implementation/RefSeq/RefSeq_RHCE_Exome_new.fasta"
        ]
    return []

def get_r2_file(fastq_gz_folder, file_name):
    prefix = file_name.split("_")
    r2_file = [f for f in os.listdir(fastq_gz_folder) if (f.find(prefix[0]) != -1 and "_R2" in f.upper())]
    return r2_file[0]

def step_01(fastq_gz_folder, r1_file, r2_file, output_folder):
    # Execute: TrimmomaticPE 
    # -threads 8 
    # $NAMEFILE"_L001_R1_001.fastq.gz" 
    # $NAMEFILE"_L001_R2_001.fastq.gz" 
    # $NAMEFILE"_L001_R1_trimmed_001.fastq.gz" 
    # $NAMEFILE"_L001_R1_unpaired_001.fastq.gz" 
    # $NAMEFILE"_L001_R2_trimmed_001.fastq.gz" 
    # $NAMEFILE"_L001_R2_unpaired_001.fastq.gz" 
    # ILLUMINACLIP:adapters.fasta:2:30:10:2:true 
    # HEADCROP:15 SLIDINGWINDOW:6:15 MINLEN:50

    r1_stripped = r1_file.strip(".fastq.gz")
    r2_stripped = r2_file.strip(".fastq.gz")
    print(r2_stripped)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    adapters = "/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/implementation/00_pipeline_Trimmomatic/Adapters/adapters.fasta"
    # SE se single-end, PE se paired-end
    command = f"trimmomatic PE -phred33 -threads 8 \
        {fastq_gz_folder}{r1_file} \
        {fastq_gz_folder}{r2_file} \
        {output_folder}{r1_stripped}_trimmed.fastq.gz \
        {output_folder}{r1_stripped}_unpaired.fastq.gz \
        {output_folder}{r2_stripped}_trimmed.fastq.gz \
        {output_folder}{r2_stripped}_unpaired.fastq.gz \
        ILLUMINACLIP:{adapters}:2:30:10:2:true \
        HEADCROP:15 SLIDINGWINDOW:6:15 MINLEN:50"
    print("Step 1: " + command)
    subprocess.call(command, shell=True)

def step_02_spades(file_name, output_folder, r1_file, r2_file):
    # Execute:
    # spades.py 
    # -1 $NAMEFILE"_L001_R1_trimmed_001.fastq.gz" 
    # -2 $NAMEFILE"_L001_R2_trimmed_001.fastq.gz" 
    # -t 8 --only-assembler -k 21,33,55,77 -o assembly
    
    r1_stripped = r1_file.strip(".fastq.gz")
    r2_stripped = r2_file.strip(".fastq.gz")
    if not os.path.exists(f"{output_folder}assembly/tmp/"):
        os.makedirs(f"{output_folder}assembly/tmp/")
    # remover -2 caso não seja paired-end
    print(output_folder)
    print(r1_stripped)
    command = f"/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/SPAdes-3.15.5/bin/spades.py \
        -1 {output_folder}{r1_stripped}_trimmed.fastq.gz \
        -2 {output_folder}{r2_stripped}_trimmed.fastq.gz \
        -t 8 -m 4 --only-assembler -k 21,33,55,77 -o {output_folder}assembly"
    print("Step 2: " + command)
    subprocess.call(command, shell=True)

def step_02_minimap2(file_name, output_folder, trimmed_folder, r1_file, r2_file, ref_seq):
    #Execute:
    #minimap2 -a -t 2 -x sr sequence_reference.fasta assembly/$FILE_R1_trim.fastq.gz assembly/FILE_R2_trim.fastq.gz
    # | samtools view -bS -F 4 - | samtools sort -o LACENMT_20240516_1173_S22.sorted.bam
    r1_stripped = r1_file.strip(".fastq.gz")
    r2_stripped = r2_file.strip(".fastq.gz")

    if not os.path.exists(f"{trimmed_folder}assembly/tmp/"):
        os.makedirs(f"{trimmed_folder}assembly/tmp/")

    command = f"minimap2 -a -t 2 -x sr {ref_seq} {trimmed_folder}{r1_stripped}_trimmed.fastq.gz \
        {trimmed_folder}{r2_stripped}_trimmed.fastq.gz \
            | samtools view -bS -F 4 - | samtools sort -o {output_folder}{r1_stripped}_sorted.bam"

    print("Step 2: " + command)
    
    subprocess.call(command, shell=True)

def step_03(ref_seq, output_folder, file_name):
    # Execute: 
    # minimap2 -x sr --frag=yes --secondary=yes 
    # -N 5 -p 0.8 -a 
    # refseq.fasta assembly/scaffolds.fasta 
    # | samtools view -bS -F 4 - 
    # | samtools sort -o $NAMEFILE"_sorted.bam" -
    command = f"minimap2 -x sr --frag=yes --secondary=yes -N 5 -p 0.8 -a \
        {ref_seq} {output_folder}assembly/scaffolds.fasta \
        | samtools view -bS -F 4 - \
        | samtools sort -o {output_folder}{file_name}_sorted.bam"
    print("Step 3: " + command)
    subprocess.call(command, shell=True)

def step_04(file_name, output_folder):
    # Execute: 
    # samtools index $NAMEFILE"_sorted.bam"
    command = f"samtools index {output_folder}{file_name}_sorted.bam"
    print("Step 4: " + command)
    subprocess.call(command, shell=True)

def step_05(ref_seq, output_folder, file_name):
    # Execute: 
    # SEQ_NAME=`basename refseq.fasta`
    # REF_NAME=`cat refseq.fasta | grep '>' | tr -d '>' | cut -d ' ' -f 1`
    # LENGTH=`tail -n +2 refseq.fasta | tr -d '\n' | wc -m | xargs`
    # echo -e "$REF_NAME\t$LENGTH" > my.genome
    # bedtools bamtobed -i $NAMEFILE"_sorted.bam" > reads.bed
    # bedtools genomecov -bga -i reads.bed -g my.genome | awk '$4 < 1' > zero.bed
    # maskFastaFromBed -fi refseq.fasta -bed zero.bed -fo masked.fasta
    # bcftools mpileup -Ou -f masked.fasta $NAMEFILE"_sorted.bam" | bcftools call --ploidy 1 -mv -Oz -o test.vcf.gz
    # bcftools index test.vcf.gz
    # cat masked.fasta | bcftools consensus test.vcf.gz > new_consensus.fasta
    # echo ">$SEQ_NAME" > $NAMEFILE"_draft.fasta"
    # tail -n +2 new_consensus.fasta >> $NAMEFILE"_draft.fasta"
    # sed 's/>'"$SEQ_NAME"'/>'"$SEQ_NAME"'/g' $NAMEFILE"_draft.fasta" > $NAMEFILE"_consensus.fasta"

    #seq_name = {os.path.basename(ref_seq)}
    command = f"cat {ref_seq} | grep '>' | tr -d '>' | cut -d ' ' -f 1"
    ref_name = subprocess.getoutput(command)
    command = f"tail -n +2 {ref_seq} | tr -d '\\n' | wc -m | xargs"
    length = subprocess.getoutput(command)
    command = f"echo {ref_name}\\\t{length} > {output_folder}{file_name}.genome"
    subprocess.check_call(command, shell=True)
    subprocess.check_call(f"bedtools bamtobed -i {output_folder}{file_name}_sorted.bam > {output_folder}{file_name}_reads.bed", shell=True)
    subprocess.check_call(f"bedtools genomecov -bga -i {output_folder}{file_name}_reads.bed -g {output_folder}{file_name}.genome | awk '$4 < 1' > {output_folder}{file_name}_zero.bed", shell=True)
    subprocess.check_call(f"maskFastaFromBed -fi {ref_seq} -bed {output_folder}{file_name}_zero.bed -fo {output_folder}{file_name}_masked.fasta", shell=True)
    subprocess.check_call(f"bcftools mpileup -Ou -f {output_folder}{file_name}_masked.fasta {output_folder}{file_name}_sorted.bam | bcftools call --ploidy 1 -mv -Oz -o {output_folder}{file_name}_test.vcf.gz", shell=True)
    subprocess.check_call(f"bcftools index {output_folder}{file_name}_test.vcf.gz", shell=True)
    subprocess.check_call(f"cat {output_folder}{file_name}_masked.fasta | bcftools consensus {output_folder}{file_name}_test.vcf.gz > {output_folder}{file_name}_new_consensus.fasta", shell=True)
    subprocess.check_call(f"echo {ref_name} > {output_folder}{file_name}_draft.fasta", shell=True)
    subprocess.check_call(f"tail -n +2 {output_folder}{file_name}_new_consensus.fasta >> {output_folder}{file_name}_draft.fasta", shell=True)
    subprocess.check_call(f"sed 's/>'\"{file_name}\"'/>'\"{file_name}\"'/g' {output_folder}{file_name}_draft.fasta > {output_folder}{file_name}_consensus.fasta", shell=True)
trimmomatic()