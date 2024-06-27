import pysam
import vcf
import os
from Scrap_Snp_Table import get_gene_mutation_table 

def find_snp_type(gene, relative_position, ref, alt):
    mutation_table = get_gene_mutation_table(gene)
    mutation_type = [mutation for mutation in mutation_table if f'c.{relative_position}{ref}>{alt}' in mutation['Nucleotide change']]
    return mutation_type


def find_exon_position(position, gene):
    if gene == 'RHD':
        if position >= 5020 and position < 5206:
            return [position - 5020 + 1, 1]
        if position >= 17084 and position <= 17270:
            return [position - 17084 + 186, 2]
        if position >= 23152 and position <= 23302:
            return [position - 23152 + 336, 3]
        if position >= 33457 and position <= 33604:
            return [position - 33457 + 483, 4]
        if position >= 34031 and position <= 34197:
            return [position - 33334 + 649, 5]
        if position >= 35833 and position <= 35970:
            return [position - 35833 + 786, 6]
        if position >= 39107 and position <= 39240:
            return [position - 39107 + 919, 6]
        if position >= 49511 and position <= 49590:
            return [position - 49511 + 998, 8]
        if position >= 54400 and position <= 54473:
            return [position - 54400 + 1071, 9]
        if position >= 61409 and position <= 62956:
            return [position - 61409 + 2618, 10]
    if gene == 'RHCE':
        if position >= 14368 and position <= 14554:
            return [position - 14368 + 1, 1]
        if position >= 26324 and position <= 26510:
            return [position - 26324 + 186, 2]
        if position >= 32447 and position <= 32597:
            return [position - 32447 + 373, 3]
        if position >= 43052 and position <= 43199:
            return [position - 43052 + 524, 4]
        if position >= 44278 and position <= 44444:
            return [position - 44278 + 671, 5]
        if position >= 46080 and position <= 46217:
            return [position - 46080 + 837, 6]
        if position >= 49349 and position <= 49482:
            return [position - 49349 + 974, 7]
        if position >= 59765 and position <= 59844:
            return [position - 59765 + 1107, 8]
        if position >= 64653 and position <= 64726:
            return [position - 64653 + 1186, 9]
        if position >= 72640 and position <= 72944:
            return [position - 72640 + 1260, 10]
    else:
        return [0,0]

def find_snps():
    result_folder = "/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/result/Trimmomatic/"
    results = os.listdir(result_folder)
    snp_counter = 0
    snp_history = []
    for result in results:
        # Caminho para o arquivo BAM
        bam_file = [f for f in os.listdir(f"{result_folder}{result}") if f.endswith("_sorted.bam")]

        # Caminho para o arquivo VCF gerado no pipeline
        vcf_file = [f for f in os.listdir(f"{result_folder}{result}") if f.endswith("_test.vcf.gz")]

        # Abrir o arquivo BAM
        bam = pysam.AlignmentFile(f"{result_folder}{result}/{bam_file[0]}", "rb")

        # Inicializar o leitor VCF
        vcf_reader = vcf.Reader(filename=f"{result_folder}{result}/{vcf_file[0]}")
        print(f"Sequência: {result}")
        # Loop através das entradas do VCF
        for record in vcf_reader:
            # Obter informações sobre a posição do SNP
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = str(record.ALT[0])  # Assumindo um único SNP
            isdeletion = record.is_deletion
            indel = record.is_indel
            qual = record.QUAL
            
            # Coletar informações do BAM na posição do SNP
            for pileupcolumn in bam.pileup(chrom, pos - 1, pos):
                if pileupcolumn.pos == pos - 1:
                    # Analisar as leituras que cobrem o SNP
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.query_position is not None:
                            # Obter a base na posição do SNP
                            base = pileupread.alignment.query_sequence[pileupread.query_position]

                            # Verificar se a base difere da referência (SNP)
                            if base != ref:
                                gene = ''
                                if 'RHD' in result:
                                    gene = 'RHD'
                                else:
                                    gene = 'RHCe'
                                
                                relative_position = find_exon_position(pos, gene)
                                #print(relative_position)
                                if relative_position is not None: 
                                    if relative_position[0]!=0:
                                        snp_types = find_snp_type(gene, relative_position[0], ref, alt)
                                        print(f"SNP encontrado em {chrom}:{pos} ( exon {relative_position[1]}:c.{relative_position[0]} ). Ref: {ref}, Alt: {alt} | Deleção? {isdeletion} | Indel? {indel} | type(s): {[snp_type['Phenotype'] for snp_type in snp_types]} | Qual: {qual}")
                                    else:
                                        print(f"SNP encontrado em {chrom}:{pos}, Ref: {ref}, Alt: {alt} | Deleção? {isdeletion} | Indel? {indel} | Qual: {qual}")
                                else:
                                    print(f"SNP encontrado em {chrom}:{pos}, Ref: {ref}, Alt: {alt} | Deleção? {isdeletion} | Indel? {indel} | Qual: {qual}")
                                if pos not in snp_history:
                                    snp_history.append(pos)
                                    snp_counter +=1
                                break

        # Fechar os arquivos
        bam.close()
    print("SNPs encontrados: ", len(snp_history))
   # print(snp_history)
find_snps()
#print(find_snp_type(361, 'RHD'))
#print(find_exon_position(61409, 'RHD'))