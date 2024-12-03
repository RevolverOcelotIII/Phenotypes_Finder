from config import ROOT_PATH
from dotenv import load_dotenv
import pysam
import vcf
import os
import re
import pandas
from openpyxl import load_workbook
from openpyxl.styles import Font, Border, Side

#from Scrap_Snp_Table import get_gene_mutation_table 

#this should be on other file
import csv

phenotype_tables_path = f"{ROOT_PATH}pipeline/Phenotypes_Finder/PhenotypeTables/"

def get_gene_mutation_table(gene):
    pass
    file_path = f'{phenotype_tables_path}{gene}PhenotypesTable.csv'
    # Ler o arquivo CSV e retornar um array de objetos
    with open(file_path, mode='r', encoding='utf-8-sig') as file:
        reader = csv.DictReader(file, delimiter=',')
        return [row for row in reader]


#this should be on other file

terminal_bg_white = "\033[47m"
terminal_bg_cyan = "\033[46m"
terminal_black = "\033[30m"
terminal_reset = "\033[0m"

def count_phenotype_snps(phenotype_list, mutation, is_Rhd):
    snp_pattern = r'c\.\d+[A-Z]>[A-Z]'
    if mutation['Phenotype'] not in phenotype_list:
        if is_Rhd:
            phenotype_list[mutation['Phenotype']] = [1, len(re.findall(snp_pattern, mutation['Nucleotide change'])), mutation['Classification']]
        else:
            phenotype_list[mutation['Phenotype']] = [1, len(re.findall(snp_pattern, mutation['Nucleotide change']))]
    else:
        phenotype_list[mutation['Phenotype']][0]+=1

def find_snp_type(gene, relative_position, ref, alt, phenotype_list):
    mutation_table = get_gene_mutation_table(gene)
    mutation_types = [mutation for mutation in mutation_table if f'c.{relative_position}{ref}>{alt}' in mutation['Nucleotide change']]
    #print(mutation_type)
    for mutation in mutation_types:
        count_phenotype_snps(phenotype_list, mutation, gene == "RHD")
    return mutation_types

def find_exome_position(file_path, full_position):
    with open(file_path, 'r') as file:
        # Lê todas as linhas do arquivo
        lines = file.readlines()

    # Ignora a primeira linha
    sequence = ''.join(lines[1:]).replace('\n', '')

    # Contabilizar quantos caracteres válidos ocorreram até a posição desejada
    contador_validos = 0
    for indice, caractere in enumerate(sequence):
        if caractere != 'N':
            contador_validos += 1
        if indice + 1 == full_position:
            return contador_validos

    return None  # Caso a posição total esteja fora do intervalo do arquivo

def get_exome_sequence(file_path):
    with open(file_path, 'r') as file:
        # Lê todas as linhas do arquivo
        lines = file.readlines()
    return ''.join(lines[1:]).replace('N', '').replace('\n', '').replace(' ', '')

def find_exon(relative_position, gene):
    if gene == 'RHD':
        if 1 <= relative_position <= 187: return 1
        if 188 <= relative_position <= 374: return 2
        if 375 <= relative_position <= 525: return 3
        if 526 <= relative_position <= 673: return 4
        if 674 <= relative_position <= 840: return 5
        if 841 <= relative_position <= 978: return 6
        if 979 <= relative_position <= 1112: return 7
        if 1113 <= relative_position <= 1192: return 8
        if 1193 <= relative_position <= 1266: return 9
        if 1267 <= relative_position <= 2814: return 10
    elif gene == 'RHCE':
        if 1 <= relative_position <= 187: return 1
        if 188 <= relative_position <= 374: return 2
        if 375 <= relative_position <= 525: return 3
        if 526 <= relative_position <= 673: return 4
        # this exon have one more Nucleotide
        if 674 <= relative_position <= 841: return 5
        # this exon have one less
        if 842 <= relative_position <= 978: return 6
        if 979 <= relative_position <= 1112: return 7
        if 1113 <= relative_position <= 1192: return 8
        if 1193 <= relative_position <= 1266: return 9
        if 1267 <= relative_position <= 2814: return 10
    return None

def format_sheet(csv_path, csv_folder):
    df = pandas.read_csv(csv_path,  sep=",", skipinitialspace=True)
    df.to_excel(f"{csv_folder}/SNP_Sheet.xlsx", index=False)
    workbook = load_workbook(f"{csv_folder}/SNP_Sheet.xlsx")
    sheet = workbook.active
    border = Border(
        left=Side(border_style="thin", color="000000"),
        right=Side(border_style="thin", color="000000"),
        top=Side(border_style="thin", color="000000"),
        bottom=Side(border_style="thin", color="000000")
    )
    for cell in sheet[1]:
        cell.font = Font(bold=True)

def find_snps(project_uuid, sample_name, sample_uuid, gene):
    result_folder = f"{ROOT_PATH}pipeline/Phenotypes_Finder/result/Trimmomatic/{sample_uuid}/{gene}/"
    sheet_folder = f"{ROOT_PATH}pipeline/Phenotypes_Finder/result/SNP_Sheet/{sample_uuid}/{gene}/"
    sheet_path = f'{sheet_folder}SNP_Sheet.csv'
    sample_export = {}
    ignore_sheet = False
    results = os.listdir(result_folder)
    snp_counter = 0
    snp_history = []
    if not os.path.exists(sheet_folder):
        os.makedirs(sheet_folder, exist_ok=True)

    if os.path.exists(sheet_path):
        ignore_sheet = True


    with open(sheet_path, mode='a', newline='', encoding='utf-8') as snp_sheet:
        sheet_columns = 'Sample, RefSeq, Genome Position, Exome Position, Exon, Reference, Mutation, Possible Phenotypes, Quality, Sample Phenotype, Phenotype SNPs, Classification'
        if not ignore_sheet:
            snp_sheet.write(sheet_columns)
        #for result in results:
        # Caminho para o arquivo BAM
        bam_file = [f for f in os.listdir(f"{result_folder}") if f.endswith("_sorted.bam")]

        if len(bam_file) == 0:
            return

        # Caminho para o arquivo VCF gerado no pipeline
        vcf_file = [f for f in os.listdir(f"{result_folder}") if f.endswith("_test.vcf.gz")]

        # Caminho do consenso gerado no pipeline
        consensus_file = [f for f in os.listdir(f"{result_folder}") if f.endswith("_new_consensus.fasta")]
        sample_export['aligned_sequence'] = get_exome_sequence(f"{result_folder}/{consensus_file[0]}")
        sample_export['snp_info'] = []

        # Abrir o arquivo BAM
        bam = pysam.AlignmentFile(f"{result_folder}/{bam_file[0]}", "rb")

        # Inicializar o leitor VCF
        vcf_reader = vcf.Reader(filename=f"{result_folder}/{vcf_file[0]}")
        phenotype_list = {}
        print(f'{terminal_bg_white}{terminal_black}Sequência: {sample_name}{terminal_reset}')
        snp_sample_count = 0
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
                                snp_sample_count += 1
                                #gene = ''
                                #if 'RHD' in result:
                                #    gene = 'RHD'
                                #else:
                                #    gene = 'RHCE'
                                
                                relative_position = find_exome_position(f"{result_folder}/{consensus_file[0]}", pos)
                                #print(relative_position)
                                if relative_position is not None and relative_position!=0: 
                                    snp_types = find_snp_type(gene, relative_position, ref, alt, phenotype_list)
                                    possible_phenotypes = [snp_type['Phenotype'] for snp_type in snp_types]
                                    phenotypes_string = ", ".join(possible_phenotypes).replace("\n", "")
                                    snp_row =  f'\n{sample_name}, { chrom }, { pos }, { relative_position }, { find_exon(relative_position, gene) }, { ref }, { alt }, "{ phenotypes_string }", {qual}'
                                    sample_export['snp_info'].append({
                                        'genome_position': pos,
                                        'exome_position': relative_position,
                                        'exon': find_exon(relative_position, gene),
                                        'reference': ref,
                                        'mutation': alt,
                                        'phenotypes_count': len(possible_phenotypes),
                                        'quality': qual
                                    })
                                    print(f"SNP encontrado em {chrom}:{pos} ( exon N/A:c.{relative_position} ). Ref: {ref}, Alt: {alt} | Deleção? {isdeletion} | Indel? {indel} | Tipos(s): {len([snp_type['Phenotype'] for snp_type in snp_types])} possíveis | Qual: {qual}")
                                else:
                                    snp_row = f'\n{sample_name}, { chrom }, { pos }, , , { ref }, { alt }, , {qual}'
                                    print(f"SNP encontrado em {chrom}:{pos}, Ref: {ref}, Alt: {alt} | Deleção? {isdeletion} | Indel? {indel} | Qual: {qual}")

                                if not ignore_sheet:
                                    snp_sheet.write(snp_row)
                                if pos not in snp_history:
                                    snp_history.append(pos)
                                    snp_counter +=1
                                break

        # Fechar os arquivos
        bam.close()
        #print(f'{terminal_bg_cyan}{terminal_black}Fenótipos encontrados em {result}:{terminal_reset}')
        #if len(phenotype_list.items()) > 0:
        # if not ignore_sheet:
        #   snp_sheet.write('\nSample Phenotypes, Phenotype SNPs')
        sample_export['phenotypes'] = []
        for phenotype, snp_data in phenotype_list.items():
            if snp_data[0]==snp_data[1]:
                no_break_phenotype = phenotype.replace("\n", "")
                if len(snp_data) > 2:
                    sample_export['phenotypes'].append({
                        'name': no_break_phenotype,
                        'snp_count': snp_data[1],
                        'classification': snp_data[2]
                    })
                    if not ignore_sheet:
                        snp_sheet.write(f'\n{sample_name}, , , , , , , , , "{no_break_phenotype}", {snp_data[1]}, {snp_data[2]}')
                else:
                    sample_export['phenotypes'].append({
                        'name': no_break_phenotype,
                        'snp_count': snp_data[1],
                        'classification': '-'
                    })
                    if not ignore_sheet:
                        snp_sheet.write(f'\n{sample_name}, , , , , , , , , "{no_break_phenotype}", {snp_data[1]}')
                print(f"{phenotype}: Foram encontradas {snp_data[0]} de {snp_data[1]} SNPs")
        if snp_sample_count == 0:
            if gene == 'RHD':
                sample_export['phenotypes'].append({
                    'name': '-',
                    'snp_count': 0,
                    'classification': 'D'
                })
                if not ignore_sheet:
                    snp_sheet.write(f'\n{sample_name}, , , , , , , , , -, No SNPs found (same as reference), D')
            else:
                sample_export['phenotype'].append({
                    'name': '-',
                    'snp_count': 0,
                    'classification': '-'
                })
                if not ignore_sheet:
                    snp_sheet.write(f'\n{sample_name}, , , , , , , , , -, No SNPs found (same as reference), -')

        print("SNPs encontrados: ", len(snp_history))
        # print(snp_history)
    if snp_sample_count!= 0:
        format_sheet(sheet_path, sheet_folder)
    print(sample_export)
    return sample_export

#find_snps()
