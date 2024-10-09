import os
import re
import requests
from bs4 import BeautifulSoup

# Diretório onde os arquivos estão localizados
directory = '/home/domdeny/src/bioinfo/pipeline-jessica/PipelineJessica/dataset/FASTQ'

# Regex para extrair a informação desejada
regex_pattern = r'[A-Za-z0-9\-]+_S\d+_L\d{3}_R[12]_001'

# Função para buscar o conteúdo no site e extrair a string com o regex
def get_renaming_pattern(srr_id):
    url = f"https://www.ncbi.nlm.nih.gov/sra/{srr_id}"
    response = requests.get(url)
    if response.status_code == 200:
        match = re.search(regex_pattern, response.text)
        if match:
            return match.group(0)  # Retorna a string encontrada
    return None

# Função principal para processar os arquivos
def process_fastq_files(directory):
    for filename in os.listdir(directory):
        if filename.startswith("SRR") and filename.endswith("_1.fastq.gz"):
            # Extrai o ID SRR
            srr_id = filename.split('_')[0]
            
            # Busca o padrão de renomeação no site NCBI
            new_pattern = get_renaming_pattern(srr_id)
            print(new_pattern)

            reverse_sra=[f for f in os.listdir(directory) if f.startswith(f"{srr_id}_2")]
            
            if new_pattern:
                # Renomear os arquivos R1 e R2
                print(filename)
                new_filename = f"{new_pattern}.fastq.gz"
                new_reverse_filename = f"{new_pattern.replace('R1', 'R2')}.fastq.gz"
                
                # Renomeia o arquivo
                old_file_path = os.path.join(directory, filename)
                new_file_path = os.path.join(directory, new_filename)
                
                print(f"Renomeando {old_file_path} para {new_file_path}")
                os.rename(old_file_path, new_file_path)

                old_file_path = os.path.join(directory, reverse_sra[0])
                new_file_path = os.path.join(directory, new_reverse_filename)
                
                print(f"Renomeando {old_file_path} para {new_file_path}")
                os.rename(old_file_path, new_file_path)
            else:
                print(f"Padrão não encontrado para {srr_id}")

# Executar o script no diretório especificado
process_fastq_files(directory)
