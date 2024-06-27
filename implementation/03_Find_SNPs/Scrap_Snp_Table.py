import csv

def get_gene_mutation_table(gene):
    pass
    file_path = f'{gene}PhenotypesTable.csv'
    # Ler o arquivo CSV e retornar um array de objetos
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter=',')
        return [row for row in reader]
