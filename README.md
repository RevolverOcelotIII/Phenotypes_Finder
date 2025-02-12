# TCC - VAMPIRE - FERRAMENTA WEB PARA IDENTIFICAÇÃO DE FENÓTIPOS DO GRUPO SANGUÍNEO RH MAPEANDO MUTAÇÕES DOS GENES RHCE E RHD
Trabalho de Conclusão do curso de Sistemas de Informação da Universidade do Estado da Bahia

## Configuração do Ambiente

Este documento fornece as instruções de instalação e configuração para o ambiente de desenvolvimento, separado em **Front-end (React)**, **Back-end (Flask API)** e **Pipeline**.

---

## 1. Requisitos Gerais
Antes de iniciar a instalação, certifique-se de ter os seguintes softwares instalados no sistema:

- **Git**: Controle de versão
- **Python 3.10+**: Necessário para a API e o Worker
- **Node.js 18+** e **npm**: Necessário para o Front-end
- **Java (JDK 11+)**: Necessário para algumas ferramentas do Pipeline
- **MySQL 8.0+**: Banco de dados

Para instalar no Ubuntu/Debian:
```bash
sudo apt update && sudo apt install -y git python3 python3-pip nodejs npm default-jdk mysql-server
```


---

## 2. Instalação do Front-end (React)

### 2.1. Clonar o repositório
```bash
git clone https://github.com/RevolverOcelotIII/GenomeInterface.git
cd GenomeInterface
```

### 2.2. Criar arquivo `.env`
Crie um arquivo `.env` dentro do diretório `GenomeInterface` com o seguinte conteúdo:

```env
REACT_APP_FETCH_API_URL=http://127.0.0.1:5000
```
Onde `REACT_APP_FETCH_API_URL` trata-se do endereço da API, que pode ser configurado posteriormente.


### 2.3. Instalar dependências
```bash
npm install
```

### 2.4. Executar o ambiente de desenvolvimento
```bash
npm start
```

O front-end ficará disponível em `http://localhost:3000/`.

---

## 3. Instalação do Banco de Dados (MySQL)

### 3.1. Iniciar o serviço do MySQL
```bash
sudo systemctl start mysql
```

### 3.2. Acessar o MySQL como root
```bash
mysql -u root -p
```

### 3.3. Criar o banco de dados e configurar usuário
```sql
CREATE DATABASE snpfinder;
CREATE USER 'genome_user'@'localhost' IDENTIFIED BY 'password';
GRANT ALL PRIVILEGES ON snpfinder.* TO 'root'@'localhost';
FLUSH PRIVILEGES;
EXIT;
```

### 3.4. Importar estrutura do banco de dados
```bash
mysql -u root -p snpfinder < Genome-API/arquivo.sql
```

## 4. Instalação do Back-end (Flask API)

### 4.1. Clonar o repositório
```bash
git clone https://github.com/RevolverOcelotIII/Genome-API.git
cd Genome-API
```

### 3.3. Criar arquivo `.env`
Crie um arquivo `.env` dentro do diretório `Genome-API` com o seguinte conteúdo:

```env
ROOT_PATH=/caminho/absoluto/para/genome-api/
DATABASE_URL=127.0.0.1
DATABASE_PORT=3306
DATABASE_USER=root
DATABASE_PASSWORD=password
DATABASE_SCHEMA=snpfinder
```

### 4.4. Instalar dependências Python
```bash
pip install -r requirements.txt
```

### 4.5. Executar a API
```bash
python src/index.py  
```

A API ficará disponível em `http://localhost:5000/`.

---

## 5. Instalação do Pipeline

O Pipeline depende de diversas ferramentas externas para processamento genômico. As instruções abaixo incluem a instalação de cada uma. Para garantir a integração com o API, é necessário clonar o repositório em um diretório dentro da pasta da api. Ele será utilizado pela API, portanto não é necessário executá-lo individualmente.

### 5.1. Entrar no diretório necessário da API

```bash
cd Genome-API
mkdir pipeline
cd pipeline
```

### 5.2. Clonar o repositório
```bash
git clone https://github.com/seu-usuario/genome-pipeline.git
cd genome-pipeline
```

### 5.3. Instalar ferramentas externas
```bash
sudo apt install -y bedtools bcftools minimap2 samtools bwa mafft
```

#### Instalar SPAdes
```bash
wget https://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz
sudo tar -xzf SPAdes-3.15.5-Linux.tar.gz -C /opt/
sudo ln -s /opt/SPAdes-3.15.5-Linux/bin/* /usr/local/bin/
```

#### Instalar Trimmomatic
```bash
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip -d /opt/trimmomatic
```

#### Instalar Pilon
```bash
wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar
sudo mv pilon-1.24.jar /usr/local/bin/pilon.jar
```


