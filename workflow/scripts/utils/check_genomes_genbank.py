import csv
import time
from Bio import Entrez

# --- 1. Configurações da API ---
Entrez.email = "seu_email@exemplo.com"  # OBRIGATÓRIO: altere para o seu e-mail
Entrez.api_key = "7c89917b1583ac920b89d9ab17f7a98e3608"     # Cole sua API key para não ser bloqueado

# --- 2. O Dicionário de Espécies ---
# 'Sungouiella pseudointermedia' e 'Wickerhamomyces sydowiorum' apareceram nas duas listas.
# Assumi como "Exata" (prioridade) para evitar duplicatas, mas você pode alterar.
dicionario_especies = {
    # --- PRÓXIMAS ---
    "Candida edax": "Próxima",
    "Debaryomyces vietnamensis": "Próxima",
    "Maudiozyma exigua": "Próxima",
    "Maudiozyma serrabonitensis": "Próxima",
    "Meyerozyma caribbica": "Próxima",
    "Ogataea salicorniae": "Próxima",
    "Ogataea trehalophila": "Próxima",
    "Ogataea degrootiae": "Próxima",
    "Scheffersomyces lignosus": "Próxima",
    "Scheffersomyces shehatae": "Próxima",
    "Scheffersomyces parashehatae": "Próxima",
    "Scheffersomyces quercinus": "Próxima",
    "Scheffersomyces insectosa": "Próxima",
    "Scheffersomyces tanahashii": "Próxima",
    "Scheffersomyces stipitis": "Próxima",
    "Scheffersomyces segobiensis": "Próxima",
    "Scheffersomyces titani": "Próxima",
    "Sugiyamaella floridensis": "Próxima",
    "Vanrija meifongana": "Próxima",
    
    # --- EXATAS ---
    "Barnettozyma californica": "Exata",
    "Blastobotrys peoriensis": "Exata",
    "Candida albicans": "Exata",
    "Candida awuae": "Exata",
    "Candida berthetii": "Exata",
    "Candida boidinii": "Exata",
    "Candida kantuleensis": "Exata",
    "Candida michaelii": "Exata",
    "Candida sinolaborantium": "Exata",
    "Candida sojae": "Exata",
    "Candida solani": "Exata",
    "Candida vrieseae": "Exata",
    "Colacogloea diffluens": "Exata",
    "Cyberlindnera mrakii": "Exata",
    "Cyberlindnera subsufficiens": "Exata",
    "Debaryomyces nepalensis": "Exata",
    "Jamesozyma lodderae": "Exata",
    "Jamesozyma spencerorum": "Exata",
    "Maudiozyma barnettii": "Exata",
    "Maudiozyma saulgeensis": "Exata",
    "Meyerozyma guilliermondii": "Exata",
    "Monosporozyma solicola": "Exata",
    "Ogataea cecidiorum": "Exata",
    "Ogataea ganodermae": "Exata",
    "Papiliotrema aurea": "Exata",
    "Papiliotrema laurentii": "Exata",
    "Papiliotrema pseudoalba": "Exata",
    "Papiliotrema tapputi": "Exata",
    "Pichia dushanensis": "Exata",
    "Pichia kudriavzevii": "Exata",
    "Pichia pseudolambica": "Exata",
    "Pseudotremella allantoinivorans": "Exata",
    "Pseudotremella lacticolor": "Exata",
    "Saturnispora silvae": "Exata",
    "Scheffersomyces illinoinensis": "Exata",
    "Scheffersomyces xylosifermentans": "Exata",
    "Sugiyamaella boreocaroliniensis": "Exata",
    "Sugiyamaella smithiae": "Exata",
    "Suhomyces ambrosiae": "Exata",
    "Suhomyces chickasaworum": "Exata",
    "Sungouiella blattae": "Exata",
    "Sungouiella intermedia": "Exata",
    "Sungouiella pseudointermedia": "Exata", # Removida da lista de "Próximas" para evitar duplicação
    "Tarzetta alpina": "Exata",
    "Vanrija humicola": "Exata",
    "Wickerhamomyces arborarius": "Exata",
    "Wickerhamomyces bovis": "Exata",
    "Wickerhamomyces sydowiorum": "Exata", # Removida da lista de "Próximas" para evitar duplicação
    "Yamadazyma dushanensis": "Exata",
    "Yamadazyma luoyangensis": "Exata",
    "Yamadazyma riverae": "Exata"
}

# --- 3. Função de Busca e Extração ---
def buscar_dados_genoma(especie):
    # Procura por qualquer genoma (removido o filtro "latest refseq" para pegar mais resultados)
    termo_busca = f'"{especie}"[Organism]'
    
    try:
        # Busca inicial (quantos existem)
        handle_busca = Entrez.esearch(db="assembly", term=termo_busca, retmode="xml")
        record_busca = Entrez.read(handle_busca)
        handle_busca.close()
        
        ids_encontrados = record_busca["IdList"]
        quantidade = len(ids_encontrados)
        
        if quantidade == 0:
            return "N/A", 0, "N/A", "N/A", "N/A", "N/A"
            
        # Pega os detalhes do primeiro genoma retornado (geralmente o mais relevante/recente)
        top_id = ids_encontrados[0]
        
        handle_sumario = Entrez.esummary(db="assembly", id=top_id, retmode="xml")
        record_sumario = Entrez.read(handle_sumario)
        handle_sumario.close()
        
        # Extração de métricas de qualidade (o XML do NCBI é aninhado)
        document_summary = record_sumario['DocumentSummarySet']['DocumentSummary'][0]
        
        accession = document_summary.get('AssemblyAccession', 'N/A')
        assembly_level = document_summary.get('AssemblyStatus', 'N/A')
        coverage = document_summary.get('Coverage', 'N/A')
        tamanho_mb = "N/A"
        biosample = document_summary.get('BioSampleAccn', 'N/A')
        
        # O tamanho às vezes vem aninhado na metadados
        meta = document_summary.get('Meta', '')
        if 'Stat type="total_length"' in meta:
            # Extração simples via split de string para evitar parseadores de XML pesados
            try:
                parte1 = meta.split('Stat type="total_length">')[1]
                tamanho_bases = int(parte1.split('</Stat>')[0])
                tamanho_mb = f"{(tamanho_bases / 1000000):.2f} Mb"
            except:
                pass
        
        return accession, quantidade, assembly_level, coverage, tamanho_mb, biosample
        
    except Exception as e:
        print(f"Erro ao buscar {especie}: {e}")
        return "Erro", "Erro", "Erro", "Erro", "Erro", "Erro"

# --- 4. Execução e Geração do CSV ---
nome_arquivo = "relatorio_genomas_ncbi.csv"

print("Iniciando busca no GenBank. Isso pode levar alguns minutos...")

with open(nome_arquivo, mode='w', newline='', encoding='utf-8') as arquivo_csv:
    writer = csv.writer(arquivo_csv, delimiter=';') # Ponto e vírgula é melhor para abrir no Excel BR
    
    # Cabeçalho atualizado com as novas métricas
    writer.writerow([
        "Espécie", 
        "Categoria", 
        "Top Accession", 
        "Quantidade Disponível", 
        "Nível de Montagem", 
        "Cobertura", 
        "Tamanho (Mb)", 
        "BioSample ID"
    ])
    
    for especie, categoria in dicionario_especies.items():
        print(f"Buscando: {especie}...")
        
        # Executa a função
        acc, qtd, status_assembly, cov, tam, b_sample = buscar_dados_genoma(especie)
        
        # Escreve a linha no CSV
        writer.writerow([especie, categoria, acc, qtd, status_assembly, cov, tam, b_sample])
        
        # Pausa obrigatória para não sobrecarregar a API do NCBI (3 a 10 req/s)
        time.sleep(0.35) 

print(f"\nConcluído! Relatório salvo como '{nome_arquivo}'.")
