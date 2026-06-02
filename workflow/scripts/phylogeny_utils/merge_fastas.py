import os
import glob
from Bio import SeqIO

def gerar_super_fasta(diretorio_entrada: str, arquivo_saida: str):
    """
    Lê todos os arquivos '*_consensus.fasta' de um diretório, 
    substitui o cabeçalho interno pelo nome do arquivo (ex: >LBCM_LS038)
    e condensa tudo num único 'Super FASTA'.
    """
    # 1. Pega todos os arquivos que terminam com _consensus.fasta na pasta
    padrao_busca = os.path.join(diretorio_entrada, "*_consensus.fasta")
    arquivos_fasta = glob.glob(padrao_busca)
    
    if not arquivos_fasta:
        print(f"⚠️ Nenhum arquivo '_consensus.fasta' encontrado na pasta: {diretorio_entrada}")
        return

    sequencias_limpas = []
    
    print(f"🔄 Processando {len(arquivos_fasta)} arquivos...")

    # 2. Iterar sobre cada arquivo encontrado
    for caminho_arquivo in arquivos_fasta:
        # Extrai apenas o nome do arquivo, sem o caminho da pasta
        nome_arquivo = os.path.basename(caminho_arquivo)
        
        # Remove a parte final para pegar só a tag (ex: LS038)
        tag_amostra = nome_arquivo.replace("_consensus.fasta", "")
        
        # Monta a nova etiqueta oficial
        novo_cabecalho = f"LBCM_{tag_amostra}"
        
        try:
            # Lê a sequência de dentro do arquivo
            # Como é um arquivo de consenso, esperamos apenas 1 sequência por arquivo
            record = SeqIO.read(caminho_arquivo, "fasta")
            
            # Substitui o ID e limpa qualquer descrição extra que estivesse lá
            record.id = novo_cabecalho
            record.description = ""
            record.name = "" # Limpa para evitar redundâncias do Biopython
            
            # Guarda na lista
            sequencias_limpas.append(record)
            
        except ValueError:
            print(f"Erro ao ler {nome_arquivo}: O arquivo está vazio ou tem mais de uma sequência.")
        except Exception as e:
            print(f"Erro inesperado em {nome_arquivo}: {e}")

    # 3. Salva todas as sequências empilhadas no novo arquivo Super FASTA
    if sequencias_limpas:
        SeqIO.write(sequencias_limpas, arquivo_saida, "fasta")
        print("\n" + "="*60)
        print(f"Sucesso! Super FASTA criado com {len(sequencias_limpas)} sequências.")
        print(f"Arquivo salvo em: {arquivo_saida}")
        print("="*60 + "\n")

# ==========================================
# COMO EXECUTAR NO COLAB
# ==========================================
# Mude "./pasta_dos_consensos" para o caminho real da sua pasta 
pasta_origem = "/home/paulo/1_LAB/YeastProspecta/results/02_consensus/fasta" 
arquivo_final = "/home/paulo/1_LAB/YeastProspecta/results/Banco_LBCM_Completo.fasta"

gerar_super_fasta(pasta_origem, arquivo_final)