from Bio import Entrez, SeqIO
import time
import os
from dotenv import load_dotenv

def buscar_sequencia_ncbi(accession: str, email: str, api_key: str = None, db: str = "nucleotide") -> str:
    """
    Busca uma sequência no NCBI a partir do seu Accession Code usando API Key.
    """
    # 1. Configurações de identificação do NCBI
    Entrez.email = email
    
    if api_key:
        Entrez.api_key = api_key  
    
    try:
        with Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text") as handle:
            registro = SeqIO.read(handle, "fasta")
            return str(registro.seq)
            
    except Exception as e:
        print(f"Erro ao buscar o accession {accession}: {e}")
        return None


load_dotenv("../../../keys/ncbi.env")

meu_email = os.getenv("EMAIL")
minha_chave = os.getenv("NCBI_API_KEY")

lista_accessions = ["NG_055130"]
dicionario_resultados = {}

#for acc in lista_accessions:
    #seq = buscar_sequencia_ncbi(accession=acc, email=meu_email, api_key=minha_chave)
    
    #if seq:
        #dicionario_resultados[acc] = seq
        
    # Uma pequena pausa é sempre uma boa prática, mesmo com a API Key
    #time.sleep(0.15) 

print(meu_email, minha_chave)

print(dicionario_resultados)