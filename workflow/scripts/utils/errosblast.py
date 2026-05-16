import os
import xml.etree.ElementTree as ET

xml_dir = "results/03_blast/raw_hits"

print(f"{'Amostra':<20} | {'Hits Encontrados':<15}")
print("-" * 40)

# Lista os arquivos em ordem alfabética para facilitar a leitura
files = sorted([f for f in os.listdir(xml_dir) if f.endswith(".xml")])

total_hits = 0
empty_files = 0

for filename in files:
    path = os.path.join(xml_dir, filename)
    sample_id = filename.replace("_blast.xml", "")
    
    try:
        # Carrega o XML
        tree = ET.parse(path)
        root = tree.getroot()
        
        # O BLAST XML organiza hits dentro de Iteration -> Iteration_hits -> Hit
        # .//Hit encontra todos os elementos Hit em qualquer nível
        hits = root.findall(".//Hit")
        count = len(hits)
        
        print(f"{sample_id:<20} | {count:<15}")
        
        total_hits += count
        if count == 0:
            empty_files += 1
            
    except Exception as e:
        print(f"{sample_id:<20} | ❌ ERRO AO LER: {e}")

print("-" * 40)
print(f"Resumo: {len(files)} arquivos analisados.")
print(f"Total de hits processados: {total_hits}")
print(f"Amostras sem nenhum hit: {empty_files}")
