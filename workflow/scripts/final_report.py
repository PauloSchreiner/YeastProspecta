import pandas as pd
import numpy as np
from Bio import SearchIO

import re

import os 

from openpyxl.formatting.rule import CellIsRule
from openpyxl.styles import PatternFill



#
# EXTRACTION
#

def parse_blast_xmls(xml_path_list: list) -> pd.DataFrame:
    """
    Lê os arquivos XML do BLAST e retorna um DataFrame 
    com todos os hits de todas as amostras.
    """
    rows = []

    for xml_file in xml_path_list:
        sample = os.path.basename(xml_file).rsplit("_", 1)[0] 
        qresult = SearchIO.read(xml_file, "blast-xml")

        for hit_index, hit in enumerate(qresult, start=1):
            hsp = hit[0] 

            rows.append({
                "Sample": sample.zfill(3),
                "Hit #": hit_index,
                "Accession": hit.accession,
                "Identity": round((hsp.ident_num / hsp.aln_span) * 100, 3),
                "Coverage": round((hsp.query_span / qresult.seq_len) * 100, 3),
                "E-value": hsp.evalue,
                "Title": hit.description
            })

    return pd.DataFrame(rows)



def extract_trim_data(trim_path: str) -> pd.DataFrame:
    """
    Lê o CSV de trimming, passa as reads F e R para colunas
    e retorna um DataFrame indexado pela amostra.
    """
    # 1. Lê o arquivo de uma vez (sep='\t' indica que é separado por Tab)
    df = pd.read_csv(trim_path, sep="\t")
    
    # 2. Divide a coluna 'sample' (ex: "1_F") em duas novas colunas
    # O expand=True faz com que o split crie colunas reais no DataFrame.
    df[['Sample_ID', 'Direction']] = df['sample'].str.rsplit('_', n=1, expand=True)
    
    # 3. Aplica o zfill(3) para padronizar o nome (ex: "1" vira "001")
    df['Sample_ID'] = df['Sample_ID'].astype(str).str.zfill(3)
    
    # 4. Transforma as linhas F e R em colunas
    df_pivot = df.pivot(index='Sample_ID', columns='Direction', values='post_trim_length')
    
    # 5. Renomeia as colunas para o padrão desejado
    df_pivot.rename(columns={'F': 'F_len', 'R': 'R_len'}, inplace=True)
    
    # Limpa o nome do eixo das colunas (estética) e retorna
    df_pivot.columns.name = None 
    return df_pivot


def extract_consensus_data(consensus_path: str) -> pd.DataFrame:
    """
    Lê o CSV de consenso, padroniza as amostras e trata scores vazios/inválidos.
    """
    df = pd.read_csv(consensus_path, sep="\t")
    
    # 1. Renomeia as colunas para o nosso padrão
    df.rename(columns={
        'Sample': 'Sample_ID', 
        'Length': 'Consensus_len', 
        'Score': 'Consensus_score'
    }, inplace=True)
    
    # 2. Padroniza o nome da amostra
    df['Sample_ID'] = df['Sample_ID'].astype(str).str.zfill(3)
    
    # 3. Trata a coluna Score (o equivalente ao try/except do seu código antigo)
    df['Consensus_score'] = pd.to_numeric(df['Consensus_score'], errors='coerce').fillna(0).astype(int)
    
    # 4. Define a Amostra como o Índice (linha) da tabela
    df.set_index('Sample_ID', inplace=True)
    
    return df



#
# PROCESSING
#


def extract_species_from_title(ncbi_title:str) -> str:
    """
    Extracts the binomial species name from a raw NCBI title.
    
    Args: ncbi_title: The raw title from BLAST.
    Returns: The clean 'Genus species' string.
    """

    pattern = r"^(\[?[A-Z][a-z]+\]?\s+[a-z]+\.?)"
    match = re.search(pattern, ncbi_title)

    if match:
        return match.group(1)
    else: 
        return ncbi_title


def is_type(ncbi_title:str) -> bool:
    """
    Given an NCBI title, checks if it is the type strain
    """
    return ("type" in ncbi_title.lower())



def build_qc_df(df_trim: pd.DataFrame, df_cons: pd.DataFrame, df_blast: pd.DataFrame) -> pd.DataFrame:
    """
    Une os dados de Trimming, Consenso e calcula os hits >99% a partir do BLAST.
    Retorna o DataFrame consolidado de Quality Control pronto para exportação.
    """
    # 1. Une Trim e Consenso (Outer join garante que não perdemos amostras isoladas)
    df_qc = df_trim.join(df_cons, how='outer')
    
    # 2. Calcula quantos hits têm mais de 99% de identidade (Vetorizado e rápido)
    hits_over_99 = df_blast[df_blast['Identity'] > 99.0].groupby('Sample').size()
    
    # Dá um nome para essa nova série para facilitar o join
    hits_over_99.name = 'Hits_>99_id_pct' 
    
    # 3. Adiciona a contagem ao DataFrame de QC
    df_qc = df_qc.join(hits_over_99, how='outer')
    
    # 4. Tratamento de NaNs e Tipagem de Dados (Evitando a "maldição do Float")
    # O tipo 'Int64' (com I maiúsculo) do Pandas aceita valores nulos (NaN) sem virar decimal.
    colunas_numericas = ['F_len', 'R_len', 'Consensus_len', 'Consensus_score', 'Hits_>99_id_pct']
    
    for col in colunas_numericas:
        if col in df_qc.columns:
            # Preenche NaNs com 0 e força o tipo para inteiro seguro
            df_qc[col] = df_qc[col].fillna(0).astype('Int64')
            
    # Remove a nomenclatura extra do índice e organiza por ordem alfabética
    df_qc.index.name = "Sample"
    df_qc.sort_index(inplace=True)
            
    return df_qc



def build_summary_df(df_qc: pd.DataFrame, df_blast: pd.DataFrame) -> pd.DataFrame:
    """
    Gera o relatório resumo cruzando dados de QC e os Top Hits do BLAST.
    Totalmente vetorizado para máxima performance.
    """
    
    # 1. Criação das Máscaras Booleanas
    is_type_mask = df_blast["Title"].apply(is_type)
    is_above_99_id_mask = df_blast["Identity"] > 98.90
    
    # Aplica ambas as condições de uma vez só usando o '&'
    df_type = df_blast[is_type_mask & is_above_99_id_mask].copy()


    df_type['Species'] = df_type['Title'].apply(extract_species_from_title)
    df_type['Species'] = df_type['Species'].str.replace(r'[\[\]]', '', regex=True)
    
    # 2. Ordena por Amostra e Identidade (Garante que o melhor hit esteja no topo)
    df_type.sort_values(by=['Sample', 'Identity'], ascending=[True, False], inplace=True)
    
    # Ao invés de iterar, mandamos o Pandas apagar as duplicatas, mantendo só a PRIMEIRA linha de cada amostra
    df_top = df_type.drop_duplicates(subset=['Sample'], keep='first').copy()
    df_top = df_top[['Sample', 'Species', 'Identity']].rename(columns={'Species': 'Top hit'})
    
    # Mantém apenas o melhor hit de CADA espécie por amostra
    df_unique_species = df_type.drop_duplicates(subset=['Sample', 'Species'], keep='first').copy()
    
    # Cria a string formatada "Espécie (99.0%)" para todo mundo de uma vez
    df_unique_species['Formatted'] = df_unique_species['Species'] + " (" + df_unique_species['Identity'].astype(str) + "%)"
    
    # Agrupa por amostra e junta as strings separadas por "; ", pulando a primeira linha (que é o Top Hit)
    df_other = df_unique_species.groupby('Sample').apply(
        lambda x: "; ".join(x['Formatted'].iloc[1:])
    ).reset_index(name='Other possible species')
    
    # --- UNINDO TUDO ---
    # Pega os dados básicos do QC
    df_final = df_qc.reset_index()[['Sample', 'Consensus_score', 'Hits_>99_id_pct']].copy()
    
    # Mescla (Left Join) os Top Hits e as Outras Espécies
    df_final = df_final.merge(df_top, on='Sample', how='left')
    df_final = df_final.merge(df_other, on='Sample', how='left')
    
    # Preenche vazios (Para amostras que não tiveram NENHUM match com Type Strain)
    df_final['Top hit'] = df_final['Top hit'].fillna("None")
    df_final['Identity'] = df_final['Identity'].fillna(0)
    df_final['Other possible species'] = df_final['Other possible species'].fillna("")
    
    # --- REGRAS DE NEGÓCIO (O STATUS) ---
    # np.select avalia múltiplas condições de uma vez só!
    condicoes = [
        df_final['Consensus_score'] < 2000,                                 # Regra 1
        df_final['Other possible species'] != "",                           # Regra 2
        (df_final['Hits_>99_id_pct'] == 0) & (df_final['Identity'] < 99)    # Regra 3
    ]
    
    resultados = [
        "Bad quality",
        "Ambiguous",
        "No match"
    ]
    
    # Se bater na Regra 1, aplica "Bad quality". Se não bater na 1 mas bater na 2, "Ambiguous", etc.
    df_final['Status'] = np.select(condicoes, resultados, default="OK")
    
    # Organiza a ordem das colunas e retorna
    colunas_finais = ["Sample", "Status", "Top hit", "Identity", "Other possible species"]
    
    return df_final[colunas_finais]



#
# PRESENTING
#


def auto_adjust_columns(writer: pd.ExcelWriter, sheet_name: str):
    """
    Adjusts the width of all columns in the specified Excel sheet 
    to fit the content.
    """
    worksheet = writer.sheets[sheet_name] 
    for column in worksheet.columns:
        max_length = 0
        column_letter = column[0].column_letter # Get 'A', 'B', etc.
        
        for cell in column:
            try:
                if len(str(cell.value)) > max_length:
                    max_length = len(str(cell.value))
            except:
                pass
        
        # Capped at 60 characters so really long descriptions don't break the UI
        adjusted_width = min(max_length + 2, 60)
        worksheet.column_dimensions[column_letter].width = adjusted_width


def apply_workbook_styling(writer: pd.ExcelWriter, df_summary: pd.DataFrame, df_blast: pd.DataFrame, df_qc: pd.DataFrame):
    """
    Applies all conditional formatting and styling to the entire Excel workbook.
    """
    # 1. DEFINE PRETTY HUES
    colors = {
        "green":      PatternFill(start_color='C6EFCE', end_color='C6EFCE', fill_type='solid'),
        "lime":       PatternFill(start_color='E3EDB5', end_color='E3EDB5', fill_type='solid'),
        "yellow":     PatternFill(start_color='FFEB9C', end_color='FFEB9C', fill_type='solid'),
        "red":        PatternFill(start_color='FFC7CE', end_color='FFC7CE', fill_type='solid'),
        "purple":     PatternFill(start_color='E4CCFF', end_color='E4CCFF', fill_type='solid'),
        "light_gray": PatternFill(start_color='F2F2F2', end_color='F2F2F2', fill_type='solid'),
        "white":      PatternFill(start_color='FFFFFF', end_color='FFFFFF', fill_type='solid')
    }

    # --- SHEET 1: SUMMARY (Status Colors) ---
    ws_sum = writer.sheets["Summary"]
    status_range = f"B2:B{len(df_summary) + 1}"
    
    ws_sum.conditional_formatting.add(status_range, CellIsRule(operator='equal', formula=['"OK"'], fill=colors["green"]))
    ws_sum.conditional_formatting.add(status_range, CellIsRule(operator='equal', formula=['"No match"'], fill=colors["purple"]))
    ws_sum.conditional_formatting.add(status_range, CellIsRule(operator='equal', formula=['"Ambiguous"'], fill=colors["yellow"]))
    ws_sum.conditional_formatting.add(status_range, CellIsRule(operator='equal', formula=['"Bad quality"'], fill=colors["red"]))

    # --- SHEET 2: BLAST DETAILS (Alternating Sample Groups) ---
    ws_blast = writer.sheets["BLAST_Details"]
    current_fill = colors["white"]
    last_sample = None
    
    for i, sample_id in enumerate(df_blast["Sample"], start=2):
        if last_sample is not None and sample_id != last_sample:
            current_fill = colors["light_gray"] if current_fill == colors["white"] else colors["white"]
        
        for col in range(1, ws_blast.max_column + 1):
            ws_blast.cell(row=i, column=col).fill = current_fill
        last_sample = sample_id

    # --- SHEET 3: QUALITY CONTROL (Consensus Score Rules) ---
    ws_qc = writer.sheets["Quality_Control"]
    
    # Identify the 'Consensus_score' column dynamically
    score_col_letter = None
    for cell in ws_qc[1]:
        if cell.value == "Consensus_score":
            score_col_letter = cell.column_letter
            break
    
    if score_col_letter:
        qc_range = f"{score_col_letter}2:{score_col_letter}{len(df_qc) + 1}"
        ws_qc.conditional_formatting.add(qc_range, CellIsRule(operator='greaterThan', formula=['2250'], fill=colors["green"]))
        ws_qc.conditional_formatting.add(qc_range, CellIsRule(operator='greaterThan', formula=['2000'], fill=colors["lime"]))
        ws_qc.conditional_formatting.add(qc_range, CellIsRule(operator='greaterThan', formula=['1750'], fill=colors["yellow"]))
        ws_qc.conditional_formatting.add(qc_range, CellIsRule(operator='lessThanOrEqual', formula=['1750'], fill=colors["red"]))


def export_to_excel(df_summary: pd.DataFrame, df_blast: pd.DataFrame, df_qc: pd.DataFrame, output_file: str):
    """
    Orquestra a exportação dos DataFrames para as abas do Excel,
    aplica os filtros visuais (Top 20 hits) e aciona as formatações.
    """
    # 1. Aplica o filtro de apresentação: Apenas Top 20 hits na aba BLAST
    df_blast_export = df_blast[df_blast['Hit #'] <= 20].copy()

    # 2. Ordenação final estética
    df_summary.sort_values(by="Sample", inplace=True)
    df_blast_export.sort_values(by=["Sample", "Hit #"], inplace=True)
    
    # Prepara o df_qc transformando o índice 'Sample' em uma coluna real para o Excel
    df_qc_export = df_qc.reset_index() 

    # 3. Escrita no arquivo
    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
        
        # Escreve os DataFrames
        df_summary.to_excel(writer, sheet_name="Summary", index=False)
        df_blast_export.to_excel(writer, sheet_name="BLAST_Details", index=False)
        df_qc_export.to_excel(writer, sheet_name="Quality_Control", index=False)

        # --- A MÁGICA DO CONGELAMENTO DE PAINÉIS ---
        # Itera sobre todas as abas criadas e congela a primeira linha
        for sheet_name in writer.sheets:
            writer.sheets[sheet_name].freeze_panes = "A2"
        # -------------------------------------------

        # Ajusta larguras das colunas
        auto_adjust_columns(writer, "Summary")
        auto_adjust_columns(writer, "BLAST_Details")
        auto_adjust_columns(writer, "Quality_Control")

        # Aplica cores baseadas no df_blast já filtrado
        apply_workbook_styling(writer, df_summary, df_blast_export, df_qc_export)

#
# MAIN LOGIC
#

def main():
    # --- 1. EXTRAÇÃO ---
    df_trim = extract_trim_data(snakemake.input.trim_summary)
    df_cons = extract_consensus_data(snakemake.input.consensus_summary)
    df_blast_completo = parse_blast_xmls(snakemake.input.xmls)
    
    # --- 2. PROCESSAMENTO ---
    df_qc = build_qc_df(df_trim, df_cons, df_blast_completo)
    df_summary = build_summary_df(df_qc, df_blast_completo)
    
    # --- 3. APRESENTAÇÃO / EXPORTAÇÃO ---
    export_to_excel(df_summary, df_blast_completo, df_qc, snakemake.output.report)



if __name__ == "__main__":
    main()


