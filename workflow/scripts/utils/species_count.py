import os
import re
import pandas as pd


def gerar_relatorio_especies(caminho_excel, caminho_saida="contagem_especies.xlsx"):
    if not os.path.exists(caminho_excel):
        print(f"❌ Erro: O arquivo '{caminho_excel}' não foi encontrado.")
        return

    print(f"📊 Lendo dados de: {caminho_excel}...")

    # 1. Ler apenas a aba 'Summary' do seu Excel de resultados
    df = pd.read_excel(caminho_excel, sheet_name="Summary")

    # 2. Filtrar amostras com status 'OK', 'Many species' ou '4 or more mutations'
    df["Status"] = df["Status"].astype(str).str.strip()
    status_validos = ["OK", "Many species", "4 or more mutations"]  # <-- ATUALIZADO AQUI
    df_filtrado = df[df["Status"].isin(status_validos)].copy()

    if df_filtrado.empty:
        print(
            "⚠️ Aviso: Nenhuma amostra com os status selecionados foi encontrada."
        )
        return

    # 3. Extrair dinamicamente o prefixo do laboratório (ex: MR, SF, AM, LS)
    df_filtrado["Prefixo"] = df_filtrado["Sample"].str.extract(r"^([A-Za-z]+)")
    df_filtrado["Prefixo"] = df_filtrado["Prefixo"].fillna("Outros")

    # 4. Criar a tabela de contagem cruzada (Crosstab)
    tabela_contagem = pd.crosstab(
        index=df_filtrado["Top hit (type)"],
        columns=df_filtrado["Prefixo"],
        values=df_filtrado["Sample"],
        aggfunc="count",
    ).fillna(0).astype(int)

    # 5. Adicionar a coluna de Total horizontal (soma das linhas)
    tabela_contagem["Total"] = tabela_contagem.sum(axis=1)

    # 6. Ordenar decrescente pelo Total
    tabela_contagem.sort_values(by="Total", ascending=False, inplace=True)

    # 7. Ajustar os nomes e a ordem das colunas para o formato final
    tabela_contagem = tabela_contagem.reset_index()
    tabela_contagem.rename(columns={"Top hit (type)": "Espécie"}, inplace=True)

    # Organiza os cabeçalhos das colunas de contagem em ordem alfabética
    colunas_grupos = sorted(
        [c for c in tabela_contagem.columns if c not in ["Espécie", "Total"]]
    )
    colunas_finais = ["Espécie"] + colunas_grupos + ["Total"]
    tabela_contagem = table_contagem = tabela_contagem[colunas_finais]

    # Renomear os cabeçalhos para o formato "Contagem no XX"
    renomear_colunas = {grupo: f"Contagem no {grupo}" for grupo in colunas_grupos}
    tabela_contagem.rename(columns=renomear_colunas, inplace=True)

    # 8. Salvar o resultado em um novo arquivo Excel
    with pd.ExcelWriter(caminho_saida, engine="openpyxl") as writer:
        tabela_contagem.to_excel(
            writer, sheet_name="Resumo_Especies", index=False
        )

        # Ajuste automático das larguras das colunas
        ws = writer.sheets["Resumo_Especies"]
        for col in ws.columns:
            max_len = max(len(str(cell.value or "")) for cell in col)
            col_letter = col[0].column_letter
            ws.column_dimensions[col_letter].width = max(max_len + 3, 12)

    print(f"✨ Relatório salvo com sucesso em: {caminho_saida}")

    # Exibe uma prévia do resultado final no terminal
    print("\n📈 Prévia do Relatório (incluindo variantes/mutações):")
    print(tabela_contagem.to_string(index=False))


if __name__ == "__main__":
    # Ajuste aqui para o caminho real do arquivo gerado pela sua pipeline principal
    arquivo_entrada = "results/final_report.xlsx"

    gerar_relatorio_especies(arquivo_entrada)
