import pandas as pd
from Bio import SearchIO
import re
import csv
import os 
from openpyxl.formatting.rule import CellIsRule
from openpyxl.styles import PatternFill

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


def auto_adjust_columns(writer, sheet_name):
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
        
        # Set width to max_length + padding (e.g., 2)
        # We cap it at 60 characters so really long descriptions don't break the UI
        adjusted_width = min(max_length + 2, 60)
        worksheet.column_dimensions[column_letter].width = adjusted_width


def extract_trim_data(trim_path:str) -> dict:
    """
    Return:
        dict or list with:
        sample, Post-trim length of F read, Post-trim length of R read
        for each sample
    """
    results = {}

    with open(trim_path, mode="r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            raw_sample_name = row['sample']
            length = int(float(row['post_trim_length']))
            
            if "_" not in raw_sample_name:
                print(f"Warning: Unexpected sample name format: {raw_sample_name}")
                continue

            sample_id, direction = raw_sample_name.rsplit('_', 1)
            
            # Initialize sample entry if it has not been iterated over yet
            if sample_id not in results:
                results[sample_id] = {"F_len":0, "R_len":0}

            if direction == "F":
                results[sample_id]["F_len"] = length
            elif direction == "R":
                results[sample_id]["R_len"] = length
    
    return results


def extract_consensus_data(consensus_path:str) -> dict:
    """
    Return:
        dict with:
        sample, consensus_score, consensus_length
        for each sample
    """
    results = {}

    with open(consensus_path, mode="r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            sample = row["Sample"]
            consensus_length = int(float(row["Length"]))
            try:
                consensus_score = int(float(row["Score"]))
            except ValueError:
                consensus_score = 0
            results[sample] = {"Consensus_length": consensus_length, 
                               "Consensus_score": consensus_score}

    return results


def get_total_hits_per_sample(xml_path_list:list) -> dict:
    """
    Input: snakemake.input.xmls
    Return:
        dict or list with
        sample, Total number of >99% identity hits
        per sample
    """
    results = {}
    
    for xml_file in xml_path_list:
        sample = os.path.basename(xml_file).rsplit("_", 1)[0]

        hits_over_99 = 0

        qresult = SearchIO.read(xml_file, "blast-xml")

        for hit in qresult:

            hsp = hit[0]

            identity_pct = (hsp.ident_num / hsp.aln_span) * 100
            if identity_pct > 99.0:
                hits_over_99 += 1
        
        results[sample] = {"Hits_>99_id_pct": hits_over_99}
        
    return results


def parse_blast_xmls(xml_path_list:list) -> pd.DataFrame:
    """
    Raw BLAST Sheet (for each of the top 10 hits):
    - Sample Name
    - Hit #
    - Title
    - Accession
    - Identity
    - Coverage
    - E-value
    """
    rows = []

    for xml_file in xml_path_list:
        sample = os.path.basename(xml_file).rsplit("_", 1)[0] 

        qresult = SearchIO.read(xml_file, "blast-xml")

        hit_index = 0 
        for hit in qresult:
            # Limit to only 10 hits 
            hit_index += 1
            if hit_index > 10:
                break

            hsp = hit[0]

            identity_pct = round((hsp.ident_num / hsp.aln_span) * 100, 3)

            query_coverage = round((hsp.query_span / qresult.seq_len) * 100, 3)
        
            row_dict = {
                "Sample": sample.zfill(3),
                "Hit #": hit_index,
                "Accession": hit.accession,
                "Identity": identity_pct,
                "Coverage": query_coverage,
                "E-value":hsp.evalue,
                "Title": hit.description
            }

            rows.append(row_dict)

    return pd.DataFrame(rows)
                

def generate_summary(df_qc:pd.DataFrame, df_blast:pd.DataFrame) -> pd.DataFrame:
    """
    Determine status based on:
    Consensus_score,
    Title of each hit,
    """
    
    rows = []

    for sample_name in df_qc.index:

        group = df_blast[df_blast["Sample"] == sample_name] 
        group = group.sort_values("Identity", ascending=False).copy()

        if group.empty:
            top_species = "None"
            top_identity = 0
            other_species = ""
        else: 
            top_hit_row = group.iloc[0]
            top_identity = round(top_hit_row["Identity"], 3)
            
            top_species = extract_species_from_title(top_hit_row["Title"])
            
            # --- START OF OTHER SPECIES LOGIC --- #
            clean_top = top_species.replace('[', '').replace(']', '')
            seen_species = {clean_top}
            
            other_species_list = []

            for index, row in group.iterrows():
                
                current_species = extract_species_from_title(row["Title"])
                clean_current = current_species.replace('[', '').replace(']', '')
                
                if clean_current not in seen_species:
                    current_identity = round(row["Identity"], 3)
                    
                    formatted_string = f"{current_species} ({current_identity}%)"
                    other_species_list.append(formatted_string)
                    
                    seen_species.add(clean_current)
            
            other_species = "; ".join(other_species_list)
            # --- END OF OTHER SPECIES LOGIC --- #

        consensus_score = df_qc.loc[sample_name, "Consensus_score"]
        status = "Pending"

        if consensus_score < 2000:
            status = "Bad quality"
        
        elif other_species:
            status = "Ambiguous"

        elif df_qc.loc[sample_name, "Hits_>99_id_pct"] == 0 and top_identity < 99:
            status = "No match"

        else:
            status = "OK"

        rows.append({
            "Sample": sample_name.zfill(3),
            "Status": status,
            "Top hit": top_species,
            "Identity": top_identity,
            "Other possible species": other_species
        })
    
    return pd.DataFrame(rows)


from openpyxl.formatting.rule import CellIsRule
from openpyxl.styles import PatternFill

def apply_workbook_styling(writer, df_summary, df_blast, df_qc):
    """
    Applies all conditional formatting and styling to the entire Excel workbook.
    """
    # 1. DEFINE PRETTY HUES
    colors = {
        "green":      PatternFill(start_color='C6EFCE', end_color='C6EFCE', fill_type='solid'), # Status OK / High Score
        "lime":       PatternFill(start_color='E3EDB5', end_color='E3EDB5', fill_type='solid'), # Good Score
        "yellow":     PatternFill(start_color='FFEB9C', end_color='FFEB9C', fill_type='solid'), # Ambiguous / Mid Score
        "red":        PatternFill(start_color='FFC7CE', end_color='FFC7CE', fill_type='solid'), # Bad quality / Low Score
        "purple":     PatternFill(start_color='E4CCFF', end_color='E4CCFF', fill_type='solid'), # No match
        "light_gray": PatternFill(start_color='F2F2F2', end_color='F2F2F2', fill_type='solid'), # BLAST grouping
        "white":      PatternFill(start_color='FFFFFF', end_color='FFFFFF', fill_type='solid')  # BLAST grouping
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
    
    # Identify the 'Consensus_score' column
    score_col_letter = None
    for cell in ws_qc[1]:
        if cell.value == "Consensus_score":
            score_col_letter = cell.column_letter
            break
    
    if score_col_letter:
        qc_range = f"{score_col_letter}2:{score_col_letter}{len(df_qc) + 1}"
        
        # Rules must be added from highest to lowest if using 'greaterThan' 
        # to ensure the first true condition met is the one Excel keeps.
        ws_qc.conditional_formatting.add(qc_range, CellIsRule(operator='greaterThan', formula=['2250'], fill=colors["green"]))
        ws_qc.conditional_formatting.add(qc_range, CellIsRule(operator='greaterThan', formula=['2000'], fill=colors["lime"]))
        ws_qc.conditional_formatting.add(qc_range, CellIsRule(operator='greaterThan', formula=['1750'], fill=colors["yellow"]))
        ws_qc.conditional_formatting.add(qc_range, CellIsRule(operator='lessThanOrEqual', formula=['1750'], fill=colors["red"]))


def main():

    # QC SHEET
    trim_dict = extract_trim_data(snakemake.input.trim_summary)
    consensus_dict = extract_consensus_data(snakemake.input.consensus_summary)
    blast_qc_dict = get_total_hits_per_sample(snakemake.input.xmls)

    df_trim = pd.DataFrame.from_dict(trim_dict, orient='index')
    df_cons = pd.DataFrame.from_dict(consensus_dict, orient='index')
    df_blast = pd.DataFrame.from_dict(blast_qc_dict, orient='index')

    df_qc_final = df_trim.join([df_cons, df_blast], how='outer')
    df_qc_final.index.name = "Sample"
    df_qc_final.index = df_qc_final.index.astype(str).str.zfill(3)
    df_qc_final.sort_index(inplace=True)

    df_qc_export = df_qc_final.reset_index() # This avoids using a different font

    # RAW BLAST SHEET
    df_blast_details_final = parse_blast_xmls(snakemake.input.xmls)

    # SUMMARY SHEET
    df_summary = generate_summary(df_qc_final, df_blast_details_final)

    # SORT SHEETS
    df_qc_final.sort_index(inplace=True)
    df_blast_details_final.sort_values(by=["Sample", "Hit #"], inplace=True)
    df_summary.sort_values(by="Sample", inplace=True)

    # WRITE SHEETS TO EXCEL
    output_file = snakemake.output.report
    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
        df_summary.to_excel(writer, sheet_name="Summary", index=False)
        auto_adjust_columns(writer, "Summary")

        df_blast_details_final.to_excel(writer, sheet_name="BLAST_Details", index=False)
        auto_adjust_columns(writer, "BLAST_Details")

        df_qc_export.to_excel(writer, sheet_name="Quality_Control", index=False)
        auto_adjust_columns(writer, "Quality_Control")

        apply_workbook_styling(writer, df_summary, df_blast_details_final, df_qc_export)


main()