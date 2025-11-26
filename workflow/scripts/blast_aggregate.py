import pandas as pd
from Bio import SearchIO


for xml_file in snakemake.input.xmls:
    qresult = SearchIO.read(xml_file, "blast-xml")
    
    hits_over_99 = 0

    for hit in qresult:

        hsp = hit[0] 

        # Get number of hits with more than 99% of identity
        identity_pct = (hsp.ident_num / hsp.aln_span) * 100
        if identity_pct > 99.0:
            hits_over_99 += 1
    
    
    results = {
        "hits_over_99" : hits_over_99
    }

