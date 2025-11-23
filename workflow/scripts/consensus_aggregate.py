import os
import re
import pandas as pd

data_rows = []

for report_path in snakemake.input.reports:
    sample = os.path.basename(report_path).removesuffix("_consensus_report.txt")
    
    metrics = {"Sample": sample}

    with open(report_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("# Length:"):
                metrics["Length"] = line.split(":", 1)[1].strip()
            elif line.startswith("# Identity:"):
                m = re.search(r"\((\s*[\d\.]+%)\)", line)
                metrics["Identity"] = m.group(1).strip() if m else ""
            elif line.startswith("# Similarity:"):
                m = re.search(r"\((\s*[\d\.]+%)\)", line)
                metrics["Similarity"] = m.group(1).strip() if m else ""
            elif line.startswith("# Gaps:"):
                m = re.search(r"\((\s*[\d\.]+%)\)", line)
                metrics["Gaps"] = m.group(1).strip() if m else ""
            elif line.startswith("# Score:"):
                m = re.search(r"# Score:\s*([\d\.]+)", line)
                metrics["Score"] = m.group(1) if m else ""

    data_rows.append(metrics)


df = pd.DataFrame(data_rows)
df.to_csv(snakemake.output[0], sep='\t', index=False)