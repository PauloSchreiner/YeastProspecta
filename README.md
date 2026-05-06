# YeastProspecta 🧬

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**YeastProspecta** is an automated, scalable, and reproducible Snakemake pipeline designed for the comprehensive analysis of Sanger sequencing data. It automates the processing of Forward (F) and Reverse (R) `.ab1` read pairs, from raw sequence trimming to species identification via automated NCBI BLAST searches.

The pipeline applies quality control thresholds, assembles consensus sequences, performs remote alignments, and aggregates all metrics into a highly readable, color-coded Excel report to facilitate taxonomic identification.

---

## 📑 Table of Contents
1. [Workflow Overview](#workflow-overview)
2. [Installation & Dependencies](#installation-dependencies)
3. [Input Data Preparation](#input-data-preparation)
4. [Configuration & Parameters](#configuration-parameters)
5. [Running the Pipeline](#running-the-pipeline)
6. [Outputs & Interpretation](#outputs-interpretation)
7. [Identification Rules Logic](#identification-rules-logic)

---

## <a id="workflow-overview"></a>⚙️ Workflow Overview

The pipeline consists of three major interconnected steps:

1. **Quality Trimming:** Reads raw `.ab1` chromatograms and performs quality trimming using Mott's algorithm based on a user-defined error probability cutoff.
2. **Consensus Assembly:** Aligns the paired trimmed reads (F and R) to generate a high-quality consensus sequence using EMBOSS `merger`.
3. **Taxonomic Identification (BLAST):** Automatically queries the consensus sequences against the specified database (e.g., NCBI `nt`). It retrieves the top hits, extracts binomial species names, filters for type strains, and evaluates the final identification status based on strict identity and consensus score rules.

---

## <a id="installation-dependencies"></a>💻 Installation & Dependencies

YeastProspecta relies on Conda/Mamba for environment management to ensure strict reproducibility. 

### 1. Clone the repository
```bash
git clone [https://github.com/](https://github.com/)[your-username]/YeastProspecta.git
cd YeastProspecta
```

### 2. Dependencies
*[Sugestão: A melhor prática em Snakemake é ter um arquivo `environment.yaml` na raiz do projeto. Se você tiver, a instalação fica apenas isso:]*

```bash
conda env create -f environment.yaml
conda activate yeastprospecta
```

*(If you are installing dependencies manually, ensure the following are available in your PATH: Python 3+, Snakemake, Biopython, Pandas, Openpyxl, EMBOSS `merger`, and NCBI BLAST+ CLI).*

---

## <a id="input-data-preparation"></a>📂 Input Data Preparation

The pipeline expects paired Sanger chromatograms (`.ab1` or `.seq` format) in the input directory.

### Naming Convention
To ensure the pipeline pairs the reads correctly, files must follow a strict naming convention:
`[SampleID]_F.ab1` and `[SampleID]_R.ab1`

*(Example: `AM003_F.ab1` and `AM003_R.ab1`)*

**Utility Script:** A pre-processing script is provided to standardize raw sequencing names coming from the facility. Run the following command to format your files automatically:
```bash
python scripts/utils/rename_samples.py --dir path/to/your/raw_data --execute
```

---

## <a id="configuration-parameters"></a>🔧 Configuration & Parameters

All pipeline parameters are controlled via the `config.yaml` file located in the root directory. You do not need to modify any Python or Shell scripts to adapt the pipeline to your specific organism or amplicon size.

### `config.yaml`
```yaml
directories:
  inputs: "data"                                      # Where the input .ab1 files must be placed
  results: "results"                                  # Results root directory
  trim_fastq:   "results/01_trim/fastq"               # Trimmed sequences (F and R separately)
  trim_reports: "results/01_trim/reports"             # Details of trimming process
  consensus_fasta:   "results/02_consensus/fasta"     # Consensus sequences (F and R reads combined)
  consensus_reports: "results/02_consensus/reports"   # Consensus .txt reports (scores and visual alignment)
  blast_raw: "results/03_blast/raw_hits"              # BLAST xml files 
  summary: "results/summary"                          # Final summaries 

parameters:
  trim_cutoff: 0.05          # Stringency on mott's algorithm for trimming
  blast_database: "nt"       # Target database for BLAST searches (e.g., nt, ITS_RefSeq_Fungi)
  identity_threshold: 99.0   # Minimum identity % required for species identification
  good_consensus_score: 2000 # Threshold for a reliable consensus alignment score (EMBOSS merger)
  max_excel_hits: 20         # Max amount of BLAST hits exported to the BLAST_Details tab
```

### NCBI API Key (Crucial for BLAST)
If you plan to run multiple BLAST jobs simultaneously, an NCBI API Key is required to prevent connection timeouts or temporary IP bans.

1. Generate an API Key at your NCBI Account settings.
2. Export it to your environment before running the pipeline:
```bash
export NCBI_KEY="your_api_key_here"
```
3. Update the `profiles/ncbi_safe/config.yaml` file to allow parallel BLAST requests:
```yaml
resources:
  - ncbi_connection=4  # Set this to the number of parallel requests allowed by your key
```
*Note: If you do not have an API key, keep `ncbi_connection=1` to ensure safe, sequential requests.*

---

## <a id="running-the-pipeline"></a>🚀 Running the Pipeline

Once your data is in the `inputs` folder and configured, run the pipeline using the provided safety profile:

**Dry-run (to check which steps will be executed):**
```bash
snakemake --profile profiles/ncbi_safe -n
```

**Execute pipeline:**
```bash
snakemake --profile profiles/ncbi_safe --cores [number_of_cpu_cores]
```

---

## <a id="outputs-interpretation"></a>📊 Outputs & Interpretation

The pipeline generates intermediate files in the `results/` directory. The primary output for researchers is located at:

📁 **`results/summary/final_report.xlsx`**

This comprehensive report contains three sheets:
1. **Summary:** The definitive identification report. Samples are automatically flagged with conditional coloring:
   * 🟩 **OK:** Solid identification matching a type strain (Identity ≥ `identity_threshold`).
   * 🟨 **Ambiguous:** Multiple species match the criteria.
   * 🟦 **No hit with type:** Has hits meeting the identity threshold, but none are recognized type strains.
   * 🟪 **Below identity threshold:** No hits with any GenBank entries met the required threshold. *(You might have found something new!)*
   * 🟥 **Bad quality:** The consensus score is below the `good_consensus_score`, making identification unreliable.
2. **BLAST_Details:** Contains the top raw hits from the database for every sample, including E-values and Coverage, limited by your `max_excel_hits` setting.
3. **Quality_Control:** Detailed metrics for sequence lengths (F, R, and Consensus) and EMBOSS alignment scores, visually color-coded based on your chosen consensus baseline.

---

## <a id="identification-rules-logic"></a>🎛️ Identification Rules Logic

The taxonomic identification status is determined dynamically based on the parameters set in your `config.yaml`. The pipeline evaluates each sample top-to-bottom using the following logical hierarchy:

1. **Rule 1 (Quality Check):** Is the `Consensus_score` strictly lower than `good_consensus_score`? 
   * *If Yes ->* Flags as **Bad quality**.
2. **Rule 2 (Ambiguity Check):** Are there multiple valid type strain species matching the criteria? 
   * *If Yes ->* Flags as **Ambiguous**.
3. **Rule 3 (Identity Check):** Is the highest identity found across ALL hits strictly lower than the `identity_threshold`? 
   * *If Yes ->* Flags as **Below identity threshold**.
4. **Rule 4 (Type Strain Check):** Did the sample pass the identity threshold but fail to match against any recognized "type material"? 
   * *If Yes ->* Flags as **No hit with type**.
5. **Rule 5 (Success):** If none of the above rules are triggered, the sample is verified.
   * *Result ->* Flags as **OK**.

---

**Maintained by:** *Paulo Schreiner / Laboratório de Biologia Computacional e Molecular (LBCM)* **Contact:** *paulorangelschreiner@gmail.com*