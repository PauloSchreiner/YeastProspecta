# YeastProspecta 🧬

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**YeastProspecta** is an automated, scalable, and reproducible Snakemake pipeline designed for the comprehensive analysis of Sanger sequencing data. It automates the processing of Forward (F) and Reverse (R) `.ab1` read pairs, from raw sequence trimming to species identification via automated NCBI BLAST searches.

The pipeline applies quality control thresholds, assembles consensus sequences, performs remote alignments, and aggregates all metrics into a highly readable, color-coded Excel report to facilitate taxonomic identification.

---

## 📑 Table of Contents
1. [Workflow Overview](#-workflow-overview)
2. [Installation & Dependencies](#-installation--dependencies)
3. [Input Data Preparation](#-input-data-preparation)
4. [Configuration](#-configuration)
5. [Running the Pipeline](#-running-the-pipeline)
6. [Outputs & Interpretation](#-outputs--interpretation)
7. [Citation](#-citation)

---

## ⚙️ Workflow Overview

The pipeline consists of three major interconnected steps:

1. **Quality Trimming:** Reads raw `.ab1` chromatograms and performs quality trimming using Mott's algorithm based on a user-defined error probability cutoff.
2. **Consensus Assembly:** Aligns the paired trimmed reads (F and R) to generate a high-quality consensus sequence using EMBOSS `merger`.
3. **Taxonomic Identification (BLAST):** Automatically queries the consensus sequences against the NCBI `nt` database. It retrieves the top hits, extracts binomial species names, filters for type strains, and evaluates the final identification status based on strict identity and consensus score rules.

---

## 💻 Installation & Dependencies

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

## 📂 Input Data Preparation

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

## 🔧 Configuration

All pipeline parameters are controlled via the `config.yaml` file located in the root directory.

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
  blast_database: "nt"       # NCBI database for BLAST searches
  identity_threshold: 99.0   # Identity threshold for species identification via D1/D2 
  good_consensus_score: 2000 # Arbitrary threshold for a "good" consensus score (EMBOSS merger)
  max_excel_hits: 20         # Max amount of BLAST hits to be included in BLAST_Details tab in final_report.xlsx


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

## 🚀 Running the Pipeline

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

## 📊 Outputs & Interpretation

The pipeline generates intermediate files in the `results/` directory. The primary output for researchers is located at:

📁 **`results/summary/final_report.xlsx`**

This comprehensive report contains three sheets:
1. **Summary:** The definitive identification report. Samples are automatically flagged with conditional coloring:
   * 🟩 **OK:** Solid identification matching a type strain (>98.90% identity).
   * 🟨 **Ambiguous:** Multiple species match the criteria.
   * 🟧 **No hit with type:** Has hits with >98.90% identity, but none are type strains.
   * 🟪 **No hit above 99:** No hits with any GenBank entries within the threshold. *(You might have found something new!)*
   * 🟥 **Bad quality:** The consensus score is too low for reliable identification.
2. **BLAST_Details:** Contains the top 20 raw hits from the NCBI database for every sample, including E-values and Coverage.
3. **Quality_Control:** Detailed metrics for sequence lengths (F, R, and Consensus) and EMBOSS alignment scores, as well as the total count of hits with more than 99% identity.

---

## 🎛️ Customizing Identification Rules

The taxonomic identification status (OK, Ambiguous, etc.) is determined by a strict set of sequential rules. If you need to adjust the thresholds for your specific organism or laboratory standards (e.g., lowering the identity threshold to 98.5% or changing the required consensus score), you can easily modify these parameters in the `scripts/final_report.py` file.

Look for the following code block inside the `build_summary_df` function:

```python
    conditions = [
        df_final['Consensus_score'] < 2000,           # Rule 1: Poor sequence quality
        df_final['Other possible species'] != "",     # Rule 2: More than one valid species found
        df_final['Highest identity overall'] < 99,    # Rule 3: No good hit found at all
        df_final['Top hit (type strain)'] == "None"   # Rule 4: Good hit exists, but no type found 
    ]
    
    labels = [
        "Bad quality",
        "Ambiguous",
        "No hit above 99",
        "No hit with type"
    ]
```

*Note: The script evaluates these rules from top to bottom. The first rule that evaluates to `True` assigns the status. If no rules are triggered, the sample is marked as **OK**.*

---

**Maintained by:** *Paulo Schreiner / Laboratório de Biologia Computacional e Molecular (LBCM)* **Contact:** *paulorangelschreiner@gmail.com*