# YeastProspecta 🧬

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**YeastProspecta** is an automated, scalable, and reproducible Snakemake pipeline designed for the comprehensive analysis of Sanger sequencing data in bioprospecting projects. It automates the processing of Forward (F) and Reverse (R) `.ab1` read pairs, from raw sequence trimming to species identification via automated NCBI BLAST searches.

The pipeline applies quality control thresholds, assembles consensus sequences, performs remote alignments against the NCBI database, and aggregates all results and metrics into a readable, color-coded Excel report to facilitate taxonomic identification.

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

The pipeline consists of three major steps:

1. **Quality Trimming:** Reads raw `.ab1` chromatograms and performs quality trimming using Mott's algorithm based on a user-defined error probability cutoff.
2. **Consensus Assembly:** Aligns the paired trimmed reads (F and R) to generate a high-quality consensus sequence using EMBOSS `merger`.
3. **Taxonomic Identification (BLAST):** Automatically queries the consensus sequences against the specified database (e.g., NCBI `nt`). It retrieves the top 500 hits, extracts binomial species names, filters for type strains, and evaluates the final identification status based on strict identity and consensus score rules. Other details and metrics are displayed in `final_report.xlsx`.

---

## <a id="installation-dependencies"></a>💻 Installation & Dependencies

YeastProspecta relies on Conda/Mamba for environment management to ensure strict reproducibility. 

### 1. Clone the repository
```bash
git clone https://github.com/PauloSchreiner/YeastProspecta.git
cd YeastProspecta
```

### 2. Dependencies

To install the necessary dependencies, run these commands from the project root directory:

```bash
conda env create -f environment.yaml
conda activate yeastprospecta
```

And you're done!

---

## <a id="input-data-preparation"></a>📂 Input Data Preparation

The pipeline expects paired Sanger chromatograms (`.ab1` or `.seq` format) in the input directory.

### Naming Convention
To ensure the pipeline pairs the reads correctly, files must follow a strict naming convention:
`[SampleID]_F.ab1` and `[SampleID]_R.ab1`

*(Example: `AM003_F.ab1` and `AM003_R.ab1`)*

**Utility Script:** A pre-processing script is provided to standardize raw sequencing names coming straight from the sequencing machine. If you wish to use it, run:
```bash
python scripts/utils/rename_samples.py --dir path/to/your/raw_data --execute
```

---

## <a id="configuration-parameters"></a>🔧 Configuration & Parameters

All pipeline parameters are controlled via the `config.yaml` file located in the root directory, so that you can easily customize them to adapt the pipeline to your needs. 

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

*Keep in mind that the `good_consensus_score` value is totally arbitrary and may vary depending on sequence length. It refers to the consensus score reported by **EMBOSS merger**, and is used only for quality control.*

### NCBI API Key (Crucial for BLAST)
If you plan to run multiple BLAST jobs simultaneously, an **NCBI API Key** is required to prevent connection timeouts or temporary IP bans.

1. Generate an API Key at your NCBI Account settings.
2. Export it to your environment **before running the pipeline**:
```bash
export NCBI_KEY="your_api_key_here"
```

---

## <a id="running-the-pipeline"></a>🚀 Running the Pipeline

Once your data is in the `inputs` folder and you have configured the parameters, you can choose between two ways of running the pipeline:

**1. Standard Run (No NCBI API Key, slower)**
If you do not have an NCBI API key, you must limit the pipeline to 1 connection to the NCBI server at a time. The command below allows Snakemake to use multiple CPU cores for trimming and consensus, but forces the BLAST step to run sequentially:

```bash
snakemake --cores all --resources ncbi_connection=1
```

You can alter the number of cores by changing `all` to any number.

**2. Accelerated Run (NCBI API Key Required)**
If you have an NCBI API key, you can speed up the pipeline by allowing multiple parallel BLAST requests. (Instructions for obtaining an API key can be found in [Configuration & Parameters](#configuration-parameters).)

First, export your key to your environment:
```bash
export NCBI_KEY="your_api_key_here"
```
Then, run the pipeline increasing the `ncbi_connection` resource (e.g., to 4 simultaneous connections):
```bash
snakemake --cores all --resources ncbi_connection=4
```

(Tip: To perform a dry-run and see which steps will be executed without actually running them, add the `-n` flag to any of the commands above).

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

These rules can be modified by altering the following code, found at line 227 of the `final_report.py` script:

```python
conditions = [
        df_final['Consensus_score'] < good_cons_score,                      # Rule 1: Poor sequence quality
        df_final['Other possible species'] != "",                           # Rule 2: More than one valid species found
        df_final['Highest identity overall'] < id_thresh,                   # Rule 3: No hit above threshold
        df_final['Top hit (type strain)'] == "None"                         # Rule 4: No type found above threshold
    ]
    
    labels = [
        "Bad quality",
        "Ambiguous",
        "No hit above threshold",
        "No hit with type"
    ]

# Note that the script evaluates these rules from top to bottom. 
# The first rule that evaluates to True assigns the status. 
# If no rules are triggered, the sample is marked as OK.
```
---

**Created and maintained by:** *Paulo Schreiner / Prof. Dr. Diego Bonatto - Laboratório de Biologia Computacional e Molecular (LBCM) - UFRGS* **Contact:** *paulorangelschreiner@gmail.com*