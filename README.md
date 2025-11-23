# YeastProspecta

## How to use

1. Quick start:
- install necessary dependencies
- export NCBI API key
  - (how to generate one), (change number of parallel jobs in profile)
- snakemake command


---

If you have an NCBI API key, run this command right before running the pipeline:

```bash 
export NCBI_KEY="your_key"
```

This way you can safely set your profile to run 4 blast jobs at once by changing the following line in `/profiles/ncbi_safe`

```bash 
resources:
  - ncbi_connection=4
```

If you do not have an NCBI API key, this will likely lead to an error due to excessive requests. 

To run the pipeline, simply run:
```bash 
snakemake --profile profiles/ncbi_safe
```



---

Pipeline intended steps:
### 1. Trim
**Input**: 
- raw .ab1 files

**Output**: 
- trimmed reads (fastq)
- trim reports per read (csv)
- an aggregated trim report (csv)

### 2. Consensus
**Input**:
- trimmed reads (fastq)

**Output**: 
- consensus sequence ()
- consensus report per read ()
- consensus metrics summary (csv)

### 3. BLAST
**Input**:
- consensus sequence

**Output**:
- BLAST results per read (use XML instead of tabular)
- BLAST full results 
- BLAST interpretation (processed full results)

---


### config.yaml
```yaml
# Do not add "/" at the end of paths.
directories:
  inputs: "data"

  results: "results"
  trim_fastq:   "results/01_trim/fastq"
  trim_reports: "results/01_trim/reports"
  consensus_fasta:   "results/02_consensus/fasta"
  consensus_reports: "results/02_consensus/reports"
  blast_raw: "results/03_blast/raw_hits"
  summary: "results/summary/"

parameters:
  trim_cutoff: 0.05
```

--- 

### profile/ncbi_safe/config.yaml
```yaml
cores: 4

# Ensures there won't be paralellism in BLASTn API, 
# which would result in banning
default-resources:
  - ncbi_connection=1

# UX settings
keep-going: True
printshellcmds: True
```

---

### Results directory structure:
```bash
results/
├── 01_trim/
│   ├── fastq/
│   └── reports/
│
├── 02_consensus/
│   ├── fasta/
│   └── reports/         
│
├── 03_blast/
│   └── raw_hits/
│
└── summary/
    ├── trim_metrics.tsv
    ├── consensus_metrics.tsv
    ├── blast_full.tsv
    └── blast_interpreted.tsv
```