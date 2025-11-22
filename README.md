# YeastProspecta

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


Results directory structure:
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