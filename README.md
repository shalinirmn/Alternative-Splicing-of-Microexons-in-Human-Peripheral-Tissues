# Alternative Splicing of Microexons in Human Peripheral Tissues

[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.05.13.725022-B31B1B?logo=arxiv&logoColor=white)](https://doi.org/10.64898/2026.05.13.725022)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Pipeline](https://img.shields.io/badge/Pipeline-Bash%20%7C%20VAST--TOOLS%20%7C%20R-2ea44f.svg)](#pipeline-overview)

A reproducible bioinformatics pipeline for quantifying and characterizing **microexon (3–27 nt) alternative splicing** from bulk RNA-seq cohorts across four human peripheral tissues: **liver, lung, kidney, and colon**.

---

## Pipeline Overview

```mermaid
flowchart LR
    A[Public SRA/ENA<br/>Bulk RNA-seq] -->|01_fetch_fastq.sh<br/>wget| B[Raw Paired-End<br/>FASTQ]
    B -->|02_qc_and_trim.sh<br/>FastQC + fastp + MultiQC| C[Trimmed & QC'd<br/>Reads]
    C -->|03_vast_align.sh<br/>vast-tools align| D[Chunked Junction<br/>Alignments]
    D -->|04_vast_combine.sh<br/>vast-tools combine| E[INCLUSION_LEVELS_FULL<br/>PSI Matrix]
    E -->|R Downstream<br/>Analysis| F[Microexon Filtering,<br/>ΔPSI & Figures]
```

---

## Repository Structure

```text
Alternative-splicing-of-microexons-in-peripheral-tissues/
├── README.md
├── LICENSE
├── environment.yml                 # Conda environment for QC & preprocessing
├── config/
│   └── samples_lung_sclc.tsv       # Cohort sample manifest (GROUP, SAMPLE, READ1, READ2, URLs)
├── scripts/
│   ├── 01_fetch_fastq.sh           # Automated bulk FASTQ download via wget
│   ├── 02_qc_and_trim.sh           # FastQC, fastp adapter/quality trimming, and MultiQC
│   ├── 03_vast_align.sh            # VAST-TOOLS align wrapper (SLURM array or sequential)
│   └── 04_vast_combine.sh          # VAST-TOOLS combine to generate merged PSI inclusion tables
├── example_data/                   # Subsampled toy data for local verification
└── results/
    └── figures/                    # Output plots
```

---

## Prerequisites & Setup

### 1. Preprocessing & QC Environment
```bash
git clone [https://github.com/shalinirmn/Alternative-splicing-of-microexons-in-peripheral-tissues.git](https://github.com/shalinirmn/Alternative-splicing-of-microexons-in-peripheral-tissues.git)
cd Alternative-splicing-of-microexons-in-peripheral-tissues

conda env create -f environment.yml
conda activate microexon-splicing
```

### 2. VAST-TOOLS & VASTDB Setup
Splicing quantification relies on `vast-tools` and the human `VASTDB` library (`hg38`). Install `vast-tools` and download the species database following the official [vast-tools repository instructions](https://github.com/vastgroup/vast-tools).

---

## Step-by-Step Execution

### Step 1: Bulk Data Retrieval
Download paired-end bulk RNA-seq FASTQ files defined in a cohort manifest (`config/samples_<tissue>.tsv`):
```bash
bash scripts/01_fetch_fastq.sh config/samples_lung_sclc.tsv data/raw/lung
```

### Step 2: Quality Control & Adapter Trimming
Run pre-trimming `FastQC`, adapter/low-complexity trimming via `fastp` (`-l 36`, `--average_qual 20`), and aggregate reports with `MultiQC`:
```bash
bash scripts/02_qc_and_trim.sh \
  config/samples_lung_sclc.tsv \
  data/raw/lung \
  results/qc_lung \
  24
```

### Step 3: Splicing Alignment (`vast-tools align`)
Align trimmed paired-end reads against `VASTDB` (`hg38`). Supports both sequential execution and SLURM job arrays (`sbatch --array=1-N`):
```bash
bash scripts/03_vast_align.sh \
  config/samples_lung_sclc.tsv \
  results/qc_lung/fastp_out \
  results/vast_lung \
  /path/to/VASTDB \
  hg38 \
  24
```

### Step 4: Merge Cohort PSI Tables (`vast-tools combine`)
Combine per-sample junction alignments into the master `INCLUSION_LEVELS_FULL` PSI matrix:
```bash
bash scripts/04_vast_combine.sh \
  results/vast_lung \
  /path/to/VASTDB \
  hg38 \
  24
```

### Step 5: Downstream R Analysis *(In Progress)*
R scripts for microexon length filtering (`3–27 nt`), read-coverage quality filtering, differential splicing ($\Delta\text{PSI}$), cross-tissue comparison, and figure generation will be documented here.

---

## Citation

If you use this pipeline or data in your work, please cite:

> Raman, S., & Gupta, I. (2026). *In silico investigation of alternative splicing of microexons in human peripheral tissues*. bioRxiv. https://doi.org/10.64898/2026.05.13.725022
