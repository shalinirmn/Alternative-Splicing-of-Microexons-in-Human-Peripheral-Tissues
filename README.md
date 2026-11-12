# Alternative Splicing of Microexons in Human Peripheral Tissues

[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.05.13.725022-B31B1B?logo=arxiv&logoColor=white)](https://doi.org/10.64898/2026.05.13.725022)

A reproducible bioinformatics pipeline for quantifying and characterizing **microexon (3–27 nt) alternative splicing** from bulk RNA-seq cohorts across four human peripheral tissues: **lung, kidney, colon, and liver**.

---

## Pipeline Overview

```mermaid
flowchart LR
    A[Public GEO / SRA<br/>Bulk RNA-seq FASTQs] -->|01_qc_and_trim.sh<br/>FastQC + fastp + MultiQC| B[Trimmed & QC'd<br/>Reads]
    B -->|02_vast_align.sh<br/>vast-tools align| C[Chunked Junction<br/>Alignments]
    C -->|03_vast_combine.sh<br/>vast-tools combine| D[INCLUSION_LEVELS_FULL<br/>PSI Matrix]
    D -->|R Downstream<br/>Analysis| E[Microexon Filtering,<br/>ΔPSI & Figures]
```

---

## Datasets & Project Accessions

Publicly available bulk RNA-seq datasets spanning ten disease cohorts across four human peripheral tissues (**lung, kidney, colon, and liver**) were retrieved from NCBI BioProject and Gene Expression Omnibus (GEO). Sample-to-accession mappings and clinical group labels (`Disease` vs. `Control`) are documented in [`config/sample_metadata.tsv`](config/sample_metadata.tsv).

| Tissue | Disease Cohort | Abbr. | Accession ID | Repository |
| :--- | :--- | :--- | :--- | :--- |
| **Lung** | Small Cell Lung Cancer | SCLC | [`PRJNA257389`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA257389) | NCBI BioProject |
| **Lung** | Idiopathic Pulmonary Fibrosis | IPF | [`PRJNA358081`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA358081) | NCBI BioProject |
| **Lung** | Chronic Obstructive Pulmonary Disease | COPD | [`PRJNA24581`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA24581) | NCBI BioProject |
| **Lung** | Coronavirus Disease 2019 | COVID-19 | [`GSE183533`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183533) | NCBI GEO |
| **Kidney** | Diabetic Nephropathy | DNE | [`PRJNA595590`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA595590) | NCBI BioProject |
| **Kidney** | Renal Cell Carcinoma | RCC | [`PRJNA704974`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA704974) | NCBI BioProject |
| **Colon** | Colorectal Cancer | CC | [`PRJNA218851`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA218851) | NCBI BioProject |
| **Colon** | Ulcerative Colitis | UC | [`PRJNA1019284`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1019284) | NCBI BioProject |
| **Liver** | Hepatic Cellular Carcinoma | avHCC; eHCC; LC; DN; CH | [`PRJNA471801`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA471801) | NCBI BioProject |
| **Liver** | Non Alcoholic Fatty Liver Disease | NAFLD | [`PRJNA558102`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA558102) | NCBI BioProject |

---

## Repository Structure


## Repository Structure

```text
Alternative-splicing-of-microexons-in-peripheral-tissues/
├── README.md
├── LICENSE
├── environment.yml                     # Conda environment for QC & preprocessing
├── scripts/
│   ├── 01_qc_and_trim.sh               # FastQC, fastp adapter/quality trimming, and MultiQC
│   ├── 02_vast_align.sh                # VAST-TOOLS align wrapper across trimmed FASTQ pairs
│   ├── 03_vast_combine.sh              # VAST-TOOLS combine to generate merged PSI inclusion tables
│   └── 04_differential_microexons.R    # Wilcoxon rank-sum differential PSI & exon length stratification
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
Splicing quantification relies on `vast-tools` and the human `VASTDB` library (`Hsa`).
1. Install `vast-tools` following the official [vast-tools repository instructions](https://github.com/vastgroup/vast-tools).
2. Download the human species database (`vastdb.hsa.*.tar.gz` or `vastdb.hs2.*.tar.gz`) directly from the [CRG VASTDB Library Archive](https://vastdb.crg.eu/libs/) and extract it into your `vast-tools/VASTDB/` directory.
---

## Step-by-Step Execution

### Step 1: Quality Control & Adapter Trimming
After retrieving paired-end FASTQs from the project accessions above, run `FastQC`, adapter/low-complexity trimming via `fastp` (`-l 36`, `--average_qual 20`), and aggregate reports with `MultiQC`:
```bash
bash scripts/01_qc_and_trim.sh /path/to/raw_fastq results/qc_out
```

### Step 2: Splicing Alignment (`vast-tools align`)
Align trimmed paired-end reads (`*_1_trim.fq.gz` and `*_2_trim.fq.gz`) against `VASTDB` (`Hs2`):
```bash
bash scripts/02_vast_align.sh results/qc_out/fastp_out vast_align
```

### Step 3: Merge Cohort PSI Tables (`vast-tools combine`)
Combine per-sample junction alignments into the master `INCLUSION_LEVELS_FULL` PSI matrix:
```bash
bash scripts/03_vast_combine.sh vast_align Hsa
```

### Step 4: Differential Microexon Analysis in R
Filter alternative splicing events (`COMPLEX` $\in \{\text{S, C1, C2, C3, MIC, ANN}\}$), enforce non-zero variance and minimum non-`NA` sample counts per condition, compute two-sided Wilcoxon rank-sum $p$-values with Benjamini–Hochberg correction, and stratify significant events ($|\Delta\text{PSI}| > 10$, $p < 0.05$) across exon length bins (`1–15 nt`, `16–27 nt`, `28–51 nt`, and `≤ 51 nt`). Supply a two-column tab-separated metadata file (`sample` and `condition`) matching the sample column headers in `INCLUSION_LEVELS_FULL*.tab`:

```bash
Rscript scripts/04_differential_microexons.R \
  --input vast_align/INCLUSION_LEVELS_FULL-Hsa.tab \
  --metadata sample_metadata.tsv \
  --baseline Control \
  --case Disease \
  --name cohort_wilcox \
  --cutoff_baseline 3 \
  --cutoff_case 20 \
  --outdir results/differential_splicing
```
---

## Citation

If you use this pipeline in your work, please cite:

> In silico investigation of alternative splicing of microexons in human peripheral tissues
Shalini Raman, Pallavi Gupta, Ishaan Gupta
bioRxiv 2026.05.13.725022; doi: https://doi.org/10.64898/2026.05.13.725022
