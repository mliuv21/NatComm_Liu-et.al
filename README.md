# Analysis Code for: ALK R1275Q mutation drives expansion of SCP-like cells during sympathoadrenal commitment and primes neuroblastoma initiation

This repository contains the complete source code and computational pipelines used for the analysis presented in the manuscript by **Liu et al.**.

The analysis is divided into two main parts: 1. **Single-cell RNA-seq (scRNA-seq)**: Seurat-based analysis pipeline. 2. **Bulk RNA-seq**: DESeq2-based differential expression analysis.

## 📁 Repository Structure

``` text
Liu et. al/
├── README.md               # Overview and instructions
├── LICENSE                 # MIT License
├── scRNA_seq.d/              # Single-cell analysis pipeline
│   ├── run_analysis_pipeline.R     # Master script for full analysis
│   ├── run_demo.R                  # Script for quick validation
│   ├── 0. Loading the data...R
│   ├── 1. QC.R
│   ├── 2. Integration.R
│   ├── 3.SourceCode_Main_Figures.html          # Code for Main Figures
│   ├── 4.SourceCode_Supplementary_Figures.html # Code for Supp Figures
│   ├── Load.R                      # Library dependencies
│   └── [Demo Data files: .rds]
│
└── Bulk_RNA_seq.d/           # Bulk RNA-seq analysis pipeline
├── 1.DataPreparation.R         # (Reference only)
├── 2.LinearModelnDESeq2.R      # Start here for reproduction
├── 3.annotation.R
├── 4.batchanalysis.R
├── 5.visualization.R
├── 6.enrichment.R
├── SampleInfo.txt              # Metadata
├── t2g_hm109_gencode.txt       # Gene annotation
└── 1.txi.count.rds             # Pre-processed count matrix
```

💻 System Requirements OS: Windows, Mac, or Linux.

R Version: 4.4.1

Key Dependencies:

scRNA-seq: Seurat (v4.4.0), SCP (v0.5.6), Monocle3, slinghot, tidyverse.

Bulk RNA-seq: DESeq2, tximport, clusterProfiler.

Note: For a full list of required packages, please refer to scRNA_seq.d/Load.R.

# 🛠️ Installation Clone this repository:

``` bash
git clone [https://github.com/](https://github.com/)[Your_Username]/[Your_Repo_Name].git
cd [Your_Repo_Name]
```

Install R Packages: Open R and install the required dependencies. You can check scRNA_seq.d/Load.R for details.

``` r
# Example installation commands
install.packages(c("Seurat", "tidyverse", "BiocManager", "devtools"))
devtools::install_github("zhanghao-njmu/SCP")
BiocManager::install("DESeq2")
```

Typical installation time: 30-60 minutes.

# 🧬 Part 1: Single-Cell RNA-seq Analysis Located in the ./scRNA_seq.d/ folder.

1.  Quick Demo (Validation) To verify the computational environment without downloading large raw datasets, we provide a Demo Suite using subsetted data.

``` r
setwd("./scRNA_seq.d")
source("run_demo.R")
```

Expected Output: Generates 3 UMAP plots in './demo_output/' verifying the installation. Run Time: \< 5 minutes.

2.  Full Pipeline Reproduction Step 0 (Download Data): Download raw 10x Genomics data from GEO (Accession: GSE310952). Place the data in a folder named data_raw inside the scRNA_seq.d directory (or modify the path in script 0.).

Run Pipeline:

``` r
setwd("./scRNA_seq.d")
source("run_analysis_pipeline.R")
# This executes Steps 0 (Loading), 1 (QC), and 2 (Integration).
```

3.  Figure Generation Refer to the HTML reports (3.SourceCode_Main_Figures.html and 4.SourceCode_Supplementary_Figures.html) for the exact R code chunks used to generate the manuscript figures.

# 🔬 Part 2: Bulk RNA-seq Analysis Located in the ./Bulk_RNA_seq.d/ folder.

Note on Reproducibility: Raw sequencing data (FASTQ/Salmon output) is not included due to size constraints. However, we provide the pre-processed count matrix (1.txi.count.rds) to allow full reproduction starting from the differential expression step. Processed count matrices have been deposited in GEO under accession number **GSE310836**. For code reproducibility within this repository, we provide the R-compatible `1.txi.count.rds` file so Step 2 can be executed immediately.

Instructions: Navigate to the folder:

``` r
setwd("../Bulk_RNA_seq.d")  # Or set absolute path to Bulk_RNA_seq
```

Start from Step 2 (Differential Expression): Script 1.DataPreparation.R is provided for methodological transparency but should be skipped. Start the analysis by running Step 2, which loads the provided 1.txi.count.rds and SampleInfo.txt.

``` r
source("2.LinearModelnDESeq2.R")
```

Note: This script calculates differential expression and saves the results for subsequent steps.

Run Downstream Analysis (Steps 3-6): Execute the remaining scripts sequentially to perform annotation, batch analysis, visualization, and enrichment.

``` r
source("3.annotation.R")
source("4.batchanalysis.R")
source("5.visualization.R")
source("6.enrichment.R")
```

📄 License This project is licensed under the MIT License - see the LICENSE file for details.

📧 Contact For any questions regarding the code or data, please contact the corresponding author as listed in the manuscript.
