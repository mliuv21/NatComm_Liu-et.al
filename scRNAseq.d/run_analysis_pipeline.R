# Script Name: run_analysis_pipeline.R
# Description: Master script to execute the analysis pipeline (Steps 0-2).
#              Step 3 & 4 (Figure Generation) are provided as HTML reports.
# Usage: Run this script from the project root directory.

message("Starting Analysis Pipeline for Liu et al., Nature Communications...")

# Step 0: Load Raw Data
message("\n>>> Step 0: Loading Raw Data...")
if(!file.exists("./data_raw/NCC_feature_bc_matrix/matrix.mtx.gz")) {
  stop("Error: Raw data not found. Please see README for download instructions.")
}
source("0. Loading the data from Cell Ranger exported data.R")
message("Step 0 Complete. Raw objects saved.")

# Step 1: Quality Control
message("\n>>> Step 1: Running Quality Control (QC)...")
source("1. QC.R")
message("Step 1 Complete. Filtered objects saved.")

# Step 2: Integration
message("\n>>> Step 2: Performing Data Integration...")
source("2.Integration.R")
message("Step 2 Complete. Integrated object saved as './output/2.Integration.d/Integrated_Seurat_Object.RDS'")

# Step 3 & 4: Visualization Guidance
message("\n>>> Step 3 & 4: Figure Generation")
message("----------------------------------------------------------------")
message("NOTE: The code for generating the Main and Supplementary Figures")
message("is provided in the following standalone HTML reports:")
message("   - 3.SourceCode_Main_Figures.html")
message("   - 4.SourceCode_Supplementary_Figures.html")
message("----------------------------------------------------------------")
message("Please open these files in a web browser to view the source code")
message("and corresponding outputs used in the manuscript.")

message("\nPipeline execution finished.")