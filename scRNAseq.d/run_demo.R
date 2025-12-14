# Script Name: run_demo.R
# Description: Minimal script to validate the computational environment for Liu et al.
#              It loads the 3 demo datasets and generates basic UMAP visualization.
# Expected Run time: < 2 minutes on a standard desktop.

message("========================================================")
message("   Starting Demo Run for Nature Communications Manuscript")
message("========================================================")

# 1. Load Required Libraries
# --------------------------------------------------------
message("\n[Step 1] Loading libraries...")
suppressPackageStartupMessages({
  if(!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is not installed!")
  if(!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is not installed!")
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})
message("   -> Libraries loaded successfully.")

# 2. Define Input Files 
# --------------------------------------------------------
demo_files <- list(
  "Ctrl_Subset"      = "demo_data_Ctrlonly.rds",
  "Integrated_Subset" = "demo_data_AllSample.rds",
  "SAP_Lineage"      = "demo_data_SAP.rds"
)

# Create output directory
output_dir <- "./demo_output/"
if(!dir.exists(output_dir)) dir.create(output_dir)

# 3. Processing Loop
# --------------------------------------------------------
message("\n[Step 2] Processing demo datasets...")

for (sample_name in names(demo_files)) {
  
  file_name <- demo_files[[sample_name]]
  
  # Check if file exists
  if(!file.exists(file_name)){
    warning(paste0("   [WARNING] File not found: ", file_name, ". Skipping..."))
    next
  }
  
  message(paste0("   -> Processing: ", sample_name, " (", file_name, ")"))
  
  # Read RDS
  obj <- readRDS(file_name)
  
  # Check if UMAP exists, if not, perform basic processing to generate it
  # (This proves that the computational linear algebra libraries are working)
  if(!"umap" %in% names(obj@reductions)){
    message("      - UMAP not found. Calculating basic UMAP...")
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:10, verbose = FALSE)
  }
  
  # Generate Plot
  p <- DimPlot(obj, reduction = "umap") + 
       ggtitle(paste0("Demo Plot: ", sample_name)) +
       theme_minimal() +
       theme(plot.title = element_text(hjust = 0.5))
  
  # Save Plot
  save_path <- paste0(output_dir, "Demo_UMAP_", sample_name, ".pdf")
  ggsave(save_path, plot = p, width = 6, height = 5)
  message(paste0("      - Plot saved to: ", save_path))
}

# 4. Final Success Message
# --------------------------------------------------------
message("\n========================================================")
message("   Demo Completed Successfully!")
message(paste0("   Please check the '", output_dir, "' folder for results."))
message("========================================================")