#Project Name: Modelling NB by using patient-derived PSCs
#Script Name: 2. Integrate Seurat Objects
#Name: Mingzhi Liu
#Date: 2023-11-11
#(Modified for GitHub Submission)

#### Loading the data from 10x ####
# rm(list = ls()) # Commented out

if(TRUE){
  require(Seurat)
  require(dplyr)
  require(tidyverse)
  require(ggExtra)
  require(RColorBrewer)
  
  path <- "./output/2.Integration.d/"
  plot.path <- "./plots.d/2.Integration.d/"  
  
  if(!dir.exists(plot.path)) dir.create(plot.path, recursive = TRUE)
  if(!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# Load Filtered Objects (using relative paths)
input_dir <- "./1.UncombinedObjects/"

# Check if files exist
files_to_read <- c(
  paste0(input_dir, "Filter_Seurat_Object_NCCsinglet.RDS"),
  paste0(input_dir, "Filter_Seurat_Object_SAPsinglet.RDS"),
  paste0(input_dir, "Filter_Seurat_Object_D8singlet.RDS"),
  paste0(input_dir, "Filter_Seurat_Object_D12singlet.RDS")
)

if(!all(file.exists(files_to_read))){
   stop("Some input files from Step 1 are missing. Please check ./1.UncombinedObjects/")
}

Filter_seurat_NCC <- readRDS(files_to_read[1])
Filter_seurat_SAP <- readRDS(files_to_read[2])
Filter_seurat_D8  <- readRDS(files_to_read[3])
Filter_seurat_D12 <- readRDS(files_to_read[4])

# Merge
ML_data <- merge(Filter_seurat_NCC, 
                 y = c(Filter_seurat_SAP, Filter_seurat_D8, Filter_seurat_D12),
                 add.cell.ids = c("NCC", "SAP", "D8", "D12"), 
                 project = "NB_Project")

# Basic Processing
ML_seurat <- NormalizeData(ML_data)
ML_seurat <- FindVariableFeatures(ML_seurat)
ML_seurat <- ScaleData(ML_seurat)
ML_seurat <- RunPCA(ML_seurat)

# Pre-integration Plot
DimPlot(ML_seurat, reduction = "pca", group.by = "orig.ident")
ggsave(paste0(plot.path, "2.UMAP_before_integration_orig.png"), width = 20, height = 15, units = "cm")


#### Integration (RPCA) ####
features <- SelectIntegrationFeatures(object.list = list(Filter_seurat_NCC, Filter_seurat_SAP, Filter_seurat_D8, Filter_seurat_D12))

# Prepare list for integration
obj_list <- list(Filter_seurat_NCC, Filter_seurat_SAP, Filter_seurat_D8, Filter_seurat_D12)
obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
  return(x)
})

anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, reduction = "rpca")
ML_seurat_integrated <- IntegrateData(anchorset = anchors)

# Post-integration Processing
DefaultAssay(ML_seurat_integrated) <- "integrated"
ML_seurat_integrated <- ScaleData(ML_seurat_integrated, verbose = F)
ML_seurat_integrated <- RunPCA(ML_seurat_integrated, npcs = 30, verbose = F)
ML_seurat_integrated <- RunUMAP(ML_seurat_integrated, reduction = "pca", dims = 1:30, verbose = F)

# Save Plot
p_final <- DimPlot(ML_seurat_integrated, reduction = "umap", group.by = "orig.ident")
ggsave(paste0(plot.path, "3.UMAP_after_integration.png"), p_final, width = 25, height = 20, units = "cm")

# Save Integrated Object for SCP.R
saveRDS(ML_seurat_integrated, paste0(path, "Integrated_Seurat_Object.RDS"))
message("Integration Complete. Object saved.")