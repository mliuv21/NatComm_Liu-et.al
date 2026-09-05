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


#### Integration (SCTransform-based standard workflow) ####
obj_list <- list(Filter_seurat_NCC, Filter_seurat_SAP, Filter_seurat_D8, Filter_seurat_D12)
obj_list <- lapply(obj_list, function(x) {
  x <- CellCycleScoring(x, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  x <- SCTransform(x, vars.to.regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"), verbose = FALSE)
  x
})

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", anchor.features = features)
ML_seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Post-integration Processing
DefaultAssay(ML_seurat_integrated) <- "integrated"
ML_seurat_integrated <- RunPCA(ML_seurat_integrated, npcs = 50, verbose = FALSE)
ML_seurat_integrated <- RunUMAP(ML_seurat_integrated, reduction = "pca", dims = 1:22, verbose = FALSE)
ML_seurat_integrated <- FindNeighbors(ML_seurat_integrated, dims = 1:22)
ML_seurat_integrated <- FindClusters(ML_seurat_integrated, resolution = 0.2)

# Save Plot
p_final <- DimPlot(ML_seurat_integrated, reduction = "umap", group.by = "orig.ident")
ggsave(paste0(plot.path, "3.UMAP_after_integration.png"), p_final, width = 25, height = 20, units = "cm")

# Save Integrated Object for SCP.R
saveRDS(ML_seurat_integrated, paste0(path, "Integrated_Seurat_Object.RDS"))
message("Integration Complete. Object saved.")
