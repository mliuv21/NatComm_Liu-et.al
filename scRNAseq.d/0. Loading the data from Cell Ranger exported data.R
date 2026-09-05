# Project Name: Modelling NB by using patient-derived PSCs
# Script Name: 0. Loading the data from Cell Ranger exported data
# Name: Mingzhi Liu
# Date: 2023-11-11
# (Modified for GitHub Submission)

#### Loading the data from 10x ####
# rm(list = ls()) # Commented out for pipeline execution

if(TRUE){
  require(Seurat)
  require(dplyr)
  require(tidyverse)
}

# Define output directory
output_dir <- "./1.UncombinedObjects/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

### Loading data from NCC sample ###
# Note: Ensure raw data is placed in './data_raw/' as per README
data_dir_NCC <- './data_raw/NCC_feature_bc_matrix/'

if(!dir.exists(data_dir_NCC)){
  stop(paste("Data directory not found:", data_dir_NCC, "\nPlease download data from GEO and organize folder structure."))
}

list.files(data_dir_NCC) 
NCC_data <- Read10X(data.dir = data_dir_NCC)
seurat_object_NCC = CreateSeuratObject(counts = NCC_data$`Gene Expression`,
                                       project = "NCC")
seurat_object_NCC[["HTO"]] = CreateAssayObject(counts = NCC_data$`Antibody Capture`)

# Fixed typo: seruat -> seurat
saveRDS(seurat_object_NCC, paste0(output_dir, "seurat_object_NCC.RDS"))
message("NCC Object Saved.")


### Loading data from SAP sample ###
data_dir_SAP <- './data_raw/SAP_feature_bc_matrix/'
if(dir.exists(data_dir_SAP)){
  SAP_data <- Read10X(data.dir = data_dir_SAP)
  seurat_object_SAP = CreateSeuratObject(counts = SAP_data$`Gene Expression`, 
                                         project = "SAP")
  seurat_object_SAP[['HTO']] = CreateAssayObject(counts = SAP_data$`Antibody Capture`)
  saveRDS(seurat_object_SAP, paste0(output_dir, "seurat_object_SAP.RDS"))
  message("SAP Object Saved.")
}


### Loading data from D8 sample ###
data_dir_D8 <- './data_raw/D8_feature_bc_matrix/'
if(dir.exists(data_dir_D8)){
  D8_data <- Read10X(data.dir = data_dir_D8)
  seurat_object_D8 = CreateSeuratObject(counts = D8_data$`Gene Expression`, 
                                        project = "D8")
  seurat_object_D8[['HTO']] = CreateAssayObject(counts = D8_data$`Antibody Capture`)
  saveRDS(seurat_object_D8, paste0(output_dir, "seurat_object_D8.RDS"))
  message("D8 Object Saved.")
}


### Loading data from D12 sample ###
data_dir_D12 <- './data_raw/D12_feature_bc_matrix/'
if(dir.exists(data_dir_D12)){
  D12_data <- Read10X(data.dir = data_dir_D12)
  D12_data_GE <- D12_data$`Gene Expression`
  D12_data_HTO <- D12_data$`Antibody Capture`
  
  # Clean potential empty rows/cols
  D12_data_GE <- D12_data_GE[, colSums(D12_data_GE) > 0]
  D12_data_HTO <- D12_data_HTO[, colnames(D12_data_GE)]
  
  seurat_object_D12 = CreateSeuratObject(counts = D12_data_GE, 
                                         project = "D12")
  seurat_object_D12[['HTO']] = CreateAssayObject(counts = D12_data_HTO)
  saveRDS(seurat_object_D12, paste0(output_dir, "seurat_object_D12.RDS"))
  message("D12 Object Saved.")
}
