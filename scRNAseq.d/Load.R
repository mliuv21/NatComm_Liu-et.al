if(T){
  require(Seurat)
  require(dplyr)
  require(tidyverse)
  require(ggExtra)
  require(RColorBrewer)
  require(cluster)
  require(cowplot)
  require(clustree)
  require(viridis)
  require(scCustomize)
  require(paletteer)
  library(slingshot)
  library(SingleCellExperiment)
  library(rafalib)
  require(patchwork)
  require(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(SingleR)
  require(reshape2)
  library(ggpubr)
  library(clusterProfiler)
  #remotes::install_github("pwwang/scplotter")
  # If needed for SCP, replace this with your local conda binary path.
  # options(reticulate.conda_binary = "/path/to/conda", SCP_env_name = "SCP_env")
  
  library(SCP)
  # Replace this with your working directory containing the intermediate objects below.
  # setwd("your_working_directory")

col <-  DiscretePalette_scCustomize(num_colors = 70, 
                                    palette = "varibow", 
                                    shuffle_pal = TRUE,
                                    seed = 42)

cluster_col <- DiscretePalette_scCustomize(num_colors = 14, 
                                           palette = "varibow", 
                                           shuffle_pal = TRUE,
                                           seed = 22)


ident_col <- c("#3D98D3FF", "#9C27B0FF", "#8F7289FF", "#4CAF50FF", "#F49600FF", "#C62828FF")
ident_col <- ident_col[c(1,2,4,6)]
sample_col <- c("#29B6F6FF" ,"#96281BFF", "#FF73B3FF" , "#891FCCFF"  )
sampleinfo_col <-  c("#4575B4", "#D73027")
label_transfer_col <- DiscretePalette_scCustomize(num_colors = 14, 
                                                  palette = "varibow", 
                                                  shuffle_pal = TRUE,
                                                  seed = 42)
inte_color <- DiscretePalette_scCustomize(num_colors = 20,
                                          palette = "alphabet2", 
                                          shuffle_pal = TRUE)
SAP_cluster_col <- c('#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FDBF6F','#FF7F00',
                     '#FB9A99','#E31A1C','#CAB2D6','#6A3D9A','#FFFF99','#B15928')
Cycle_col <- c("#66C2A5FF", "#FC8D62FF", "#8DA0CBFF")
}

# The following objects are intermediate outputs from the original analysis workflow.
ML_seurat_Ctrl <- readRDS("./Ctrl.d/3.Embedding.d/Name.d/ML_seurat_normalized_withSCT_Name_C10_pub.RDS")

ML_seurat_SAP <- readRDS("./All_Sample.d/ML_seurat_SAP_pub.RDS")
ML_seurat_all <- readRDS("./All_Sample.d/3.Embedding.d/ML_seurat_normalized_withSCT_NBClassification_pub.RDS")
ML_seurat_all@meta.data$SampleInfo <- factor(ML_seurat_all@meta.data$SampleInfo, levels = c("Ctrl", "NB"))
ML_seurat_all@meta.data$Individuals <- factor(ML_seurat_all@meta.data$Individuals, levels = c("Ctrl10", "NB1", "NB2", "NB4"))

ML_seurat_Ctrl@meta.data$cluster_ids_internal <- factor(ML_seurat_Ctrl@meta.data$cluster_ids_internal, 
                                          levels = c("NCC-1", "NCC-2\n(MT^high)", "SAP-1\n(Progenitors)", "SAP-2", "SAM-1\n(MES)",  "SAM-2\n(SYM)", 
                                                     "Lineage-1\n(SA)", "Lineage-2\n(Immature SA)", "Lineage-3\n(NE)", "Other"))
ML_seurat_Ctrl@meta.data$cluster_ids <- factor(ML_seurat_Ctrl@meta.data$cluster_ids, 
                                          levels = c("NCC-1", "NCC-2", "SAP-1", "SAP-2", "SAM-1",  "SAM-2", "Lineage-1", "Lineage-2", "Lineage-3", "Other"))


ALK_REACTOME_genes <- c(
  "PIK3CB", "PRDM1", "EP300", "HIF1A", "JAK3", "PIK3R2", "PTN", "PTPRZ1", "MDK",
  "PTPN6", "HDAC1", "CD274", "PIK3CA", "PLCG1", "DNMT1", "MYCN", "MYC",
  "PIK3R1", "IL2RG", "SHC1", "FRS2", "STAT3", "IRS1", "SIN3A", "ALK", "HDAC3",
  "ALKAL2", "HDAC2", "ALKAL1", "SRC"
)
ML_seurat_SAP <- AddModuleScore(ML_seurat_SAP, 
                                 features = list(ALK_REACTOME_genes), 
                                 name = "ALK_REACTOME", assay = "SCT")

nb_cells <- subset(ML_seurat_SAP, subset = SampleInfo == "NB")

# Add new colume score 

cluster5_markers <- read.csv("./Cluster5_Markers_all.csv", header = T, row.names =1)
cluster5_markers %>% head
# Set up the Marker gene for cluster 5
cluster5marker_score <- cluster5_markers[cluster5_markers$p_val_adj < 0.01 & cluster5_markers$avg_log2FC > 0.5,] %>% rownames() 
cluster5marker_score
