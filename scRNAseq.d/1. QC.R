'''
Project Name: Modelling NB by using patient-derived PSCs
Script Name: 1. QC.R
Description: Quality Control using adaptive thresholds (MAD-based) and HTO demultiplexing.
Name: Mingzhi Liu
Date: 2023-11-11
(Optimized for Publication/Reproduction)
'''

# 1. Load Libraries
# ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(ggExtra)
})

# 2. Setup Directories and Variables
# ------------------------------------------------------------------
samples <- c("NCC", "SAP", "D8", "D12")
input_dir <- "./1.UncombinedObjects/"
output_dir <- "./1.UncombinedObjects/" 
plot_root <- "./plots.d/1.QC.d/"

# Ensure output directory exists
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if(!dir.exists(plot_root)) dir.create(plot_root, recursive = TRUE)

message(">>> Starting QC Pipeline with Adaptive Thresholds (MAD)...")

# 3. Main Processing Loop
# ------------------------------------------------------------------
for (name in samples) {
  
  message(paste0("\nProcessing sample: ", name))
  
  # --- 3.1 Load Data ---
  file_path <- paste0(input_dir, "seurat_object_", name, ".RDS")
  if(!file.exists(file_path)) {
    warning(paste0("File not found: ", file_path, ". Skipping.")); next
  }
  
  data <- readRDS(file_path)
  data[["Barcode"]] <- rownames(data@meta.data) # Preserve original barcodes
  
  # Setup plot directory for this sample
  sample_plot_path <- paste0(plot_root, name, ".d/")
  if (!dir.exists(sample_plot_path)) dir.create(sample_plot_path, recursive = TRUE)
  
  # --- 3.2 Calculate QC Metrics ---
  # Mitochondrial genes (pattern: MT-)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  # Ribosomal genes (pattern: RPS or RPL)
  data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")
  
  # Initial Violin Plot
  feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
  p_vln <- VlnPlot(data, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 4) + NoLegend()
  ggsave(paste0(sample_plot_path, "QC_1_Vln_BeforeFilter_", name, ".png"), p_vln, width = 30, height = 20, units = "cm")
  
  
  # --- 3.3 Adaptive Filtering (MAD-based) ---
  meta_data <- data@meta.data
  
  # A. Mitochondrial Filter Limits
  max.mito.thr <- median(meta_data$percent.mt) + 5*mad(meta_data$percent.mt)
  min.mito.thr <- median(meta_data$percent.mt) - 5*mad(meta_data$percent.mt)
  
  # Plot: Mito vs nFeature
  p_mito <- ggplot(meta_data, aes(x=nFeature_RNA, y=percent.mt)) +
    geom_point(alpha=0.6) +
    geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
    annotate("text", x = mean(range(meta_data$nFeature_RNA)), y = max(meta_data$percent.mt), 
             label = paste0("Max Threshold: ", round(max.mito.thr, 2), "%")) +
    theme_bw()
  
  ggsave(paste0(sample_plot_path, "QC_2_Mito_nFeature_Histo_Scatter_", name, ".png"),
         ggMarginal(p_mito, type = "histogram", fill="lightgrey", bins=100), width = 20, height = 20, units = "cm")
  
  # B. Ribosome Filter Limits
  max.ribo.thr <- median(meta_data$percent.ribo) + 5*mad(meta_data$percent.ribo)
  min.ribo.thr <- median(meta_data$percent.ribo) - 5*mad(meta_data$percent.ribo)
  
  # Plot: Ribo vs nFeature
  p_ribo <- ggplot(meta_data, aes(x=nFeature_RNA, y=percent.ribo)) +
    geom_point(alpha=0.6) +
    geom_hline(aes(yintercept = max.ribo.thr), colour = "red", linetype = 2) +
    theme_bw()
  
  ggsave(paste0(sample_plot_path, "QC_3_Ribo_nFeature_Histo_Scatter_", name, ".png"),
         ggMarginal(p_ribo, type = "histogram", fill="lightgrey", bins=100), width = 20, height = 20, units = "cm")
  
  # C. Apply Mito/Ribo Filters
  meta_data <- meta_data %>% 
    filter(percent.mt < max.mito.thr & percent.mt > min.mito.thr) %>% 
    filter(percent.ribo < max.ribo.thr & percent.ribo > min.ribo.thr)
  
  
  # --- 3.4 nFeature & nCount Outlier Detection (Linear Regression) ---
  
  # Calculate thresholds for Genes (nFeature) and Counts (nCount)
  min.Genes.thr <- median(log10(meta_data$nFeature_RNA)) - 5*mad(log10(meta_data$nFeature_RNA))
  max.Genes.thr <- median(log10(meta_data$nFeature_RNA)) + 5*mad(log10(meta_data$nFeature_RNA))
  max.nCount.thr <- median(log10(meta_data$nCount_RNA)) + 5*mad(log10(meta_data$nCount_RNA))
  
  # Pre-filter for regression
  meta_data <- meta_data %>% 
    filter(log10(nFeature_RNA) > min.Genes.thr) %>% 
    filter(log10(nCount_RNA) < max.nCount.thr)
  
  # Build Linear Model (log-log relationship)
  lm.model <- lm(data = meta_data, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
  
  # Define valid cells (Residual cutoff: -0.09)
  # Logic: Cells falling significantly below the regression line are low quality
  meta_data$valideCells <- log10(meta_data$nFeature_RNA) > (log10(meta_data$nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - 0.09))
  
  # Plot: Regression and Cutoff
  p_reg <- ggplot(meta_data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point(aes(colour = valideCells), alpha=0.5) +
    geom_smooth(method="lm", color="blue") +
    geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange", linetype="dashed") + 
    labs(title = paste0("QC Model - ", name), subtitle = "Orange line: Lower bound threshold") +
    theme_bw() + theme(legend.position="bottom")
  
  ggsave(paste0(sample_plot_path, "QC_7_Regression_Filter_", name, ".png"),
         ggMarginal(p_reg, type = "histogram", fill="lightgrey"), width = 20, height = 20, units = "cm")
  
  
  # --- 3.5 Apply Final Filtering ---
  meta_data <- meta_data %>% filter(valideCells)
  data_filtered <- subset(data, subset = Barcode %in% meta_data$Barcode)
  
  message(paste0("   -> Cells before QC: ", ncol(data)))
  message(paste0("   -> Cells after QC:  ", ncol(data_filtered)))
  
  
  # --- 3.6 HTO Demultiplexing (Doublet Removal) ---
  # Note:  HTO data exists in the object
  if("HTO" %in% names(data_filtered@assays)) {
    data_filtered <- NormalizeData(data_filtered, assay = "HTO", normalization.method = "CLR", verbose = FALSE)
    data_filtered <- HTODemux(data_filtered, assay = "HTO", positive.quantile = 0.99)
    
    # Save HTO Ridge Plot
    ggsave(paste0(sample_plot_path, "QC_10_HTO_RidgePlot_", name, ".png"),
           RidgePlot(data_filtered, assay = "HTO", features = rownames(data_filtered[["HTO"]])[1:4], ncol = 4),
           width = 40, height = 20, units = "cm")
    
    # Subset to Singlets
    data_singlet <- subset(data_filtered, idents = "Singlet")
    message(paste0("   -> Singlets found:  ", ncol(data_singlet)))
    
  } else {
    warning("HTO assay not found! Skipping doublet removal.")
    data_singlet <- data_filtered
  }
  
  
  # --- 3.7 Save Final Object ---
  # Naming convention matches '2. Integration.R' requirement
  save_path <- paste0(output_dir, "Filter_Seurat_Object_", name, "singlet.RDS")
  saveRDS(data_singlet, file = save_path)
  message(paste0("   -> Saved: ", save_path))
}

message("\nQC Pipeline Finished Successfully.")