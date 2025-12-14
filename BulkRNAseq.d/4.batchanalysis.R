##########################################################################
#                                                                        #
#                  !Batch Analysis!                                      #
#                    Mingzhi Liu                                         #
#                   2022 - 03 - 02                                       #
#                                                                        #
#                                                                        #
##########################################################################

if (TRUE) {
  # Remove everything in the environment
  rm(list = ls())
  
  
  # Load essential packages
  library(DESeq2)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
}

################################# Batch analyze #########################


## INPUT FACTOR
dds.diff_group <- readRDS("2.ddsObj_diff-group.Rds")
results.NBvCtrl.10 <- results(dds.diff_group, 
                              contrast = list(c("Group_NB_vs_Ctrl", "DiffD10.GroupNB")))

sum(results.NBvCtrl.10$padj < 0.05, na.rm = TRUE)
table(results.NBvCtrl.10$padj < 0.05 & results.NBvCtrl.10$log2FoldChange < -1.5)

head(results.NBvCtrl.10)

analysis <- function(dds, results, compare_name){
  ## Set working directory
  ### Load packages
  library(DESeq2)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(ashr)
  
  ah <- AnnotationHub()
  hmEnsem <- ah[["AH109336"]]
  annotations <- genes(hmEnsem, return.type = "data.frame")
  
  annot <- as.data.frame(results) %>%
    rownames_to_column("gene_id") %>% 
    left_join(annotations, "gene_id") %>%
    dplyr::rename(logFC = log2FoldChange, FDR = padj,GeneID = gene_id,
                  Entrez = entrezid, Symbol = gene_name)
  
  write_tsv(annot, paste0("./4.batchanalysis.d/", compare_name, "_Results_Annotated.tsv"))
  
  print("Data loading successful, going to do analysis")
  
  ### Visualisation
  ## P-value histogram
  
  par(mfrow = c(1,1))
  
  pdf(paste0("./4.batchanalysis.d/plot.d/P-value", "_", compare_name, ".pdf"))
  hist(annot$pvalue)
  dev.off()
  
  print(" P value histogram have been saved!")
  ## Shrink the log2FC
  
  ddsShrink <- lfcShrink(dds,
                         res = results,
                         type = "ashr")
  
  shrinkTab <- as.data.frame(ddsShrink) %>%
    rownames_to_column("gene_id") %>% 
    left_join(annotations, "gene_id") %>%
    dplyr::rename(logFC = log2FoldChange, FDR = padj,GeneID = gene_id,
                  Entrez = entrezid, Symbol = gene_name)
  resultTab <- as.data.frame(results) %>%
    rownames_to_column("gene_id") %>% 
    left_join(annotations, "gene_id") %>%
    dplyr::rename(logFC = log2FoldChange, FDR = padj,GeneID = gene_id,
                  Entrez = entrezid, Symbol = gene_name)
  ## MA plots
  pdf(paste0("./4.batchanalysis.d/plot.d/MAplot_muilt", "_", compare_name, ".pdf"))
  par(mfrow = c(1,2))
  plotMA(results)
  plotMA(ddsShrink)
  dev.off()
  
  print("MA plots has been saved")
  
  ## Volcano plot
  volcanoTab <- shrinkTab %>%
    dplyr::mutate("-log10(pvalue)" = - log10(pvalue))
  volcanoTab$change = as.factor(
    ifelse(volcanoTab$FDR < 0.05 & abs(volcanoTab$logFC) >= 1, 
           ifelse(volcanoTab$logFC> 1 ,'Up','Down'),'NoSignifi'))
  
  print("Data to draw volcano plot has been organized, next to draw valcano plot in shrink")
  
  volcanoPlot <- ggplot(volcanoTab, aes(x = logFC, y =  `-log10(pvalue)`,colour = change, label = Symbol)) + 
    geom_label_repel(data = ~top_n(.x, 100, wt=-FDR), size = 5)+ 
    geom_point(size = 5) +
    geom_label_repel(data = volcanoTab %>% 
    dplyr::filter(Symbol %in% c("ALK","CDKN1A")), aes(label = Symbol)) + 
    scale_color_manual(values=c(Down = "#313695", 
                                NoSignifi = "grey",
                                Up = "#A50026")) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",linewidth=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",linewidth=0.8) +
    labs(x="log2(Fold Change)",y="-log10 (p-value)",title= compare_name) +  
    theme(panel.background = element_blank(),
          panel.grid = element_line(color="grey"),
          plot.title = element_text(face = 'bold',size =25, hjust = 0.5),
          axis.ticks.x=element_blank(), 
          axis.title.y = element_text(face = 'bold',size =25),
          axis.text.y  = element_text(colour="black", size=20, face="bold"),
          axis.title.x = element_text(face = 'bold',size =25),
          axis.text.x  = element_text(colour="black", size=20, face="bold", hjust = 1),
          legend.text  = element_text(colour="black", size=20, face="bold"),
          legend.text.align = 0,
          legend.title = element_text(colour="black", size=20, face="bold"),
          # legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
          legend.title.align = 1,
          legend.position = "bottom") + 
    guides(color = guide_legend(title = "Change", title.hjust = 0))
  volcanoPlot
  ggsave(paste0("./4.batchanalysis.d/plot.d/VolcanoPlot_shrink_", compare_name, ".pdf"), 
         dpi = 600, width = 20, height = 20, units = "in")
  print("Volcano Plot in shrink has been saved!")
  
  
  volcanoTab <- resultTab %>%
    dplyr::mutate("-log10(pvalue)" = - log10(pvalue))
  volcanoTab$change = as.factor(
    ifelse(volcanoTab$pvalue < 0.05 & abs(volcanoTab$logFC) >= 1, 
           ifelse(volcanoTab$logFC> 1 ,'Up','Down'),'NoSignifi'))
  
  print("Data to draw volcano plot has been organized, next to draw valcano plot")
  
  volcanoPlot <- ggplot(volcanoTab, aes(x = logFC, y =  `-log10(pvalue)`,colour = change, label = Symbol)) + 
    geom_label_repel(data = ~top_n(.x, 100, wt=-FDR), size = 5)+ 
    geom_point(size = 5) +
    geom_label_repel(data = volcanoTab%>% 
    dplyr::filter(Symbol %in% c("ALK", "CDKN1A")), aes(label = Symbol)) + 
    scale_color_manual(values=c(Down = "#313695", 
                                NoSignifi = "grey",
                                Up = "#A50026")) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",linewidth=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",linewidth=0.8) +
    labs(x="log2(Fold Change)",y="-log10 (p-value)",title= compare_name) +  
    theme(panel.background = element_blank(),
          panel.grid = element_line(color="grey"),
          plot.title = element_text(face = 'bold',size =25, hjust = 0.5),
          axis.ticks.x=element_blank(), 
          axis.title.y = element_text(face = 'bold',size =25),
          axis.text.y  = element_text(colour="black", size=20, face="bold"),
          axis.title.x = element_text(face = 'bold',size =25),
          axis.text.x  = element_text(colour="black", size=20, face="bold", hjust = 1),
          legend.text  = element_text(colour="black", size=20, face="bold"),
          legend.text.align = 0,
          legend.title = element_text(colour="black", size=20, face="bold"),
          # legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
          legend.title.align = 1,
          legend.position = "bottom") + 
    guides(color = guide_legend(title = "Change", title.hjust = 0))
  ggsave(paste0("./4.batchanalysis.d/plot.d/VolcanoPlot_", compare_name, ".pdf"), dpi = 600, width = 20, height = 20, units = "in")
  print("Volcano Plot in shrink has been saved!")
  
  saveRDS(results, file = paste0("./4.batchanalysis.d/ResultsObject_", compare_name, ".Rds"))
  saveRDS(ddsShrink, file = paste0("./4.batchanalysis.d/ddshrinkObject_", compare_name, ".Rds"))
  saveRDS(shrinkTab, file = paste0("./4.batchanalysis.d/ddshrinkResult_", compare_name, ".Rds"))
  saveRDS(annot, file = paste0("./4.batchanalysis.d/Results_", compare_name, "_annotation.Rds"))
  
  print("DONE!! CONGRATULATIONS!!")
  
}


results.NBvCtrl.12 <- results(dds.diff_group, 
                              contrast = list(c("Group_NB_vs_Ctrl", "DiffD12.GroupNB")))


analysis(dds = dds.diff_group,
         results = results.NBvCtrl.12,
         compare_name = "NBvsCtrl_12")

results.NBvCtrl.10 <- results(dds.diff_group, 
                              contrast = list(c("Group_NB_vs_Ctrl", "DiffD10.GroupNB")))


analysis(dds = dds.diff_group,
         results = results.NBvCtrl.10,
         compare_name = "NBvsCtrl_10")




results.NBvCtrl.8 <- results(dds.diff_group, 
                             contrast = list(c("Group_NB_vs_Ctrl", "DiffD8.GroupNB")))


analysis(dds = dds.diff_group,
         results = results.NBvCtrl.8,
         compare_name = "NBvsCtrl_8")



results.NBvCtrl.6 <- results(dds.diff_group, 
                             contrast = list(c("Group_NB_vs_Ctrl", "DiffD6.GroupNB")))


analysis(dds = dds.diff_group,
         results = results.NBvCtrl.6,
         compare_name = "NBvsCtrl_6")


results.NBvCtrl.SAP <- results(dds.diff_group, 
                               contrast = list(c("Group_NB_vs_Ctrl", "DiffSAP.GroupNB")))


analysis(dds = dds.diff_group,
         results = results.NBvCtrl.SAP,
         compare_name = "NBvsCtrl_SAP")


results.NBvCtrl.NCC <- results(dds.diff_group, 
                               contrast = list(c("Group_NB_vs_Ctrl")))


analysis(dds = dds.diff_group,
         results = results.NBvCtrl.NCC,
         compare_name = "NBvsCtrl_NCC")

