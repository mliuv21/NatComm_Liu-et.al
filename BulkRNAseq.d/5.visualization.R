##########################################################################
#                                                                        #
#                  !Visualization!                                       #
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
  library(genefilter)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggVennDiagram)
}

dds.diff_group <- readRDS("2.ddsObj_diff-group.Rds")
txi <- readRDS("1.txi.count.rds")
rawCounts <- round(txi$counts, 0)
keep <- rowSums(rawCounts) > 5
rawCounts <- rawCounts[keep, ]


###### Heatmaps ####

### all genes
## Calculate vst for the whole dataset


gene <- vst(dds.diff_group) %>% assay()
gene <- rlog(dds.diff_group) %>% assay() #!#

topgenenum <- 1000
topVarGenes <- head(order(rowVars(gene),decreasing=TRUE),topgenenum)
scaled_plotDat = t(scale(t(gene[topVarGenes,]), center = T, scale = T))

# If you need show Symbol on your HM
# annt_scaled_plotDat <- annotations %>% 
#   dplyr::select(gene_id, gene_name) %>%
#   dplyr::filter(gene_id %in% rownames(scaled_plotDat))
# annt_scaled_plotDat <- as.data.frame(scaled_plotDat) %>%
#   rownames_to_column("gene_id") %>% 
#   left_join(annt_scaled_plotDat, "gene_id") %>%
#   rename(GeneID = gene_id, Symbol = gene_name) %>%
#   dplyr::select(!"GeneID") %>%
#   dplyr::filter(!is.na(Symbol)) 
# rownames(annt_scaled_plotDat) <- make.unique(annt_scaled_plotDat$Symbol)
# annt_scaled_plotDat %>% dplyr::select(!"Symbol") -> annt_scaled_plotDat
# dim(annt_scaled_plotDat)

#Color Palette
myPalette <- c("#4575B4", "#FFFFBF", "#D73027")
myRamp <- colorRamp2(c(-2,0, 2), myPalette)
brewer.pal(n = 7, name = "Greens")
colData(dds.diff_group)[,"Diff"]  <- factor(colData(dds.diff_group)[,"Diff"], levels = c("NCC","SAP", "D6", "D8", "D10", "D12"))
ha1 <- HeatmapAnnotation(df = colData(dds.diff_group)[,c("Diff", "Group")],
                         col = list(Group = c("Ctrl" = "#4575B4",
                                              "NB" = "#D73027"),
                                    Diff = c("NCC" = "#C7E9C0",
                                             "SAP" = "#A1D99B",
                                             "D6"  = "#74C476",
                                             "D8"  = "#41AB5D",
                                             "D10" = "#238B45",
                                             "D12" = "#005A32")),
                         "log2 ALK" = anno_barplot(as.numeric(log2(rawCounts[rownames(rawCounts) == "ENSG00000171094",] + 1))))
hmap_tpm <- Heatmap(scaled_plotDat, 
                name = "z-score",
                col = myRamp, 
                clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                column_names_side = "bottom", column_dend_side = "top",
                show_row_names = F,
                row_split = 8, column_split = 4,
                #row_title = c("Healthy related", "SAP related", "NB related", "NCC related"),
                column_title = NULL,
                border = TRUE, 
                top_annotation = ha1,
                column_dend_height = unit(5, "mm"),
                use_raster = TRUE, raster_quality = 1)

png("./5.visualization.d/Heatmap_top1kMVG_unsupervised_TPM_r8c4.png", width=10,height=8,units="in",res=2000)
draw(hmap_tpm)
dev.off()


HM = draw(hmap_tpm)
saveRDS(hmap_tpm, "./5.visualization_HM.Rds")
############### HM Gene Extraction #############################
# extract genes from clusters
r.dend <- row_dend(HM)  #Extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters

# loop to extract genes for each cluster.
for (i in 1:length(row_order(HM))){
  if (i == 1) {
    clu <- t(t(row.names(scaled_plotDat[row_order(HM)[[i]],])))
    gene_cluster <- cbind(clu, paste("cluster", i, sep=""))
    colnames(gene_cluster) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(scaled_plotDat[row_order(HM)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    gene_cluster <- rbind(gene_cluster, clu)
  }
}
#check
gene_cluster <- gene_cluster %>% 
  as.data.frame() %>%
  mutate(Cluster = str_replace(Cluster, "cluster1", "Healthy"),
         Cluster = str_replace(Cluster, "cluster2", "SAP"),
         Cluster = str_replace(Cluster, "cluster3", "NB"),
         Cluster = str_replace(Cluster, "cluster4", "NCC"))




NB_related <- gene_cluster %>% 
  dplyr::filter(Cluster == "NB")
SAP_related <- gene_cluster %>% 
  dplyr::filter(Cluster == "SAP")
Healthy_related <- gene_cluster %>% 
  dplyr::filter(Cluster == "Healthy")
NCC_related <- gene_cluster %>% 
  dplyr::filter(Cluster == "NCC")
HMgene <- list(NB_related, NCC_related, Healthy_related,SAP_related)

saveRDS(HMgene, "HMgene.Rds")




##### HM without NCC ######

sampleinfo <- read_tsv("./SampleInfo.txt", col_types = c("cfcfc"))
sampleinfo_noNCC <- sampleinfo %>% dplyr::filter(!grepl('NCC', SampleName))

coldata <- mutate(sampleinfo_noNCC, Diff = fct_relevel(Diff, "SAP"),
                  Group = fct_relevel(Group, "Ctrl"))


dds.diff_group_noNCC <- dds.diff_group[,colnames(dds.diff_group) %in% sampleinfo_noNCC$SampleName]
gene <- rlog(dds.diff_group_noNCC) %>% assay() #!#

topgenenum <- 1000
topVarGenes <- head(order(rowVars(gene),decreasing=TRUE),topgenenum)


scaled_plotDat = t(scale(t(gene[topVarGenes,]), center = T, scale = T))
myPalette <- c("#4575B4", "#FFFFBF", "#D73027")
myRamp <- colorRamp2(c(-2,0, 2), myPalette)
brewer.pal(n = 7, name = "Greens")
colData(dds.diff_group_noNCC)[,"Diff"]  <- factor(colData(dds.diff_group_noNCC)[,"Diff"], levels = c("NCC","SAP", "D6", "D8", "D10", "D12"))
rawCounts <- rawCounts %>% as.data.frame %>% select(!contains("NCC"))
ha1 <- HeatmapAnnotation(df = colData(dds.diff_group_noNCC)[,c("Diff", "Group")],
                         col = list(Group = c("Ctrl" = "#4575B4",
                                              "NB" = "#D73027"),
                                    Diff = c("SAP" = "#A1D99B",
                                             "D6"  = "#74C476",
                                             "D8"  = "#41AB5D",
                                             "D10" = "#238B45",
                                             "D12" = "#005A32")),
                         "log2 ALK" = anno_barplot(as.numeric(log2(rawCounts[rownames(rawCounts) == "ENSG00000171094",] + 1))))
hmap_tpm <- Heatmap(scaled_plotDat, 
                    name = "z-score",
                    col = myRamp, 
                    clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                    column_names_side = "bottom", column_dend_side = "top",
                    show_row_names = F,
                    row_split = 5, column_split = 3,
                    #row_title = c("Healthy related", "SAP related", "NB related", "NCC related"),
                    column_title = NULL,
                    border = TRUE, 
                    top_annotation = ha1,
                    column_dend_height = unit(5, "mm"),
                    use_raster = TRUE, raster_quality = 1)

png("./5.visualization.d/Heatmap_top1kMVG_unsupervised_TPM_noNCC_r5c3.png", width=10,height=8,units="in",res=2000)
draw(hmap_tpm)
dev.off()


############### HM Gene Extraction #############################
# extract genes from clusters
r.dend <- row_dend(HM)  #Extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters

# loop to extract genes for each cluster.
for (i in 1:length(row_order(HM))){
  if (i == 1) {
    clu <- t(t(row.names(scaled_plotDat[row_order(HM)[[i]],])))
    gene_cluster <- cbind(clu, paste("cluster", i, sep=""))
    colnames(gene_cluster) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(scaled_plotDat[row_order(HM)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    gene_cluster <- rbind(gene_cluster, clu)
  }
}
#check
gene_cluster <- gene_cluster %>% 
  as.data.frame() %>%
  mutate(Cluster = str_replace(Cluster, "cluster1", "Healthy"),
         Cluster = str_replace(Cluster, "cluster2", "SAP"),
         Cluster = str_replace(Cluster, "cluster3", "NB"),
         Cluster = str_replace(Cluster, "cluster4", "NCC"))




NB_related <- gene_cluster %>% 
  dplyr::filter(Cluster == "NB")
SAP_related <- gene_cluster %>% 
  dplyr::filter(Cluster == "SAP")
Healthy_related <- gene_cluster %>% 
  dplyr::filter(Cluster == "Healthy")
NCC_related <- gene_cluster %>% 
  dplyr::filter(Cluster == "NCC")
HMgene <- list(NB_related, NCC_related, Healthy_related,SAP_related)

saveRDS(HMgene, "HMgene.Rds")


############################ Venn plot ###############################
library(ggvenn)
results.NBvCtrl.12 <- readRDS("4.batchanalysis.d/ResultsObject_NBvsCtrl_12.Rds")
results.NBvCtrl.10 <- readRDS("4.batchanalysis.d/ResultsObject_NBvsCtrl_10.Rds")
results.NBvCtrl.8 <- readRDS("4.batchanalysis.d/ResultsObject_NBvsCtrl_8.Rds")
results.NBvCtrl.6 <- readRDS("4.batchanalysis.d/ResultsObject_NBvsCtrl_6.Rds")
results.NBvCtrl.SAP <- readRDS("4.batchanalysis.d/ResultsObject_NBvsCtrl_SAP.Rds")
results.NBvCtrl.NCC <- readRDS("4.batchanalysis.d/ResultsObject_NBvsCtrl_NCC.Rds")

annotations <- readRDS("4.batchanalysis.d/Results_NBvsCtrl_12_annotation.Rds")

vennDat.10_12 <- tibble(GeneID = rownames(results.NBvCtrl.12)) %>% 
  mutate(Upregualted_NB_12 = results.NBvCtrl.12$padj < 0.05 &
           !is.na(results.NBvCtrl.12$padj) & 
           results.NBvCtrl.12$log2FoldChange > 1.5) %>%
  mutate(Downregualted_NB_12 = results.NBvCtrl.12$padj < 0.05 &
           !is.na(results.NBvCtrl.12$padj) & 
           results.NBvCtrl.12$log2FoldChange < -1.5) %>%
  mutate(Upregualted_NB_10 = results.NBvCtrl.10$padj < 0.05 &
           !is.na(results.NBvCtrl.10$padj) & 
           results.NBvCtrl.10$log2FoldChange > 1.5) %>%
  mutate(Downregualted_NB_10 = results.NBvCtrl.10$padj < 0.05 &
           !is.na(results.NBvCtrl.10$padj) & 
           results.NBvCtrl.10$log2FoldChange < -1.5)
ggvenn(vennDat.10_12, set_name_size = 8,
       text_size = 8)
ggsave("./5.visualization.d/Vennplot_D10vsD12.png",dpi = 600, width = 30, height = 20, units = "in")

vennDat.8_10 <- tibble(GeneID = rownames(results.NBvCtrl.8)) %>% 
  mutate(Upregualted_NB_8 = results.NBvCtrl.8$padj < 0.05 &
           !is.na(results.NBvCtrl.8$padj) & 
           results.NBvCtrl.8$log2FoldChange > 1.5) %>%
  mutate(Downregualted_NB_8 = results.NBvCtrl.8$padj < 0.05 &
           !is.na(results.NBvCtrl.8$padj) & 
           results.NBvCtrl.8$log2FoldChange < -1.5) %>%
  mutate(Upregualted_NB_10 = results.NBvCtrl.10$padj < 0.05 &
           !is.na(results.NBvCtrl.10$padj) & 
           results.NBvCtrl.10$log2FoldChange > 1.5) %>%
  mutate(Downregualted_NB_10 = results.NBvCtrl.10$padj < 0.05 &
           !is.na(results.NBvCtrl.10$padj) & 
           results.NBvCtrl.10$log2FoldChange < -1.5)
ggvenn(vennDat.8_10, set_name_size = 8,
       text_size = 8)
ggsave("./5.visualization.d/Vennplot_D8vsD10.png",dpi = 600, width = 30, height = 20, units = "in")


vennDat.6_8 <- tibble(GeneID = rownames(results.NBvCtrl.8)) %>% 
  mutate(Upregualted_NB_8 = results.NBvCtrl.8$padj < 0.05 &
           !is.na(results.NBvCtrl.8$padj) & 
           results.NBvCtrl.8$log2FoldChange > 1.5) %>%
  mutate(Downregualted_NB_8 = results.NBvCtrl.8$padj < 0.05 &
           !is.na(results.NBvCtrl.8$padj) & 
           results.NBvCtrl.8$log2FoldChange < -1.5) %>%
  mutate(Upregualted_NB_6 = results.NBvCtrl.6$padj < 0.05 &
           !is.na(results.NBvCtrl.6$padj) & 
           results.NBvCtrl.6$log2FoldChange > 1.5) %>%
  mutate(Downregualted_NB_6 = results.NBvCtrl.6$padj < 0.05 &
           !is.na(results.NBvCtrl.6$padj) & 
           results.NBvCtrl.6$log2FoldChange < -1.5)
ggvenn(vennDat.6_8, set_name_size = 8,
       text_size = 8)
ggsave("./5.visualization.d/Vennplot_D6vsD8.png",dpi = 600, width = 30, height = 20, units = "in")


vennDat.SAP_6 <- tibble(GeneID = rownames(results.NBvCtrl.SAP)) %>% 
  mutate(Upregualted_NB_SAP = results.NBvCtrl.SAP$padj < 0.05 &
           !is.na(results.NBvCtrl.SAP$padj) & 
           results.NBvCtrl.SAP$log2FoldChange > 1.5) %>%
  mutate(Downregualted_NB_SAP = results.NBvCtrl.SAP$padj < 0.05 &
           !is.na(results.NBvCtrl.SAP$padj) & 
           results.NBvCtrl.SAP$log2FoldChange < -1.5) %>%
  mutate(Upregualted_NB_6 = results.NBvCtrl.6$padj < 0.05 &
           !is.na(results.NBvCtrl.6$padj) & 
           results.NBvCtrl.6$log2FoldChange > 1.5) %>%
  mutate(Downregualted_NB_6 = results.NBvCtrl.6$padj < 0.05 &
           !is.na(results.NBvCtrl.6$padj) & 
           results.NBvCtrl.6$log2FoldChange < -1.5)
ggvenn(vennDat.SAP_6, set_name_size = 8,
       text_size = 8)
ggsave("./5.visualization.d/Vennplot_SAPvsD6.png",dpi = 600, width = 30, height = 20, units = "in")


vennDat.up.SAP_D12 <- tibble(GeneID = rownames(results.NBvCtrl.12)) %>% 
  mutate(Upregulated_NB_SAP= results.NBvCtrl.SAP$padj < 0.05 &
           !is.na(results.NBvCtrl.SAP$padj) & 
           results.NBvCtrl.SAP$log2FoldChange > 1.5) %>%
  mutate(Upregulated_NB_D6 = results.NBvCtrl.6$padj < 0.05 &
           !is.na(results.NBvCtrl.6$padj) & 
           results.NBvCtrl.6$log2FoldChange > 1.5) %>%
  mutate(Upregulated_NB_D8 = results.NBvCtrl.8$padj < 0.05 &
           !is.na(results.NBvCtrl.8$padj) & 
           results.NBvCtrl.8$log2FoldChange > 1.5) %>%
  mutate(Upregulated_NB_D10 = results.NBvCtrl.10$padj < 0.05 &
           !is.na(results.NBvCtrl.10$padj) & 
           results.NBvCtrl.10$log2FoldChange > 1.5) %>% 
  mutate(Upregulated_NB_D12 = results.NBvCtrl.12$padj < 0.05 &
           !is.na(results.NBvCtrl.12$padj) & 
           results.NBvCtrl.12$log2FoldChange > 1.5) %>% column_to_rownames("GeneID")
ggvenn(vennDat.up.SAP_D12, set_name_size = 8,
       text_size = 8)
venn(vennDat.up.SAP_D12, 
     ilabels = "counts", zcolor = "style", box = F, size = 1.5, linetype ="dashed",
     sncs = 1, ilcs = 1, ggplot = T) #+
 # text(485,450, "PAX8 \nPPP1R11 \n MICA \n SERF1A",  cex = .5 )
ggsave("./5.visualization.d/Vennplot_SAP-D12_up_5.png",dpi = 600, width = 20, height = 20, units = "cm")

row_sums <- apply(vennDat.up.SAP_D12, 1, sum)

# Convert row_sums to a data frame if it's not already
row_sums_df <- data.frame(GeneID = rownames(vennDat.up.SAP_D12), row_sums)

# Merge the row_sums data with the annotations data frame by "GeneID"
merged_data <- merge(row_sums_df, annotations, by.x = "GeneID", by.y = "GeneID", all.x = TRUE)

# Extract the "Symbol" column along with row_sums
result <- merged_data[, c("GeneID", "row_sums", "Symbol")]

# If you are only interested in rows where the sum is 5
result[result$row_sums == 5, ]




ggsave("./5.visualization.d/Vennplot_SAP-D12_up.png",dpi = 600, width = 30, height = 20, units = "in")


vennDat.down.SAP_D12 <- tibble(GeneID = rownames(results.NBvCtrl.12)) %>% 
  mutate(Downregulated_NB_D6 = results.NBvCtrl.6$padj < 0.05 &
           !is.na(results.NBvCtrl.6$padj) & 
           results.NBvCtrl.6$log2FoldChange < -1.5) %>%
  mutate(Downregulated_NB_D8 = results.NBvCtrl.8$padj < 0.05 &
           !is.na(results.NBvCtrl.8$padj) & 
           results.NBvCtrl.8$log2FoldChange < -1.5) %>%
  mutate(Downregulated_NB_D10 = results.NBvCtrl.10$padj < 0.05 &
           !is.na(results.NBvCtrl.10$padj) & 
           results.NBvCtrl.10$log2FoldChange < -1.5) %>% 
  mutate(Downregulated_NB_D12 = results.NBvCtrl.12$padj < 0.05 &
           !is.na(results.NBvCtrl.12$padj) & 
           results.NBvCtrl.12$log2FoldChange < -1.5)
ggvenn(vennDat.down.SAP_D12, 
       c("Downregulated_NB_D6", "Downregulated_NB_D8", "Downregulated_NB_D10", "Downregulated_NB_D12"), 
       set_name_size = 8,
       text_size = 8)
ggsave("./5.visualization.d/Vennplot_SAP-D12_dw.png",dpi = 600, width = 30, height = 20, units = "in")

save.image("./5.visualization.Rds")




x <- list(
  SAP = results.NBvCtrl.SAP[ results.NBvCtrl.SAP$padj < 0.05 &
           !is.na(results.NBvCtrl.SAP$padj) & 
           results.NBvCtrl.SAP$log2FoldChange > 1.5,] %>% rownames(),
  D6 = results.NBvCtrl.6[ results.NBvCtrl.6$padj < 0.05 &
                               !is.na(results.NBvCtrl.6$padj) & 
                               results.NBvCtrl.6$log2FoldChange > 1.5,] %>% rownames(),
  D8 = results.NBvCtrl.8[ results.NBvCtrl.8$padj < 0.05 &
                            !is.na(results.NBvCtrl.8$padj) & 
                            results.NBvCtrl.8$log2FoldChange > 1.5,] %>% rownames(),
  D10 = results.NBvCtrl.10[ results.NBvCtrl.10$padj < 0.05 &
                            !is.na(results.NBvCtrl.10$padj) & 
                            results.NBvCtrl.10$log2FoldChange > 1.5,] %>% rownames(),
   D12 = results.NBvCtrl.12[ results.NBvCtrl.12$padj < 0.05 &
                               !is.na(results.NBvCtrl.12$padj) & 
                               results.NBvCtrl.12$log2FoldChange > 1.5,] %>% rownames()
)
ggVennDiagram(x, label = "count", set_color = "black", label_alpha = 0) + 
  scale_fill_distiller(palette = "Reds", direction = 1) + labs(title = "Upregulated in NBs")
ggsave("./5.visualization.d/Vennplot_SAP-D12_up_2024.png",dpi = 300, width = 20, height = 20, units = "cm")


# Create a named list of vectors containing genes (assuming x is your list)
gene_list <- lapply(names(x), function(name) data.frame(Gene = x[[name]], Source = name))

# Combine all the lists into one data frame
combined_genes <- do.call(rbind, gene_list)

# Find genes that appear only in one list
unique_genes_up <- combined_genes %>%
  group_by(Gene) %>%
  filter(n() == 1) %>%  # Keep only genes that appear once
  ungroup()

# View the unique genes and their sources
print(unique_genes_up)

y <- list(
  SAP = results.NBvCtrl.SAP[ results.NBvCtrl.SAP$padj < 0.05 &
                               !is.na(results.NBvCtrl.SAP$padj) & 
                               results.NBvCtrl.SAP$log2FoldChange < -1.5,] %>% rownames(),
  D6 = results.NBvCtrl.6[ results.NBvCtrl.6$padj < 0.05 &
                            !is.na(results.NBvCtrl.6$padj) & 
                            results.NBvCtrl.6$log2FoldChange < -1.5,] %>% rownames(),
  D8 = results.NBvCtrl.8[ results.NBvCtrl.8$padj < 0.05 &
                            !is.na(results.NBvCtrl.8$padj) & 
                            results.NBvCtrl.8$log2FoldChange < -1.5,] %>% rownames(),
  D10 = results.NBvCtrl.10[ results.NBvCtrl.10$padj < 0.05 &
                              !is.na(results.NBvCtrl.10$padj) & 
                              results.NBvCtrl.10$log2FoldChange < -1.5,] %>% rownames(),
  D12 = results.NBvCtrl.12[ results.NBvCtrl.12$padj < 0.05 &
                              !is.na(results.NBvCtrl.12$padj) & 
                              results.NBvCtrl.12$log2FoldChange < -1.5,] %>% rownames()
)
ggVennDiagram(x, label = "count", set_color = "black", label_alpha = 0) + 
  scale_fill_distiller(palette = "Blues", direction = 1) + labs(title = "Downregulated in NBs")
ggsave("./5.visualization.d/Vennplot_SAP-D12_dw_2024.png",dpi = 300, width = 20, height = 20, units = "cm")


# Create a named list of vectors containing genes (assuming x is your list)
gene_list <- lapply(names(y), function(name) data.frame(Gene = y[[name]], Source = name))

# Combine all the lists into one data frame
combined_genes <- do.call(rbind, gene_list)

# Find genes that appear only in one list
unique_genes_dw <- combined_genes %>%
  group_by(Gene) %>%
  filter(n() == 1) %>%  # Keep only genes that appear once
  ungroup()

# View the unique genes and their sources
print(unique_genes_dw)




unique_genes_up
unique_genes_dw

findunique <- function(unique_genes, change, stage){

change <- change
unique_genes <- unique_genes
stage <- stage
# Assuming your dataframe is named 'unique_genes', and it has columns 'Gene' and 'Source'
sap_genes <- unique_genes %>%
  filter(Source == stage) %>%
  select(Gene)  # Keep only the Gene column

# If you want to save it as a vector for further analysis:
sap_gene_vector <- sap_genes$Gene

# Assuming sap_gene_vector contains your Ensembl IDs
gene_ensembl_ids <- sap_gene_vector  # Your vector of Ensembl IDs

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_ensembl_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Symbols <- bitr(gene_ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Perform GO enrichment analysis
go_enrichment <- enrichGO(gene = entrez_ids$ENTREZID, 
                          OrgDb = org.Hs.eg.db, 
                          ont = "BP",  # Biological Process, can change to MF or CC
                          pAdjustMethod = "BH",  # Benjamini-Hochberg adjustment for p-values
                          pvalueCutoff = 0.5, 
                          qvalueCutoff = 1,
                          readable = TRUE)  # Convert IDs back to gene symbols for readable results

write_tsv(as.data.frame(go_enrichment), 
          file = paste0("./5.visualization.d/go_enrichment_results_",change,"_", stage,".tsv"))

# View the GO enrichment results
dotplot(go_enrichment, showCategory = 10,label_format = 50)  # Show top 10 categories
ggsave(paste0("./5.visualization.d/go_enrichment_results_",change, "_", stage,".png"),dpi = 300, width = 20, height = 20, units = "cm")

kegg_enrichment <- enrichKEGG(gene = entrez_ids$ENTREZID,
                              organism = "hsa",  # "hsa" for human, "mmu" for mouse
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 1)

# View the KEGG enrichment results
write_tsv(as.data.frame(kegg_enrichment), 
          file = paste0("./5.visualization.d/KEGG_enrichment_results_",change,"_", stage,".tsv"))
# Visualize KEGG pathway enrichment
dotplot(kegg_enrichment, showCategory = 10)
ggsave(paste0("./5.visualization.d/KEGG_enrichment_results_",change, "_", stage,".png"),dpi = 300, width = 20, height = 20, units = "cm")
}
findunique(unique_genes = unique_genes_dw, "dw", "D12")
findunique(unique_genes = unique_genes_dw, "dw", "SAP")



