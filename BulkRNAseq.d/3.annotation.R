##########################################################################
#                                                                        #
#               !!Annotation!!                                           #
#                Mingzhi Liu                                             #
#               2022 - 03 - 02                                           #
#                                                                        #
#                                                                        #
##########################################################################

#### Analysis environment setup #### 
if (TRUE) {
  # Remove everything in the environment
  rm(list = ls())
  
  
  # Load essential packages
  library(DESeq2)
  library(dplyr)
  library(tidyverse)
  library(AnnotationHub)
  library(AnnotationDbi)
  library(ensembldb)
}

ddsObj <- readRDS("./2.ddsObj_diff-group.Rds")

### Using ddsObj comparision set as the reference
resultsNames(ddsObj)
### Extract the genes are differnt betweeen NB vs Ctrl in D12
results.NBvCtrl.12 <- results(ddsObj, 
                              contrast = list(c("Group_NB_vs_Ctrl", "DiffD12.GroupNB")),
                              alpha=0.05)

sum(results.NBvCtrl.12$padj < 0.05, na.rm = TRUE)
table(results.NBvCtrl.12$padj < 0.05 & results.NBvCtrl.12$log2FoldChange < -2)

head(results.NBvCtrl.12)

## Annotation!
# Query the database
ah <- AnnotationHub()
ah
saveRDS(ah,"./3.annotationhub_dataset.Rds")
ah[1]
# Since our sequencing is based on human, filter the genome database of human 
# Meantime, since later we need to analysis with clusterProfiler which need OrgDb type gene info
# So together filter out with OrgDb like dataset
hm <- query(ah, c("OrgDb", "Homo Sapiens"))[[1]] # There is only 1 record here which is what we need.
hmEnsem <- query(ah, c("EnsDb", "Homo Sapiens", "108", "GTF")) # AH109336 is the name of the data in annotationHub, download it
hmEnsem <- ah[["AH109336"]]

saveRDS(hm,"OrgDb_HomoSapiens.Rds")
saveRDS(hmEnsem, "Ensembel_Homosapiens.Rds")
# Extract the gene out
annotations <- genes(hmEnsem, return.type = "data.frame")
colnames(annotations)

annot <- annotations %>% 
  dplyr::select(gene_id, gene_name, entrezid) %>%
  dplyr::filter(gene_id %in% rownames(results.NBvCtrl.12))

head(annot)
length(annot$entrezid)
length(unique(annot$entrezid))
table(is.na(annot$entrezid))

annot.NBvCtrl.12 <- as.data.frame(results.NBvCtrl.12) %>%
  rownames_to_column("gene_id") %>% 
  left_join(annotations, "gene_id") %>%
  rename(logFC = log2FoldChange, FDR = padj,GeneID = gene_id,
         Entrez = entrezid, Symbol = gene_name)

head(annot.NBvCtrl.12)
write_tsv(annot.NBvCtrl.12, "./3.annotation.d/NBvsCtrl.D12_Results_Annotated.tsv")

saveRDS(annot.NBvCtrl.12, "annot.NBvCtrl.12.Rds")

