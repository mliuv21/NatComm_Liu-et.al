##########################################################################
#                                                                        #
#                  !Pathway enrichment analysis!                         #
#                    Mingzhi Liu                                         #
#                   2022 - 03 - 02                                       #
#                                                                        #
#                                                                        #
##########################################################################

if (TRUE) {
  # Remove everything in the environment
  rm(list = ls())
  
  # Load essential packages
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(gridExtra)
}



# Load the data 
shrink.NBvsCtrl.d12 <- readRDS("./4.batchanalysis.d/ddshrinkResult_NBvsCtrl_12.Rds")
result.NBvsCtrl.d12 <- readRDS("./4.batchanalysis.d/Results_NBvsCtrl_12_annotation.Rds")
result.NBvsCtrl.d10 <- readRDS("./4.batchanalysis.d/Results_NBvsCtrl_10_annotation.Rds")
result.NBvsCtrl.d8 <- readRDS("./4.batchanalysis.d/Results_NBvsCtrl_8_annotation.Rds")
result.NBvsCtrl.d6 <- readRDS("./4.batchanalysis.d/Results_NBvsCtrl_6_annotation.Rds")
result.NBvsCtrl.SAP <- readRDS("./4.batchanalysis.d/Results_NBvsCtrl_SAP_annotation.Rds")
result.NBvsCtrl.NCC <- readRDS("./4.batchanalysis.d/Results_NBvsCtrl_NCC_annotation.Rds")
resultsList <- list(NCC = result.NBvsCtrl.NCC, SAP = result.NBvsCtrl.SAP,
                    D6 = result.NBvsCtrl.d6, D8 = result.NBvsCtrl.d8,
                    D10 = result.NBvsCtrl.d10, D12 = result.NBvsCtrl.d12)

# Find sig genes
sigGenes <- shrink.NBvsCtrl.d12 %>% 
  drop_na(Entrez, FDR) %>%
  filter(FDR < 0.05 & abs(logFC) > 1) %>%
  pull(Entrez)

for (i in 1:length(resultsList)) {
  list.name <- names(resultsList[i])
  
  print(list.name)
  
  sigGenes.df <- resultsList[[i]] %>% 
    dplyr::filter(!is.na(Entrez)) %>%
    filter(FDR < 0.05 & abs(logFC) > 1)
  sigGenes.df2 <- data.frame(lapply(sigGenes.df, as.character), stringsAsFactors=FALSE)
  
  write.table(sigGenes.df2, file = paste0("./6.enrichment.d/sigGenes_NBvsCtrl", list.name, "_ALL.txt"), 
              sep = "\t", quote = FALSE, row.names = F)
  
  sigGenes <- sigGenes.df %>%
    pull(Entrez)
  
  print("sigGenes extract successfully!")
  
  goRes <- enrichGO(gene = sigGenes,
                   ont = "BP",
                   OrgDb = "org.Hs.eg.db",
                   pvalueCutoff = 1)
  dotplot(goRes) + ggtitle(paste0("GO BP NB vs Ctrl ",list.name, " ALL"))
  ggsave(paste0("./6.enrichment.d/GOBP", list.name, "_ALL.pdf"), 
         height = 9, width = 9, units = "in",dpi = 1200 )
  
  print(paste0(list.name, " has been finished!"))
}
  

sigGenes <- result.NBvsCtrl.d12 %>% 
  dplyr::filter(!is.na(Entrez)) %>%
  filter(FDR < 0.05 & abs(logFC) > 1) %>%
  pull(Entrez)

sigGenes <- result.NBvsCtrl.d12 %>% 
  dplyr::filter(!is.na(Entrez)) %>%
  filter(FDR < 0.05 & logFC > 1) %>%
  pull(Entrez)


enrichment_analysis <- function(day = "12", change = "up") {
  print(paste0("Now we are doing the enrichment analysis on day ", day, " with ", change, " condition!"))
  
  result <- readRDS(paste0("./4.batchanalysis.d/Results_NBvsCtrl_",day,"_annotation.Rds"))
  if (change == "up") {
    sigGenes <- result %>% 
      dplyr::filter(!is.na(Entrez)) %>%
      filter(FDR < 0.05 & logFC> 1) %>%
      pull(Entrez)
  }else{
    if (change == "down") {
      sigGenes <- result %>% 
        dplyr::filter(!is.na(Entrez)) %>%
        filter(FDR < 0.05 & logFC< -1) %>%
        pull(Entrez)
    }else{
      sigGenes <- result %>% 
        dplyr::filter(!is.na(Entrez)) %>%
        filter(FDR < 0.05 & abs(logFC) > 1) %>%
        pull(Entrez)
    }
  }
  print("Now, GO analysis is doing!")
  goRes <- enrichGO(gene = sigGenes,
                    ont = "BP", pAdjustMethod = "fdr",
                    OrgDb = "org.Hs.eg.db",
                    pvalueCutoff = 0.2,
                    readable = TRUE)
  write_tsv(as.data.frame(goRes), paste0("./6.enrichment.d/GOBPNBvsCtrl_",day,"_", change, ".tsv"))
  
  dotplot(goRes, showCategory = 15, label_format = 50) + 
    ggtitle(paste0("GO BP NB vs Ctrl Day ", day," ", change))
  ggsave(paste0("./6.enrichment.d/GOBPNBvsCtrl",day,change,"_2024.pdf"), 
         height = 12, width = 10, units = "in",dpi = 600 )
  
  goRes <- enrichGO(gene = sigGenes,
                    ont = "CC", pAdjustMethod = "BH",
                    OrgDb = "org.Hs.eg.db",
                    pvalueCutoff = 0.05,
                    readable = TRUE)
  write_tsv(as.data.frame(goRes), paste0("./6.enrichment.d/GOCCNBvsCtrl_",day,"_", change, "_CC.tsv"))
  
  dotplot(goRes, showCategory = 15, label_format = 50) + 
    ggtitle(paste0("GO CC NB vs Ctrl Day ", day," ", change))
  ggsave(paste0("./6.enrichment.d/GOCCNBvsCtrl",day,change,"_2024.pdf"), 
         height = 12, width = 10, units = "in",dpi = 600 )
  goplot(goRes, showCategory = 15) + 
    ggtitle(paste0("GO CC NB vs Ctrl Day ", day," ", change))
  ggsave(paste0("./6.enrichment.d/GOPLOTCCNBvsCtrl",day,change,"_2024.pdf"), 
         height = 12, width = 10, units = "in",dpi = 600 )
  print("GO analysis has been done!")
  print("Now KEGG enrichment is ongoing!")
  KEGGRes <- enrichKEGG(gene =sigGenes,
                        organism = "hsa",
                        qvalueCutoff = 1, 
                        pvalueCutoff = 0.5)
  write_tsv(as.data.frame(KEGGRes), paste0("./6.enrichment.d/KEGGNBvsCtrl_",day,"_", change, ".tsv"))
  
  dotplot(KEGGRes, showCategory = 20, label_format = 50) + 
    ggtitle(paste0("KEGG NB vs Ctrl Day ", day," ", change))
  ggsave(paste0("./6.enrichment.d/KEGGNBvsCtrl_",day,change,"_2024.pdf"), 
         height = 12, width = 10, units = "in",dpi = 600 )
  
}

days <- list("SAP", "6", "8", "10", "12")
changes <- list("up", "down", "all")

for (d in days) {
  tryCatch({
    print(d)
    day <- d
    for (c in changes) {
      tryCatch({
        print(c)
        change <- c
        enrichment_analysis(day, change)
      }, error=function(e){})
    }
  }, error=function(e){})
}




enrichment_analysis( change = "down")

# GO enrichment

goRes <- enrichGO(gene = sigGenes,
                  ont = "BP",
                  OrgDb = "org.Hs.eg.db",
                  pvalueCutoff = 0.5)
summary(goRe)

dotplot(goRe, showCategory = 15) + ggtitle(paste0("GO BP NB vs Ctrl D12 UP"))
ggsave("./6.enrichment.d/GOBPNBvsCtrl12UP.pdf", 
       height = 9, width = 9, units = "in",dpi = 1200 )


KEGGRes <- enrichKEGG(gene =sigGenes,
                      organism = "hsa",qvalueCutoff = 1, pvalueCutoff = 0.5)
dotplot(KEGGRes) 






###### Check in the pathway what kinds of genes are contribute to ######
entrezidList <- goRes %>% as.data.frame() %>% 
  dplyr::filter(grepl("axon", Description)) %>% 
  dplyr::pull(geneID) %>% 
  str_split( "/") %>% 
  unlist()
selectedGenes <- result.NBvsCtrl.d12 %>% 
  dplyr::filter(Entrez %in% entrezidList) %>% 
  dplyr::select(1:8)

write.table(selectedGenes, "./6.enrichment.d/gsea.d/P53_related_genes_D12_DW.txt", 
            sep = "\t", quote = F, row.names = F)


# GSEA

gsea.analysis <- function(shrink.result, gsea.category, gsea.subcategory=NA, regulation ="up" ){
  gseaCat = gsea.category
  sub = gsea.subcategory
  comparename <- deparse(match.call()$shrink.result)
  
  print(!is.na(sub))
  if (!is.na(sub)){
    term2gene <- msigdbr(species = "Homo sapiens",
                         category = gseaCat,
                         subcategory = sub) %>%
      dplyr::select(gs_name, entrez_gene) %>%
      distinct()
    term2name <- msigdbr(species = "Homo sapiens",
                         category = gseaCat,
                         subcategory = sub) %>%
      dplyr::select(gs_name, gs_description) %>%
      distinct()
    outdir <- paste0("./6.enrichment.d/gsea.d/", comparename, "/", gseaCat, "/", sub, "/")
    print(outdir)  
  } else {
    term2gene <- msigdbr(species = "Homo sapiens",
                         category = gseaCat) %>%
      dplyr::select(gs_name, entrez_gene) %>%
      distinct()
    term2name <- msigdbr(species = "Homo sapiens",
                         category = gseaCat) %>%
      dplyr::select(gs_name, gs_description) %>%
      distinct()
    outdir <- paste0("./6.enrichment.d/gsea.d/", comparename, "/", gseaCat, "/")
    print(outdir)
  }
  print("TERMs are ok")
  if (regulation =="up") {
    rankedGenes <- shrink.result %>% 
      drop_na(Entrez) %>%
      mutate(rank = logFC) %>%
      arrange(desc(rank)) %>%
      pull(rank, Entrez)
  } else {
    rankedGenes <- shrink.result %>% 
      drop_na(Entrez) %>%
      mutate(rank = logFC) %>%
      arrange(desc(-rank)) %>%
      pull(rank, Entrez)
    outdir <- paste0(outdir, "/", "down", "/")
    print(outdir)
  }

  gseaRes <- GSEA(rankedGenes,
                  TERM2GENE = term2gene,
                  TERM2NAME = term2name)
  print(outdir)
  outdir <- str_replace_all(outdir, ":", "_")
  
  
  if (!file.exists(outdir)){
    dir.create(file.path(outdir), recursive = T)
  } 
  gsea_table <- tibble(as.data.frame(gseaRes)) %>%
    dplyr::arrange(desc(abs(NES))) %>%
    dplyr::top_n(50, wt = -p.adjust) %>%
    dplyr::select(-core_enrichment) %>%
    dplyr::mutate(across(c("enrichmentScore", "NES"), round, digits = 3)) %>%
    dplyr::mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))
  write_tsv(gsea_table, file = paste0(outdir, "GSEA_list_top10.txt"))
  
  print("oho") 
  for (i in 1:10) {
    gs <- as.data.frame(gseaRes)$ID[i]
    print(paste0("oho",i)) 
    #pdf(filename = paste0(outdir,"/", gs, "_", i, ".pdf"))
    gsp <- gseaplot2(gseaRes,
                     geneSetID = gs,
                     title = gs,
                     base_size = 20)
    ggsave(paste0(outdir,"/", gs, "_", i, ".pdf"),
           plot = marrangeGrob(gsp, nrow=3, ncol=1),
           width = 18, height = 20, units = "in")
    #dev.off()
    Sys.sleep(3)
  }
}

gsea.analysis(shrink.result = result.NBvsCtrl.d12, "H")
gsea.analysis(shrink.result = result.NBvsCtrl.d10, "H")
gsea.analysis(shrink.result = result.NBvsCtrl.d8, "H")
gsea.analysis(shrink.result = result.NBvsCtrl.d6, "H")
gsea.analysis(shrink.result = result.NBvsCtrl.SAP, "H")


gsea.analysis(shrink.result = result.NBvsCtrl.d10, "C2", gsea.subcategory = "CP:KEGG", regulation = "up")
gsea.analysis(shrink.result = result.NBvsCtrl.NCC, "H")
gsea.analysis(shrink.result = result.NBvsCtrl.NCC, "C2", gsea.subcategory = "CP:REACTOME", regulation = "up")

REACTOME
KEGG


##### goCombine #### 20250616

# Have been save the UP regulated goRes from day 8, 10 and 12 as goRes_8, goRes_10 and goRes_12
goRes_df_8 <- goRes_8@result
goRes_df_10 <- goRes_10@result
goRes_df_12 <- goRes_12@result

top_num <- 8

top8_sel <- goRes_df_8 %>% arrange(p.adjust) %>% slice_head(n = top_num)%>%
  mutate(Group = "D8")
top10_sel <- goRes_df_10 %>% arrange(p.adjust) %>% slice_head(n = top_num) %>%
  mutate(Group = "D10")
top12_sel <- goRes_df_12 %>% arrange(p.adjust) %>% slice_head(n = top_num)%>%
  mutate(Group = "D12")
top_pathways <- unique(c(top8_sel$Description, top10_sel$Description, top12_sel$Description))

top8_filtered <- goRes_df_8 %>%
  filter(Description %in% top_pathways) %>%
  mutate(Group = "D8")

top10_filtered <- goRes_df_10 %>%
  filter(Description %in% top_pathways) %>%
  mutate(Group = "D10")

top12_filtered <- goRes_df_12 %>%
  filter(Description %in% top_pathways) %>%
  mutate(Group = "D12")

combined_df <- bind_rows(top8_filtered, top10_filtered, top12_filtered)


ggplot(combined_df, aes(y = Group, x = Description, size = Count, color = qvalue)) +
  geom_point() +
  scale_color_distiller("GnBu",direction = -1, name = "-log10(FDR)") +
  scale_size_continuous(name = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle =45, hjust = 1)) +
  labs(x = "Group", y = "Pathway")

combined_df <- bind_rows(top8_sel, top10_sel, top12_sel)
combined_df$geneID <- gsub("\\/",",", combined_df$geneID)
combined_df$CollectgeneID <- sapply(strsplit(combined_df$geneID, ","), function(x) paste(head(x, 10), collapse = ","))
combined_df$Size <- combined_df$Count
combined_df <- combined_df %>% dplyr::arrange(desc(Count))
combined_df$Group <- factor(combined_df$Group, levels = c("D8", "D10", "D12"))
ggplot (combined_df)+
  geom_bar( aes(x = -log10(qvalue), y= interaction(Description, Group),
                            fill = Group),stat="identity") +
  scale_fill_manual(values = c("#4CAF50FF", "#F49600FF", "#C62828FF"), name="Days") + 
  geom_text(aes(x = 0.1, y = interaction(Description, Group), label=Description),
            size =3, hjust = 0, color= "black") +
  #geom_text(aes(x = 0.1, y = interaction(Description, Group), label=CollectgeneID),
  #          size =2, hjust = 0, vjust=2.5, color= "black") +
  geom_point(aes(x = -3, y = interaction(Description, Group), size = Count, fill = Group), shape = 21) +
  geom_text(aes(x=-1,y=interaction(Description, Group), label=Count), size =3) +
  scale_size(range=c(4,8), guide=guide_legend(override.aes=list(fill="black")))+
  #scale_x_continuous(expand=expansion(mult=c(0,0.2)),limits = c(-2,48))+
  guides(fill=guide_legend(reverse=T))+
  labs(x="-log10(FDR)",y="Description") + theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        legend.frame  = element_rect(colour = "black"))
  
ggsave(paste0("./6.enrichment.d/KEGGNBvsCtrl_CombineD8-12_5","_2025.pdf"), 
       height = 10, width = 8, units = "in")
ggsave(paste0("./6.enrichment.d/KEGGNBvsCtrl_CombineD8-12_5","_2025_Width.pdf"), 
       height = 8, width = 12, units = "in")


# 1.
go_simplified <- simplify(goRes_8, cutoff = 0.7, by = "qvalue", select_fun = min)
dotplot(goRes_8) + dotplot(go_simplified)
  

#2.
library(GOSemSim)

go_terms <- goRes_8@result$ID

goSimMatrix <- termSim(go_terms, go_terms, semData = godata('org.Hs.eg.db', ont = "BP"))

hc <- hclust(as.dist(1 - goSimMatrix), method = "ward.D2")

clusters <- cutree(hc, k = 5) 

goRes_df_8$Cluster <- clusters

goRes_df_8$Cluster %>% table

goRes_df_8 %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  select(Cluster, Description, qvalue) %>%
  filter(Cluster==5)
