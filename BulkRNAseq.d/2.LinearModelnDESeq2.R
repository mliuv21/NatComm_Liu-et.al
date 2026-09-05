##########################################################################
#                                                                        #
#               Linear Model & DESeq2                                    #
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
  library(tximport)
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(genefilter)
}
#### Data loading ####
# Load the sample information 
sampleinfo <- read_tsv("./SampleInfo.txt", col_types = c("cfcfc"))
txi <- readRDS("./1.txi.count.rds")
groupCols <- str_replace_all(sampleinfo$Group, c("Ctrl"="#4575B4", "NB"="#D73027"))
diffCols <- str_replace_all(sampleinfo$Diff,  c("NCC" = "#C7E9C0",
                                                "SAP" = "#A1D99B",
                                                "D6"  = "#74C476",
                                                "D8"  = "#41AB5D",
                                                "D10" = "#238B45",
                                                "D12" = "#005A32"))
# Check all counts data are all in and in correct order in sampleinfo
all(colnames(txi$counts)==sampleinfo$SampleName)


#### Create the linear model for the analysis ####
# As same as as.formula( ~ Days + Group + Days:Group)
diff_group.model <- as.formula( ~ Diff * Group)
model.matrix(diff_group.model, data = sampleinfo)


coldata <- mutate(sampleinfo, Diff = fct_relevel(Diff, "NCC"), 
                  Group = fct_relevel(Group, "Ctrl"))


#### Creat DESeq2 object ####
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = coldata,
                                       design = diff_group.model)
# Filter out the unexpressed genes
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

# ## Three main step of DESeq2 to analysis
# ddsObj <- estimateSizeFactors(ddsObj.filt)
# 
# ddsObj <- estimateDispersions(ddsObj)
# plotDispEsts(ddsObj)
# 
# ddsObj <- nbinomWaldTest(ddsObj)
ddsObj <- DESeq(ddsObj.filt)

# Generate a result table
results.simple <- results(ddsObj, alpha=0.05)
results.simple

# Check the contrast of the matirx
ddsObj %>% resultsNames()

### PCA Function##
plotPCAml = function(object, intgroup="condition",
                                  ntop=500, returnData=FALSE, pcsToUse=1:2)
{
  message(paste0("using ntop=",ntop," top features by variance"))
  
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]), scale.=T)
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  pcs <- paste0("PC", pcsToUse)
  d <- data.frame(V1=pca$x[,pcsToUse[1]],
                  V2=pca$x[,pcsToUse[2]],
                  group=group, name=colnames(object), colData(object))
  colnames(d)[1:2] <- pcs
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pcsToUse]
    return(d)
  }
  
  ggplot(data=d, aes_string(x=pcs[1], y=pcs[2], color="group")) +
    geom_point(size=3) + 
    xlab(paste0(pcs[1],": ",round(percentVar[pcsToUse[1]] * 100),"% variance")) +
    ylab(paste0(pcs[2],": ",round(percentVar[pcsToUse[2]] * 100),"% variance")) +
    coord_fixed()
}


### PCA analysis 
if (T) {
  vstcounts <- vst(ddsObj, blind = TRUE)
  pca_dt <- plotPCAml(vstcounts, ntop=500,intgroup = c("Group", "Diff"), returnData = T,pcsToUse = c(1,2))
  percentVar <- round(100 * attr(pca_dt, "percentVar"))
  theme_set(theme_minimal(base_size = 8))
  pca <- ggplot(pca_dt, aes(PC1, PC2)) +
    #geom_label_repel(aes(color=Group, label = name), size = 4,   show.legend = F,
    #                nudge_x = 0, nudge_y = 2) +
    geom_point(aes(color=Group),size = 6) +
    scale_color_manual(values=c("#4575B4","#D73027")) + 
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_minimal()+
    theme(#panel.background = element_rect(fill='transparent', color="grey"),
      plot.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill='transparent'),
      legend.box.background = element_rect(fill='transparent')) ;pca
  
  ###
# Alternative PCA visualization retained for reference.
# pca <- ggplot(pca_dt) +
#   geom_arc_bar(
#     aes(
#       x0 = PC1, y0 = PC2,
#       r0 = 0, r = 1,
#       start = pi/2, end = 3*pi/2,
#       fill = Group), color = NA,alpha = 1) +
#   geom_arc_bar(
#     aes(
#       x0 = PC1, y0 = PC2,
#       r0 = 0, r = 1,
#       start = -pi/2, end = pi/2,
#       fill = Diff
#     ), color = NA,alpha = 0.6) +coord_fixed()+
#   scale_fill_manual(values=c("Ctrl"= "#4575B4","NB" = "#D73027", "NCC" = "#3D98D3FF",
#                              "SAP" = "#9C27B0FF",
#                              "D6"  = "#8F7289FF",
#                              "D8"  = "#4CAF50FF",
#                              "D10" = "#F49600FF",
#                              "D12" = "#C62828FF")) + 
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_minimal()+
#   theme(#panel.background = element_rect(fill='transparent', color="grey"),
#     plot.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_rect(fill='transparent'),
#     legend.box.background = element_rect(fill='transparent')) ;pca
  ###
  
  pc1.density <- ggplot(pca_dt) +
    geom_density(aes(PC1, group = Diff, fill = Diff),
                 color= "black", alpha = 0.6, position = "identity", show.legend = T) + 
    scale_fill_manual(values=c("#3D98D3FF", "#9C27B0FF", "#8F7289FF", "#4CAF50FF", "#F49600FF", "#C62828FF"))+ 
    scale_linetype_discrete() + coord_cartesian(xlim = c(-35,20))+
    scale_y_discrete(expand = c(0,0.001)) + 
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),    
          plot.background = element_blank(),
          panel.grid.minor = element_blank());pc1.density
  pc2.density <- ggplot(pca_dt) +
    geom_density(aes(PC2, group = Group, fill = Group),
                 color= "black", alpha = 0.6, position = "identity", show.legend = T) + 
    scale_fill_manual(values=c("#4575B4","#D73027")) + 
    scale_linetype_discrete() + coord_cartesian(xlim = c(-25,15))+
    scale_y_discrete(expand = c(0,0.001)) + 
    theme(axis.title.y = element_blank(), 
          axis.ticks.y = element_blank(),    
          plot.background = element_blank(),
          panel.grid.minor = element_blank()) + coord_flip();pc2.density

  
  p1 <- pca %>% insert_top(pc1.density, height = 0.3) %>% insert_right(pc2.density, width = 0.3) %>% as.ggplot();p1
  ggsave("./2.DESeq2.d/PCA_plot_Density_NoLabel_half.pdf", height = 10, width = 10, units = "in", bg='transparent')
 }

pca_3d <- pca_dt[,1:2]
pca_dt <- plotPCA(vstcounts, ntop=1000,intgroup = c("Group", "Diff"), returnData = T,pcsToUse = c(2,3))
pca_dt -> a
a$PC1 <- pca_3d$PC1
fig <-plotly::plot_ly(a, text = a$name, x = a$PC1, y = a$PC2, z=a$PC3, 
                color = a$Group, colors = c("#4575B4","#D73027"))




saveRDS(ddsObj,"2.ddsObj_diff-group.Rds")
