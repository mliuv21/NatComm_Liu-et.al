##########################################################################
#                                                                        #
#               !Data Preparation! --Additation                          #
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
  library(ggfortify)
  library(ggrepel)
}
#### Data loading ####
# Load the sample information 
sampleinfo <- read_tsv("./SampleInfo.txt", col_types = c("cfcfc"))
txi <- readRDS("./1.txi.count.rds")
rawCounts <- round(txi$counts, 0)
groupCols <- str_replace_all(sampleinfo$Group, c("Ctrl"="#4575B4", "NB"="#D73027"))
diffCols <- str_replace_all(sampleinfo$Diff,  c("NCC" = "#C7E9C0",
                                                "SAP" = "#A1D99B",
                                                "D6"  = "#74C476",
                                                "D8"  = "#41AB5D",
                                                "D10" = "#238B45",
                                                "D12" = "#005A32"))

#### Additional Data Expolre ####
### Count density
keep <- rowSums(rawCounts) > 5
logCounts <- log2(rawCounts[keep, ] + 1)
logCounts %>% 
  as.data.frame() %>% 
  pivot_longer(names_to = "SampleName", values_to = "logCounts", everything()) %>% 
  ggplot(aes(x=logCounts, group = SampleName)) +
  geom_density(aes(colour = SampleName)) +
  labs(x = "log2(Counts)", title = "Raw count density")
ggsave("./1.DataPreparation.d/Rawcount_density_samples.png",
       width = 2800, height = 2800, units = "px")

logCounts %>% 
  as.data.frame() %>% 
  pivot_longer(names_to = "SampleName", values_to = "logCounts", everything()) %>% 
  left_join(sampleinfo) %>% 
  ggplot(aes(x=logCounts, group = SampleName )) +
  geom_density(aes(colour = Group)) +
  scale_colour_manual(values = c("#4575B4","#D73027")) +
  labs(x = "log2(Counts)", title = "Raw count density (log2)")
ggsave("./1.DataPreparation.d/log2count_density_samples.png",
       width = 2800, height = 2800, units = "px")

rlogCnts <- rlog(rawCounts[keep, ])
rlogCnts %>% 
  as.data.frame() %>% 
  pivot_longer(names_to = "SampleName", values_to = "logCounts", everything()) %>% 
  left_join(sampleinfo) %>% 
  ggplot(aes(x=logCounts, group = SampleName)) +
  geom_density(aes(colour = Group)) +
  scale_colour_manual(values = c("#4575B4","#D73027")) +
  labs(x = "log2(Counts)", title = "Raw count density_ rlog Normalized")
ggsave("./1.DataPreparation.d/rlogCnt_density.png",
       width = 2800, height = 2800, units = "px")

vstCnts <- vst(rawCounts[keep, ])
vstCnts%>% 
  as.data.frame() %>% 
  pivot_longer(names_to = "SampleName", values_to = "logCounts", everything()) %>% 
  left_join(sampleinfo) %>% 
  ggplot(aes(x=logCounts, group = SampleName)) +
  geom_density(aes(colour = Group)) +
  scale_colour_manual(values = c("#4575B4","#D73027")) +
  labs(x = "log2(Counts)", title = "Raw count density_VST Normalized")
ggsave("./1.DataPreparation.d/Vstcount_density.png",
       width = 2800, height = 2800, units = "px")

### Hierachical clustering
library(ggdendro)
hclDat <-  t(rlogCnts) %>%
  dist(method = "euclidean") %>%
  hclust()
ggdendrogram(hclDat, rotate=TRUE)
ggsave("./1.DataPreparation.d/Dendrogram_samples.png",
       width = 2800, height = 2800, units = "px")

dendro.dat <-as.dendrogram(hclDat) %>% dendro_data()
dendro.dat$labels <- dendro.dat$labels %>%
  left_join(sampleinfo, by = c(label = "SampleName"))

ggplot(dendro.dat$segment) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_label(data = dendro.dat$labels,
             aes(x = x,
                 y = y,
                 label = label,
                 fill = Group),
             hjust = 0,
             nudge_y = 1) +
  scale_fill_manual(values = c("#4575B4","#D73027")) +
  coord_flip() +
  labs(x = NULL, y = "Distance", title = NULL, fill = "Group") +
  scale_y_reverse(expand = c(0.3, 0)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank())
ggsave("./1.DataPreparation.d/Dendrogram_samples_groups.png",
       width = 2800, height = 2800, units = "px")

### Correlation plot
library(corrplot)
corDat <- cor(rlogCnts)
rownames(corDat) <- sampleinfo$SampleName
colnames(corDat) <- sampleinfo$SampleName
corrplot(corDat, 
         col = rev(COL2('RdYlBu', 100)),
         method = "color", 
         addCoef.col = "black", 
         number.digits= 3,
         order = "hclust",
         is.corr = FALSE,
         tl.cex = 5,
         number.cex = 1,
         cl.cex	= 3,
         cl.length = 7,
         col.lim = c(0.94, 1),
         sig.level = 0.05,
         tl.col = "black")
png("./1.DataPreparation.d/CorrelationPlot_.png",
       width = 2800, height = 2800, units = "px")
dev.off()



#################################################################
sessionInfo()
