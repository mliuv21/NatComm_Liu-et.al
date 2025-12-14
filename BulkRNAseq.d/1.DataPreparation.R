##########################################################################
#                                                                        #
#               !Data Preparation!                                       #
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
sampleinfo <- read_tsv("./SampleInfo.txt", col_types = c("cfcfc")) # Group and Diff set as factor
head(sampleinfo)

# Load the quantification data from Salmon
  # Note: The raw Salmon output files are not included in the GitHub repo due to size constraints.
  # For reproduction, please skip the tximport step and start with the provided '1.txi.count.rds'.
files <- file.path("..","5.2_salmon_TrimGalore.d","salmon_results.d", 
                   paste0("salmon_output_",sampleinfo$SampleName), "quant.sf")
files <- set_names(files, sampleinfo$SampleName)

# Load the annotation file of transcripts and genes
tx2gene <- read_tsv("./t2g_hm109_gencode.txt",col_names = c("TxID", "GeneID"))

# Read the Salmon files as count or TPM
txi.count <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T)
txi.TPM <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T,
                countsFromAbundance = "lengthScaledTPM")
# Check all the loaded data
str(txi)
head(txi.count$counts)
# Save the counts data and TPM data for later use
saveRDS(txi.count, file = "./1.txi.count.rds")
saveRDS(txi.TPM, file = "./1.txi.TPM.rds")


#### Data exploration ####
# Using counts data 
if (T){
  # Even DESeq2 could load txi object, but if we need to do the data exploration
  # we still need to do it using counts matrix
  
  # Generate raw counts data from txi
  rawCounts <- round(txi.count$counts, 0) 
  dim(rawCounts)
  # for each gene, compute total count and compare to threshold
  # keeping outcome in vector of 'logicals' (ie TRUE or FALSE, or NA)
  keep <- rowSums(rawCounts) > 5
  # summary of test outcome: number of genes in each class:
  table(keep, useNA="always") 
  
  # subset genes where test was TRUE
  filtCounts <- rawCounts[keep,]
  # check dimension of new count matrix
  dim(filtCounts)
  
  summary(filtCounts)
  
  # few outliers affect distribution visualization
  boxplot(filtCounts,
          xlab="",
          ylab="Raw counts",
          las=2,
          col=groupCols,
          main="Raw counts")
  png("./1.DataPreparation.d/RawCountsDistribution.png", bg = "transparent",
      width = 960, height = 960, units = "px")
  boxplot(filtCounts,
          xlab="",
          ylab="Raw counts",
          las=2,
          col=groupCols,
          main="Raw counts")
  dev.off()
  

  
  # Raw counts mean expression Vs standard Deviation (SD)
  plot(rowMeans(filtCounts), rowSds(filtCounts), 
       main='Raw counts: sd vs mean', 
       xlim=c(0,10000),
       ylim=c(0,5000))
  png("./1.DataPreparation.d/RawCountsSDvsMean.png", bg = "transparent",
      width = 960, height = 960, units = "px")
  plot(rowMeans(filtCounts), rowSds(filtCounts), 
       main='Raw counts: sd vs mean', 
       xlim=c(0,10000),
       ylim=c(0,5000))
  dev.off()
  
  # Get log2 counts
  logcounts <- log2(filtCounts + 1)
  # summary(logcounts[,1]) # summary for first column
  # summary(logcounts) # summary for each column
  
  # make a colour vector
  groupCols <- str_replace_all(sampleinfo$Group, c(Ctrl="#4575B4", NB="#D73027"))
  diffCols <- str_replace_all(sampleinfo$Diff,  c("NCC" = "#C7E9C0",
                                                  "SAP" = "#A1D99B",
                                                  "D6"  = "#74C476",
                                                  "D8"  = "#41AB5D",
                                                  "D10" = "#238B45",
                                                  "D12" = "#005A32"))
  
  # Check distributions of samples using boxplots
  boxplot(logcounts,
          xlab="",
          ylab="Log2(Counts)",
          las=2,
          col=groupCols,
          main="Log2(Counts)")
  png("./1.DataPreparation.d/Log2CountsDistribution.png", bg = "transparent",
      width = 960, height = 960, units = "px")
  boxplot(logcounts,
          xlab="",
          ylab="Log2(Counts)",
          las=2,
          col=groupCols,
          main="Log2(Counts)")
  # Let's add a blue horizontal line that corresponds to the median
  abline(h=median(logcounts), col="blue")
  dev.off()
  # Log2 counts standard deviation (sd) vs mean expression
  plot(rowMeans(logcounts), rowSds(logcounts), 
       main='Log2 Counts: sd vs mean')
  png("./1.DataPreparation.d/Log2CountsSDvsMean.png", bg = "transparent",
      width = 960, height = 960, units = "px")
  plot(rowMeans(logcounts), rowSds(logcounts), 
       main='Log2 Counts: sd vs mean')
  dev.off()
  
  
  ## VST ##
  vst_counts <- vst(filtCounts, blind = TRUE)
  
  
  # Check distributions of samples using boxplots
  png("./1.DataPreparation.d/VSTCountsDistribution.png", bg = "transparent",
      width = 960, height = 960, units = "px")
  boxplot(vst_counts, 
          xlab="", 
          ylab="VST counts",
          las=2,
          col=groupCols)
  
  # Let's add a blue horizontal line that corresponds to the median
  abline(h=median(vst_counts), col="blue")
  dev.off()
  
  # VST counts standard deviation (sd) vs mean expression
  png("./1.DataPreparation.d/VSTCountsSDvsMean.png", bg = "transparent",
      width = 960, height = 960, units = "px")
  plot(rowMeans(vst_counts), rowSds(vst_counts), 
       main='VST counts: sd vs mean')
  dev.off()
  
  
  ## PCA using vst
  # run PCA
  pcDat.vst <- prcomp(t(vst_counts))
  # plot PCA
  autoplot(pcDat.vst,
           data = sampleinfo, 
           colour=groupCols, 
           size=5,
           main = "PCA plot (VST)") +
    geom_label_repel(aes(label = SampleName))
  
  ggsave("./1.DataPreparation.d/PCA_VST.png",
         width = 2800, height = 2800, units = "px")
  
  ## PCA using rlog
  rlogcounts <- rlog(filtCounts)
  # run PCA
  pcDat.rlog <- prcomp(t(rlogcounts))
  # plot PCA
  autoplot(pcDat.rlog,
           data = sampleinfo, 
           colour=groupCols, 
           size=5,
           main = "PCA plot (rlog)") +
    geom_label_repel(aes(label = SampleName))
  ggsave("./1.DataPreparation.d/PCA_rlog.png",
         width = 2800, height = 2800, units = "px")
}


