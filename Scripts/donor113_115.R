library(DESeq2)
library(apeglm)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(ReactomePA)
library(pheatmap)
library(AnnotationDbi)

countData <- read.csv('/projects/b1042/HultquistLab/Daphne/Hultquist02_9.1.2022/gene_count_matrix.csv',row.names = 1, check.names = FALSE)

colData <- read.csv('/projects/b1042/HultquistLab/Daphne/Hultquist02_9.1.2022/coldata.csv',row.names = 1)

#extracting gene symbol
new_rn <- c()
for (rn in row.names(countData)){
  tmp <- strsplit(rn, "|", fixed = TRUE)[[1]][2]
  new_rn <- c(new_rn, tmp)
}



#removing duplicated genes and also a missing name
countData<-countData[!duplicated(new_rn),]

new_rn <- new_rn[!duplicated(new_rn)]

todrop <- which(is.na(new_rn))
print(todrop)
new_rn <- new_rn[!is.na(new_rn)]

countData <- countData[-todrop,] #missing gene name

row.names(countData) <- new_rn

#CPSF6 vs NT for all 3 donors individually
# for (i in c("113","114","115")){

  #data filtering
#   tmpcount <- dplyr::select(countData, matches(i))
#   tmpcount <- dplyr::select(tmpcount, matches("CPSF6|NT"))
#   tmpcoldata <- colData[colData$design!="TRIM5a",]
#   tmpcoldata = tmpcoldata[(tmpcoldata$donor==as.numeric(i)),] #|(tmpcoldata$donor==114)
#   tmpcoldata$donor <- as.factor(tmpcoldata$donor)
  
  #creating a DESeqDataSet subclass from a given matrix

#   dds <- DESeqDataSetFromMatrix(countData = tmpcount,
#                                 colData = tmpcoldata,
#                                 design= ~ design)
# 
    #keeping proteins that have at least a count of 10 for at least 3 treated samples
#   smallestGroupSize <- 3
#   keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
#   dds <- dds[keep,]
# 
# 
#   #performing variance stabilizing transformation and plotting PCA
#   vsd <- vst(dds, blind=FALSE)
#   plotPCA(vsd, intgroup=c( "design"))
# 
#   #running DESeq
#   dds <- DESeq(dds)
#   res <- results(dds, contrast = c("design","CPSF6","NT") )
# 
#   #getting the top X genes by adjusted p-value
#   # topgenes<- data.matrix(tmpcount[rownames(res[order(res$padj, decreasing = FALSE),])[0:25],])
#   # logtopgenes <- log(topgenes)
#   # pheatmap(logtopgenes)
# 
  #plotting the volcano plot
#   EnhancedVolcano(res,
#                   lab = rownames(res),
#                   x = 'log2FoldChange',
#                   y = 'padj')
#   ggsave(paste(i, "volcanoplot.pdf", sep = '_'))
# 
# 
#   #filtering for the up and downregulated genes
#   dds_results_filtered <-res[complete.cases(res),]
#   upreg <- dds_results_filtered[ dds_results_filtered$log2FoldChange > 0,]
#   downreg <- dds_results_filtered[ dds_results_filtered$log2FoldChange < 0,]
#   diffreg <- dds_results_filtered[ abs(dds_results_filtered$log2FoldChange) > 0,]
#   datalist <- list(upreg, downreg, diffreg)
#   datalist_names <- c("upreg", "downreg","diffreg")
#   for (d in 1:length(datalist)){
#     data <- datalist[d]
#     name <- datalist_names[d]
#     write.csv(data, paste(paste(i, name, sep = "notsig_"), ".csv", sep = ""))
#   }
# }

## Comparing CPSF6 or TRIM5a to NT for all 3 patients combined
# tmpcount <- dplyr::select(countData, matches("CPSF6|NT"))
# tmpcount <- tmpcount[,1:12] #only donors 113,114,115
# tmpcoldata <- colData[colData$design!="TRIM5a",]
# tmpcoldata <- tmpcoldata[1:12,]
# tmpcoldata$donor <- as.factor(tmpcoldata$donor)

#including trim5
# tmpcount <- countData[,1:12] #only donors 113,114,115
# 
# tmpcoldata <- colData[1:18,]
# tmpcoldata$donor <- as.factor(tmpcoldata$donor)

#NT vs TRIM5
tmpcount <- dplyr::select(countData, matches("TRIM5a|NT"))
tmpcount <- tmpcount[,1:12] #only donors 113,114,115
tmpcoldata <- colData[colData$design!="CPSF6",]
tmpcoldata <- tmpcoldata[1:12,]
tmpcoldata$donor <- as.factor(tmpcoldata$donor)

#creating DESeqDataSet while controlling for effect of donors

dds <- DESeqDataSetFromMatrix(countData = tmpcount,
                              colData = tmpcoldata,
                              design= ~ donor+design ) #controlling for patients

#keeping proteins that have at least a count of 10 for at least 3 treated samples
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


#   #performing variance stabilizing transformation and plotting PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c( "design"))
# ggsave("all3donors_withtrim5_pca.pdf")


dds <- DESeq(dds)
# res <- results(dds, contrast = c("design","CPSF6","NT") )
res <- results(dds, contrast = c("design","TRIM5a","NT") )

# topgenes<- data.matrix(tmpcount[rownames(res[order(res$padj, decreasing = FALSE),])[0:25],])
# logtopgenes <- log(topgenes)
# pheatmap(logtopgenes)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj')
# ggsave("all3donors_trim5vNT_volcanoplot.pdf")


dds_results_filtered <-res[complete.cases(res),]
upreg <- dds_results_filtered[(dds_results_filtered$padj < 0.05 & dds_results_filtered$log2FoldChange > 0),]
downreg <- dds_results_filtered[(dds_results_filtered$padj < 0.05 & dds_results_filtered$log2FoldChange < 0),]
# diffreg <- dds_results_filtered[(dds_results_filtered$padj < 0.05 & abs(dds_results_filtered$log2FoldChange) > 0),]

write.csv(upreg, "all3donors_trim5vNT_upreg.csv")
# write.csv(diffreg, "all3donors_cpsf6vNT_diffreg.csv")
write.csv(downreg, "all3donors_trim5vNT_downreg.csv")



