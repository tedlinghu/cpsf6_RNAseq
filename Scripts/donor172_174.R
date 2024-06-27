library(DESeq2)
library(apeglm)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(ReactomePA)
library(pheatmap)
library(AnnotationDbi)

countData <- read.csv('/projects/b1042/HultquistLab/Daphne/Novogene_NVUS2023080428/gene_count_matrix.csv',row.names = 1, check.names = FALSE)

colData <- read.csv('/projects/b1042/HultquistLab/Daphne/Novogene_NVUS2023080428/condition2.csv',row.names = 1)

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

countData <- countData[,order(colnames(countData))]


colData <- colData[order(rownames(colData)),]

#different conditions that we test for
# tmpcount <- dplyr::select(countData, matches("CPSF6-6|NT-4")) #condition 1

tmpcount <- dplyr::select(countData, matches("CPSF5-3|NT-4")) #condition 2

#tmpcount <- dplyr::select(countData, matches("A77V Sorted|WT Sorted")) #condition 3

#tmpcount <- dplyr::select(countData, matches("N74D Sorted|WT Sorted")) #condition 4

# tmpcount <- dplyr::select(countData, matches("N74D Sorted|Uninfected")) #condition 5

# tmpcount <- dplyr::select(countData, matches("A77V Sorted|Uninfected")) #condition 6

# tmpcount <- dplyr::select(countData, matches("WT Sorted|Uninfected")) #condition 7

#creating DESeqDataSet while controlling for effect of donors
colData$Donor <- as.factor(colData$Donor)
dds <- DESeqDataSetFromMatrix(countData = tmpcount,
                              colData = colData,
                              design= ~ Donor + Treatment) #controlling for donors

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


#   #performing variance stabilizing transformation and plotting PCA

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c( "Treatment"))

#getting results from different conditions
dds <- DESeq(dds)
# res <- results(dds, contrast = c("Treatment","CPSF6-6","NT-4") ) #condition 1

res <- results(dds, contrast = c("Treatment","CPSF5-3","NT-4") ) #condition 2

#res <- results(dds, contrast = c("Treatment","A77V Sorted","WT Sorted") ) #condition 3

#res <- results(dds, contrast = c("Treatment","N74D Sorted","WT Sorted") ) #condition 4

# res <- results(dds, contrast = c("Treatment","N74D Sorted","Uninfected") ) #condition 5

# res <- results(dds, contrast = c("Treatment","A77V Sorted","Uninfected") ) #condition 6

# res <- results(dds, contrast = c("Treatment","WT Sorted","Uninfected") ) #condition 7


# topgenes<- data.matrix(tmpcount[rownames(res[order(res$padj, decreasing = FALSE),])[0:25],])
# logtopgenes <- log(topgenes)
# pheatmap(logtopgenes)
#setEPS()
#postscript("n74d_v_WT.eps") 

#plotting the volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                colAlpha = 1)
#dev.off() 

#filtering for the up or downregulated genes
dds_results_filtered <-res[complete.cases(res),]
# upreg <- dds_results_filtered[dds_results_filtered$log2FoldChange > 0,]
# downreg <- dds_results_filtered[dds_results_filtered$log2FoldChange < 0,]

upreg <- dds_results_filtered[dds_results_filtered$padj < 0.05 & dds_results_filtered$log2FoldChange > 0,]
downreg <- dds_results_filtered[dds_results_filtered$padj < 0.05 & dds_results_filtered$log2FoldChange < 0,]

#matching_rows <- grep("TRIM", rownames(dds_results_filtered), value = TRUE)


write.csv(upreg, "condition2_upreg.csv")
write.csv(downreg, "condition2_downreg.csv")

