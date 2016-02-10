## practice project number one
getwd()

setwd("~/Bioinformatics/myRprojects/")
library(DESeq2)


directory <-"~/Bioinformatics/myRprojects/counts"

sampleFiles <- list.files(directory)
sampleFiles

sampleCondition <- c("P1", "P1","P1","P1","P3","P3", "P3", "P3", "V","V","V","V")
factor(sampleCondition, levels = c("V", "P1", "P3"))



SampleTable <- data.frame(sampleNames = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
SampleTable

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable =  SampleTable,
                                       directory = directory,
                                       design = ~condition)

ddsHTSeq

ddsHTSeq <- DESeq(ddsHTSeq)
res <- results(ddsHTSeq)
res

res1 <- results(ddsHTSeq, contrast = c("condition", "P3", "P1"))
res1; dim(res1)
head(res1) 
head(res1[order(rownames(res1)), ])

mean(rownames(res1) == rownames(res1[order(rownames(res1)), ]))

library(biomaRt) ## to use this get the gene names

mart <- useMart(host="www.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL")
mart2<- useDataset("mmusculus_gene_ensembl", mart)

ensembl_genes<- rownames(res1)

head(ensembl_genes)

checkuptable <- getBM(
    filters= "ensembl_gene_id", values= ensembl_genes,
    attributes= c("ensembl_gene_id",  "external_gene_name"),
    mart= mart2)

dim(checkuptable)
head(checkuptable)
tail(checkuptable)
class(checkuptable)

res1Table <- as.data.frame(res1)
res1Table$ensembl_gene_id <- rownames(res1Table)
head(res1Table)
class(res1Table)

res1Final <- merge(res1Table, checkuptable, by = "ensembl_gene_id")
head(res1Final)
tail(res1Final)


res1FinalOrdered <- res1Final[order(res1Final$padj), ]
head(res1FinalOrdered)

write.csv(res1FinalOrdered, "P3versusP1orderedbyPadj.csv")





getwd()

res2 <- results(ddsHTSeq, contrast = c("condition", "P3", "V"))
res2

res2Table <- as.data.frame(res2)
res2Table$ensembl_gene_id <- rownames(res2Table)
head(res2Table)
class(res2Table)

res2Final <- merge(res2Table, checkuptable, by = "ensembl_gene_id")
head(res2Final)
tail(res2Final)

res2FinalOrdered <- res2Final[order(res2Final$padj), ]
head(res2FinalOrdered)

write.csv(res2FinalOrdered, "P3versusVorderedbyPadj.csv")


##################


res3 <- results(ddsHTSeq, contrast = c("condition", "P1", "V"))
res3
summary(res3)

res3Table <- as.data.frame(res3)
res3Table$ensembl_gene_id <- rownames(res3Table)
head(res3Table)
class(res3Table)

res3Final <- merge(res3Table, checkuptable, by = "ensembl_gene_id")
head(res3Final)
tail(res3Final)

res3FinalOrdered <- res3Final[order(res3Final$padj), ]
head(res3FinalOrdered)

write.csv(res3FinalOrdered, "P1versusVorderedbyPadj.csv")


## to subset 

res3sig <- subset(res3FinalOrdered, padj <0.05)
res3sig
dim(res3sig)


## try out visualization
rld <- rlog(ddsHTSeq, blind=FALSE) ### rlog function provides transformation needed other than differential expression 
head(assay(rld), 3)

##gene clustering
library("genefilter")
library("pheatmap")
library("ggplot2")
library("RColorBrewer")

topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)

plotPCA(rld, intgroup = "condition") ## this works out a PCA graph.

??rowVars


mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("condition", "condition")])

pdf("TryingFeb18.pdf" )
pheatmap(mat) #, annotation_col=df)



tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)

dist1 = dist(mat)
colramp = colorRampPalette(c(3, "white", 2))(9)
heatmap(as.matrix(dist1), col = colramp, cColv = NA, Rowv =  NA)



pheatmap(mat, annotation_col=df)

dev.off()

plotDispEsts(ddsHTSeq)

select <- order(rowMeans(counts(ddsHTSeq,normalized=TRUE)),decreasing=TRUE)[1:20]

nt <- normTransform(ddsHTSeq) # defaults to log2(x+1) ## function not available
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(ddsHTSeq)[,c("condition")])
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE) #, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


plotMA(res1, main = "P3 versus P1", ylim=c(-2,2))

plotMA(res2, main = "P3 versus V", ylim=c(-2,2))

plotMA(res3, main = "P1 versus V", ylim=c(-2,2))

library(biomaRt)
library(biomaRt) ## to use this get the gene names

mart <- useMart(host="www.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL")
useDataset("mmusculus_gene_ensembl", mart)

library(biomaRt)
mart2<- useDataset("mmusculus_gene_ensembl", mart)

ensembl_genes<- rownames(res3sig)

head(ensembl_genes)
checkuptable <- getBM(
    filters= "ensembl_gene_id", values= ensembl_genes,
    attributes= c("ensembl_gene_id",  "entrezgene", "external_gene_name"),
    mart= mart2)

dim(checkuptable)
head(checkuptable)

checkuptable <- getBM(
    filters= "ensembl_gene_id", values= ensembl_genes,
    attributes= c("ensembl_gene_id",  "external_gene_name"),
    mart= mart2)

dim(checkuptable)
class(checkuptable)

rownames(res3sig)

res3sigTable <- as.data.frame(res3sig)
# to make it in order
res3sigTableOrderedbyID <-res3sigTable[order(rownames(res3sigTable)), ]
res3sigTableOrderedbyID$external_gene_name <- checkuptable$external_gene_name

res3sigTableOrderedbyID
write.csv(res3sigTableOrderedbyID, "P1vsVpadjLess05.csv")

## plotting some results
par(mfrow=c(2, 2))
par(mar = c(5, 4, 3, 1))
pdf("topGenes.pdf")
topGene <- rownames(res1)[which.min(res1$padj)]
plotCounts(ddsHTSeq, gene=topGene, intgroup=c("condition"))

counts(ddsHTSeq)[topGene, ]
plot(counts(ddsHTSeq)[topGene, ])

boxplot(counts(ddsHTSeq)[topGene, ]~SampleTable$condition) ### THIS IS GREAT TO KNOW

?boxplot
plotCounts(ddsHTSeq, gene= "What", intgroup=c("condition"), pch = 19, col = "blue")
plotCounts(ddsHTSeq, gene= "What", intgroup=c("condition"), pch = 19, col = "blue")

topGene <- rownames(res2)[which.min(res2$padj)]
plotCounts(ddsHTSeq, gene=topGene, intgroup=c("condition"), pch = 19, col = "blue")

topGene <- rownames(res3)[which.min(res3$padj)]
plotCounts(ddsHTSeq, gene=topGene, intgroup=c("condition"), pch = 19, col = "red")

dev.off()

plotCounts(ddsHTSeq, gene=topGene, intgroup=c("condition"))
plotCounts(ddsHTSeq, gene=topGene, intgroup=c("condition"), pch = 19, col = "blue")


