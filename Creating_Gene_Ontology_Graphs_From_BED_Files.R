library('installr')
library(tidyverse)
library(ChIPQC)
library(ChIPseeker)
library(DiffBind)
library(clusterProfiler)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)

setwd("C:/Users/cathe/Desktop/NordenLab/Project/07102020 Trying GO Graphs with OSTN Motif/WorkingDirectory")
samplefiles <- list.files("data", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("DAY12", "DAY3", "DAY6", "DAY9", "iPSC", "MEF")
samplefiles <- samplefiles[c("MEF", "DAY3", "DAY6", "DAY9","DAY12", "iPSC")]
samplefiles
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)
peakAnnoList
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")


nanog_annot <- data.frame(peakAnnoList[["DAY3"]]@anno)
entrez <- nanog_annot$geneId
ego <- enrichGO(gene = entrez, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "none", pvalueCutoff = 0.002, qvalueCutoff = 1, readable = TRUE)
dotplot(ego, showCategory=20) + ggtitle("Day3 OSTN Motif GO Terms ")

nanog_annot6 <- data.frame(peakAnnoList[["DAY6"]]@anno)
entrez6 <- nanog_annot6$geneId
ego6 <- enrichGO(gene = entrez6, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "none", pvalueCutoff = 0.002, qvalueCutoff = 1, readable = TRUE)
dotplot(ego6, showCategory=20) + ggtitle("Day6 OSTN Motif GO Terms ")

nanog_annot9 <- data.frame(peakAnnoList[["DAY9"]]@anno)
entrez9 <- nanog_annot9$geneId
ego9 <- enrichGO(gene = entrez9, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "none", pvalueCutoff = 0.002, qvalueCutoff = 1, readable = TRUE)
dotplot(ego9, showCategory=15) + ggtitle("Day9 OSTN Motif GO Terms ")

nanog_annot12 <- data.frame(peakAnnoList[["DAY12"]]@anno)
entrez12 <- nanog_annot12$geneId
ego12 <- enrichGO(gene = entrez12, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "none", pvalueCutoff = 0.002, qvalueCutoff = 1, readable = TRUE)
dotplot(ego12, showCategory=20) + ggtitle("Day12 OSTN Motif GO Terms ")

nanog_annotIPSC <- data.frame(peakAnnoList[["iPSC"]]@anno)
entrezIPSC <- nanog_annotIPSC$geneId
egoIPSC <- enrichGO(gene = entrezIPSC, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "none", pvalueCutoff = 0.002, qvalueCutoff = 1, readable = TRUE)
dotplot(egoIPSC, showCategory=20) + ggtitle("IPSC OSTN Motif GO Terms ")

nanog_annotMEF <- data.frame(peakAnnoList[["MEF"]]@anno)
entrezMEF <- nanog_annotMEF$geneId
egoMEF <- enrichGO(gene = entrezMEF, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "none", pvalueCutoff = 0.002, qvalueCutoff = 1, readable = TRUE)
dotplot(egoMEF, showCategory=20) + ggtitle("MEF OSTN Motif GO Terms ")

ekegg <- enrichKEGG(gene = entrez, organism = 'mouse', pvalueCutoff = 0.2)
dotplot(ekegg)
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
compKEGG <- compareCluster(geneCluster = genes, fun = "enrichKEGG", organism = "mouse", pvalueCutoff  = 0.05, pAdjustMethod = "none")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")

compKEGG <- compareCluster(geneCluster = genes, fun = "enrichGO", organism = "mouse", pvalueCutoff  = 0.05, pAdjustMethod = "none")
dotplot(compKEGG, showCategory = 20, title = "GO Pathway Enrichment Analysis")