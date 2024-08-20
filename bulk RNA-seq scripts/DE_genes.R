#Set working directory
setwd("~/Documents/TomoSeq_project/DESeq2/")

library("DESeq2")

#Load meta and expression data
metadata <- read.csv(paste0("Tomobulk_metadata.csv"), header =T)
metadata$Tissue <- factor(metadata$Tissue)
rownames(metadata) <- metadata$Section

expressiondata <- read.csv("Tomobulk_expressiondata.csv", header = T)
rownames(expressiondata) <- expressiondata$X
expressiondata <- expressiondata[,-1]


#Make DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = expressiondata,
  colData = metadata,
  design = ~ Tissue + Batch
)


#Compute differentially expressed genes per oral site
HPTvsAllT <- DESeqDataSetFromMatrix(countData = expressiondata,colData = metadata, design = ~ HPvsAll)
BMTvsAllT <- DESeqDataSetFromMatrix(countData = expressiondata,colData = metadata, design = ~ BMvsAll)
VTTvsAllT <- DESeqDataSetFromMatrix(countData = expressiondata,colData = metadata, design = ~ VTvsAll)
DTTvsAllT <- DESeqDataSetFromMatrix(countData = expressiondata,colData = metadata, design = ~ DTvsAll)

mcols(dds) <- cbind(mcols(dds))
comparisons = list (dds, HPTvsAllT, BMTvsAllT, VTTvsAllT, DTTvsAllT)
vsdlist = list ()

for (i in 1:length(comparisons) ) {
  comparisons[[i]] <- estimateSizeFactors(comparisons[[i]])
  vsdlist[[i]] = vst(comparisons[[i]])
  comparisons[[i]] <- DESeq(comparisons[[i]])
}

#Extract differentially expressed genes per oral site
HPTvsAllT_res <- results(comparisons[[2]], contrast=c("HPvsAll","HP","All"), alpha = 0.05) 
BMTvsAllT_res <- results(comparisons[[3]], contrast=c("BMvsAll","BM","All"), alpha = 0.05)
VTTvsAllT_res <- results(comparisons[[4]], contrast=c("VTvsAll","VT","All"), alpha = 0.05) 
DTTvsAllT_res <- results(comparisons[[5]], contrast=c("DTvsAll","DT","All"), alpha = 0.05) 

HPTvsAllT_resSig <- subset(HPTvsAllT_res, padj < 0.05)
HPTvsAllT_resSig <- HPTvsAllT_resSig[ order(HPTvsAllT_resSig$log2FoldChange, decreasing = TRUE), ]

BMTvsAllT_resSig <- subset(BMTvsAllT_res, padj < 0.05)
BMTvsAllT_resSig <- BMTvsAllT_resSig[ order(BMTvsAllT_resSig$log2FoldChange, decreasing = TRUE), ]

VTTvsAllT_resSig <- subset(VTTvsAllT_res, padj < 0.05)
VTTvsAllT_resSig <- VTTvsAllT_resSig[ order(VTTvsAllT_resSig$log2FoldChange, decreasing = TRUE), ]

DTTvsAllT_resSig <- subset(DTTvsAllT_res, padj < 0.05)
DTTvsAllT_resSig <- DTTvsAllT_resSig[ order(DTTvsAllT_resSig$log2FoldChange, decreasing = TRUE), ]

HPTvsAllT_resSigDF <- as.data.frame(HPTvsAllT_resSig)%>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange))
BMTvsAllT_resSigDF <- as.data.frame(BMTvsAllT_resSig)%>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange))
VTTvsAllT_resSigDF <- as.data.frame(VTTvsAllT_resSig)%>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange))
DTTvsAllT_resSigDF <- as.data.frame(DTTvsAllT_resSig)%>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange))


write.csv(HPTvsAllT_resSigDF,file="HPTvsAllT_resSigDF.csv",row.names=TRUE)
write.csv(BMTvsAllT_resSigDF,file="BMTvsAllT_resSigDF.csv",row.names=TRUE)
write.csv(VTTvsAllT_resSigDF,file="VTTvsAllT_resSigDF.csv",row.names=TRUE)
write.csv(DTTvsAllT_resSigDF,file="DTTvsAllT_resSigDF.csv",row.names=TRUE)


#HP spatially reproducible genes
reprGenesHP <- read.csv(file="corr_genes_HP.tsv", sep = "\t")
reprGenesHP <- as.data.frame(reprGenesHP)%>% filter(adj.pv < 0.05)

genesHP <- reprGenesHP$gene
HPT <- as.data.frame(HPTvsAllT_res)
HPT$gene <- rownames(HPT)
rownames(reprGenesHP) <- reprGenesHP$gene
reprGenesHP <- as.data.frame(reprGenesHP)

HPgenes <- merge(HPT,reprGenesHP, by = "gene", all = TRUE)
HPgenes <- HPgenes[HPgenes$gene %in% reprGenesHP$gene,]

HPgenesSigDF <- as.data.frame(HPgenes) %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange))
HPgenesSigUP <- as.data.frame(HPgenes) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange))
HPgenesSigDOWN <- as.data.frame(HPgenes) %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0) %>% arrange(desc(log2FoldChange))

write.csv(HPgenesSigDF,file="reprGenesHP_HPvsAllgenesSigDF.csv",row.names=TRUE)
write.csv(HPgenesSigUP,file="reprGenesHP_HPvsAllgenesSigUP.csv",row.names=TRUE)
write.csv(HPgenesSigDOWN,file="reprGenesHP_HPvsAllgenesSigDOW.csv",row.names=TRUE)


#BM spatially reproducible genes
reprGenesBM <- read.csv(file="corr_genes_BM.tsv", sep = "\t")
reprGenesBM <- as.data.frame(reprGenesBM)%>% filter(adj.pv < 0.05)

genesBM <- reprGenesBM$gene
BMT <- as.data.frame(BMTvsAllT_res)
BMT$gene <- rownames(BMT)
rownames(reprGenesBM) <- reprGenesBM$gene
reprGenesBM <- as.data.frame(reprGenesBM)

BMgenes <- merge(BMT,reprGenesBM, by = "gene", all = TRUE)
BMgenes <- BMgenes[BMgenes$gene %in% reprGenesBM$gene,]

BMgenesSigDF <- as.data.frame(BMgenes) %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange))
BMgenesSigUP <- as.data.frame(BMgenes) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange))
BMgenesSigDOWN <- as.data.frame(BMgenes) %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0) %>% arrange(desc(log2FoldChange))

write.csv(BMgenesSigDF,file="reprGenesBM_BMvsAllgenesSigDF.csv",row.names=TRUE)
write.csv(BMgenesSigUP,file="reprGenesBM_BMvsAllgenesSigUP.csv",row.names=TRUE)
write.csv(BMgenesSigDOWN,file="reprGenesBM_BMvsAllgenesSigDOW.csv",row.names=TRUE)


#VT spatially reproducible genes
reprGenesVT <- read.csv(file="corr_genes_VT.tsv", sep = "\t")
reprGenesVT <- as.data.frame(reprGenesVT)%>% filter(adj.pv < 0.05)

genesVT <- reprGenesVT$gene

VTT <- as.data.frame(VTTvsAllT_res)
VTT$gene <- rownames(VTT)
rownames(reprGenesVT) <- reprGenesVT$gene
reprGenesVT <- as.data.frame(reprGenesVT)

VTgenes <- merge(VTT,reprGenesVT, by = "gene", all = TRUE)
VTgenes <- VTgenes[VTgenes$gene %in% reprGenesVT$gene,]

VTgenesSigDF <- as.data.frame(VTgenes) %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange))
VTgenesSigUP <- as.data.frame(VTgenes) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange))
VTgenesSigDOWN <- as.data.frame(VTgenes) %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0) %>% arrange(desc(log2FoldChange))

write.csv(VTgenesSigDF,file="reprGenesVT_VTvsAllgenesSigDF.csv",row.names=TRUE)
write.csv(VTgenesSigUP,file="reprGenesVT_VTvsAllgenesSigUP.csv",row.names=TRUE)
write.csv(VTgenesSigDOWN,file="reprGenesVT_VTvsAllgenesSigDOW.csv",row.names=TRUE)


#DT spatially reproducible genes
reprGenesDT <- read.csv(file="corr_genes_DT.tsv", sep = "\t")
reprGenesDT <- as.data.frame(reprGenesDT)%>% filter(adj.pv < 0.05)

genesDT <- reprGenesDT$gene

DTT <- as.data.frame(DTTvsAllT_res)
DTT$gene <- rownames(DTT)
rownames(reprGenesDT) <- reprGenesDT$gene
reprGenesDT <- as.data.frame(reprGenesDT)

DTgenes <- merge(DTT,reprGenesDT, by = "gene", all = TRUE)
DTgenes <- DTgenes[DTgenes$gene %in% reprGenesDT$gene,]

DTgenesSigDF <- as.data.frame(DTgenes) %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange))
DTgenesSigUP <- as.data.frame(DTgenes) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange))
DTgenesSigDOWN <- as.data.frame(DTgenes) %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0) %>% arrange(desc(log2FoldChange))

write.csv(DTgenesSigDF,file="reprGenesDT_DTvsAllgenesSigDF.csv",row.names=TRUE)
write.csv(DTgenesSigUP,file="reprGenesDT_DTvsAllgenesSigUP.csv",row.names=TRUE)
write.csv(DTgenesSigDOWN,file="reprGenesDT_DTvsAllgenesSigDOW.csv",row.names=TRUE)

