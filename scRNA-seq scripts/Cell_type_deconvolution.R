##Load packages
library(clustermole)
library(SCpubr)
library(Seurat)
library(UCell)


##Set working directory
setwd("~/Documents/TomoSeq_project/Reference_single_cell/")

##Load of replicate experiments, merge data, and remove individual data files from R environment
Data <- readRDS("Mouse_adult_gingiva.rds")

rownames(Data@assays[["RNA"]]@counts) <- Data@assays[["RNA"]]@meta.features$feature_name
rownames(Data@assays[["RNA"]]@data) <- Data@assays[["RNA"]]@meta.features$feature_name
rownames(Data@assays[["RNA"]]@meta.features) <- Data@assays[["RNA"]]@meta.features$feature_name

Data <- NormalizeData(Data)

##Cell cycle analysis
mouse_cell_cycle_genes <- read.csv("cellcyclegenes.csv", header = T)

s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

Data.cc <- CellCycleScoring(Data, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

head(x = Data.cc@meta.data)

Data$S.Score <- Data.cc$S.Score
Data$G2M.Score <- Data.cc$G2M.Score
Data$Phase <- Data.cc$Phase
head(x = Data@meta.data)

table(Data$Phase)

Data$CC.Difference <- Data$S.Score - Data$G2M.Score

rm(Data.cc)

Data <- FindVariableFeatures(Data)
Data <- ScaleData(Data, vars.to.regress = c("nFeaturess_RNA","nCounts_RNA","pct_counts_mito","pct_counts_ribo","CC.Difference"))
Data <- RunPCA(Data)

ElbowPlot(Data, ndims = 50, reduction = "pca")

Data <- FindNeighbors(Data, reduction = "pca", dims = 1:15)
Data <- FindClusters(Data, resolution = 0.5, algorithm = 1)
Data <- RunUMAP(Data, reduction = "pca", dims = 1:15, reduction.name = "umap")

pdf("ClusterUMAP.pdf", width = 5, height = 5)
SCpubr::do_DimPlot(Data, reduction = "umap", group.by = "seurat_clusters", label = TRUE, legend.position = "none")
dev.off()

pdf("ClusterUMAP_AuthorCellTypes.pdf", width = 8, height = 5)
SCpubr::do_DimPlot(Data, reduction = "umap", group.by = "author_cell_type", label = F, legend.position = "right")
dev.off()

genes <- list("EPI" = c("Cdh1","Epcam"),
              "Keratinocyte" = c("Krt5","Krt14","Krtdap","Krt4","Krt13"),
              "Salivary" = c("Krt8","Krt18"),
              "MES" = c("Vim"),
              "Fibroblasts" = c("Pdgfra","Col1a1","Lum"),
              "Mural" = c("Acta2","Myh11"),
              "BECs" = c("Cdh5","Pecam1"),
              "LECs" = c("Lyve1","Prox1"), 
              "Muscle" = c("Des","Tnnt1"), 
              "Glial" = c("Ncam1","Mpz"), 
              "IMM" = c("Ptprc"),
              "Myeloid" = c("Cd74","Mrc1"),
              "Lymphoid" = c("Cd3e","Cd5"))

pdf("LineageMarkerGenes_DotPlot_AuthorCellTypes.pdf", width = 16, height = 5)
DotPlot(Data, group.by = "author_cell_type", features = genes, cluster.idents = F, cols = c("#4575b4","#d73027")) + scale_colour_gradient2(low="#4575b4", mid="white", high="#d73027") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + RotatedAxis()
dev.off()

pdf("LineageMarkerGenes_DotPlot_res05.pdf", width = 16, height = 10)
DotPlot(Data, group.by = "seurat_clusters", features = genes, cluster.idents = T, cols = c("#4575b4","#d73027")) + RotatedAxis()
dev.off()

##Find markers for every cluster compared to all remaining cells, report only the positive ones
Data.markers <- FindAllMarkers(object = Data, only.pos = TRUE, min.pct = 0.25, 
                                  thresh.use = 0.25, test.use = "wilcox")
write.csv(Data.markers,"Cluster_markers_Wilcox_res05_AllCells.csv")

pdf("TOP5_genes_per_cluster_res05.pdf", width = 10, height = 15)
Data.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(Data, features = top5$gene, size = 3, angle = 0) + NoLegend() + scale_fill_gradientn(colors = viridis(10))
dev.off()

##Cell annotation
agg_exp_mat <- AggregateExpression(Data, assay = "RNA", slot = "counts", verbose = TRUE)

agg_exp_mat <- as.matrix(agg_exp_mat$RNA)
agg_exp_mat <- log1p(agg_exp_mat)

enrich_tblagg <- clustermole_enrichment(expr_mat = agg_exp_mat, species = "mm")

write.csv(enrich_tblagg,"CellTypeAnnotationAllcells_res05_Aggregate.csv")

Idents(Data)

Data$cellclusters <- Idents(Data)

Idents(Data) <- Data$cellclusters

Idents(Data)

celltypes <- c("Fibroblasts","Fibroblasts","BECs","BECs","Keratinocytes","LECs","Fibroblasts","Mural","Glial","Myeloid","Salivary","Muscle","Glial","Myeloid","Fibroblasts","BECs","BECs","Lymphoid","Muscle")

names(celltypes) <- levels(Data)

Data <- RenameIdents(Data, celltypes)

Data$celltypes <- Idents(Data)

table(Data$celltypes)

pdf("ClusterUMAP_CellTypes.pdf", width = 7, height = 5)
SCpubr::do_DimPlot(Data, reduction = "umap", group.by = "celltypes", label = F, legend.position = "right")
dev.off()

pdf("CellTypesUMAP.pdf", width = 10, height = 10)
SCpubr::do_DimPlot(Data, reduction = "umap", group.by = "celltypes", label = F, legend.position = "bottom", shuffle = F)
dev.off()

Idents(Data) <- Data$celltypes

my_levels <- c("Keratinocytes","Salivary","Fibroblasts", "Mural", "BECs", "LECs", "Muscle", "Glial", "Myeloid", "Lymphoid")

Idents(Data) <- factor(Idents(Data), levels= my_levels)

Idents(Data)

pdf("LineageMarkerGenes_DotPlot_CellTypes.pdf", width = 16, height = 4)
DotPlot(Data, features = genes, cluster = F) + RotatedAxis() + scale_colour_gradient2(low="#4575b4", mid="white", high="#d73027") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + RotatedAxis()
dev.off()


HP_cluster_genes <- read.csv ("HP_cluster_marker.csv", header = TRUE)

HPmarkers <- list()
HPmarkers$cl0 <- HP_cluster_genes$cl0[1:100]
HPmarkers$cl1 <- HP_cluster_genes$cl1[1:100]
HPmarkers$cl2 <- HP_cluster_genes$cl2[1:100]
HPmarkers$cl3 <- HP_cluster_genes$cl3[1:100]
HPmarkers$cl4 <- HP_cluster_genes$cl4[1:100]

Data <- AddModuleScore_UCell(Data, features = HPmarkers)
HPsignature.names <- paste0(names(HPmarkers), "_UCell")

pdf("UCell_DotPlot_HP.pdf", width = 5, height = 3.5)
DotPlot(Data, features = HPsignature.names, group.by = "celltypes") + scale_colour_gradient2(low="#4575b4", mid="white", high="#d73027") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + RotatedAxis()
dev.off()



BM_cluster_genes <- read.csv ("BM_cluster_marker.csv", header = TRUE)

BMmarkers <- list()
BMmarkers$cl0 <- BM_cluster_genes$cl0[1:100]
BMmarkers$cl1 <- BM_cluster_genes$cl1[1:100]
BMmarkers$cl2 <- BM_cluster_genes$cl2[1:100]
BMmarkers$cl3 <- BM_cluster_genes$cl3[1:100]
BMmarkers$cl4 <- BM_cluster_genes$cl4[1:100]

Data <- AddModuleScore_UCell(Data, features = BMmarkers)
BMsignature.names <- paste0(names(BMmarkers), "_UCell")

pdf("UCell_DotPlot_BM.pdf", width = 5, height = 3.5)
DotPlot(Data, features = BMsignature.names, group.by = "celltypes") + scale_colour_gradient2(low="#4575b4", mid="white", high="#d73027") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + RotatedAxis()
dev.off()

VT_cluster_genes <- read.csv ("VT_cluster_marker.csv", header = TRUE)

VTmarkers <- list()
VTmarkers$cl0 <- VT_cluster_genes$cl0[1:100]
VTmarkers$cl1 <- VT_cluster_genes$cl1[1:100]
VTmarkers$cl2 <- VT_cluster_genes$cl2[1:100]
VTmarkers$cl3 <- VT_cluster_genes$cl3[1:100]
VTmarkers$cl4 <- VT_cluster_genes$cl4[1:100]
VTmarkers$cl5 <- VT_cluster_genes$cl5[1:100]
VTmarkers$cl6 <- VT_cluster_genes$cl6[1:100]
VTmarkers$cl7 <- VT_cluster_genes$cl7[1:100]
VTmarkers$cl8 <- VT_cluster_genes$cl8[1:100]

Data <- AddModuleScore_UCell(Data, features = VTmarkers)
VTsignature.names <- paste0(names(VTmarkers), "_UCell")

pdf("UCell_DotPlot_VT.pdf", width = 6, height = 3.5)
DotPlot(Data, features = VTsignature.names, group.by = "celltypes") + scale_colour_gradient2(low="#4575b4", mid="white", high="#d73027") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + RotatedAxis()
dev.off()


DT_cluster_genes <- read.csv ("DT_cluster_marker.csv", header = TRUE)

DTmarkers <- list()
DTmarkers$cl0 <- DT_cluster_genes$cl0[1:100]
DTmarkers$cl1 <- DT_cluster_genes$cl1[1:100]
DTmarkers$cl2 <- DT_cluster_genes$cl2[1:100]
DTmarkers$cl3 <- DT_cluster_genes$cl3[1:100]
DTmarkers$cl4 <- DT_cluster_genes$cl4[1:100]
DTmarkers$cl5 <- DT_cluster_genes$cl5[1:100]
DTmarkers$cl6 <- DT_cluster_genes$cl6[1:100]
DTmarkers$cl7 <- DT_cluster_genes$cl7[1:100]

Data <- AddModuleScore_UCell(Data, features = DTmarkers)
DTsignature.names <- paste0(names(DTmarkers), "_UCell")


pdf("UCell_DotPlot_DT.pdf", width = 5.5, height = 3.5)
DotPlot(Data, features = DTsignature.names, group.by = "celltypes") + scale_colour_gradient2(low="#4575b4", mid="white", high="#d73027") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + RotatedAxis()
dev.off()

