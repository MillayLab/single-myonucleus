options(future.globals.maxSize = 10000 * 1024^2)
library(Seurat)
library(dplyr)
library(ggplot2)
#fivemonth
fivemonth.data <- Read10X(data.dir = "E:/10X-datasets/10X_Millay_WT_Muscle/outs/mm10/filtered_feature_bc_matrix")
fivemonth_soupx <- CreateSeuratObject(counts = fivemonthcounts, project = "A", min.cells = 3, min.features = 200)
fivemonth_soupx[["percent.mt"]] <- PercentageFeatureSet(fivemonth_soupx, pattern = "^MT-")

VlnPlot(fivemonth_soupx, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
featurescatterplot <- FeatureScatter(fivemonth_soupx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot(featurescatterplot)
fivemonth_soupx <- subset(fivemonth_soupx, subset = nFeature_RNA > 200 & nFeature_RNA < 3200)
fivemonth_soupx <- NormalizeData(fivemonth_soupx)
fivemonth_soupx <- FindVariableFeatures(fivemonth_soupx, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(fivemonth_soupx), 10)
plot1variable <- VariableFeaturePlot(fivemonth_soupx)
plot2variable <- LabelPoints(plot = plot1variable, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1variable, plot2variable))
all.genes5mo <- rownames(fivemonth_soupx)
fivemonth_soupx <- ScaleData(fivemonth_soupx, features = all.genes5mo)
fivemonth_soupx <- RunPCA(fivemonth_soupx, features = VariableFeatures(object = fivemonth_soupx))
VizDimLoadings(fivemonth_soupx, dims = 1:2, reduction = "pca")
DimPlot(fivemonth_soupx, reduction = "pca")
DimHeatmap(fivemonth_soupx, dims = 1:15, cells = 500, balanced = TRUE)
fivemonth_soupx <- FindNeighbors(fivemonth_soupx, dims = 1:12)
fivemonth_soupx <- FindClusters(fivemonth_soupx, resolution = 0.5)
fivemonth_soupx <- RunUMAP(fivemonth_soupx, dims = 1:12)
fivemonth_soupx@reductions[["umap"]] <- fivemonth_reductions
DimPlot(fivemonth_soupx, reduction = "umap", label = TRUE)
fivemonth_soupx <- RenameIdents(fivemonth_soupx, "0" = "Type IIb Myonuclei", "1" = "Type IIx Myonuclei", "2" = "Type IIb Myonuclei #2", "3" = "Type IIx Myonuclei #2", "4" = "FAPs", "5" = "Endothelial Cells", "6" = "Myotendinous Junction", "7" = "Smooth Muscle", "8" = "Satellite Cells", "9" = "Immune Cells", "10" = "Smooth Muscle #2", "11" = "Neuromuscular Junction", "12" = "Tenocytes") 
fivemonth_muscle <- subset(fivemonth_soupx, idents = c("Type IIx Myonuclei", "Type IIb Myonuclei", "Neuromuscular Junction", "Myotendinous Junction", "Type IIb Myonuclei #2", "Type IIx Myonuclei #2"))
fivemonth_soupx$CellType <- Idents(fivemonth_soupx)

#Twentyfour
library(Seurat)
library(dplyr)
twentyfour.data <- Read10X(data.dir = "E:/10X-datasets/10X_Millay_Aging/outs/filtered_feature_bc_matrix")
twentyfour_soupx <- CreateSeuratObject(counts = twentyfourcounts, project = "B", min.cells = 3, min.features = 200)
VlnPlot(twentyfour_soupx, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
featurescatterplot <- FeatureScatter(twentyfour_soupx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot(featurescatterplot)
twentyfour_soupx <-  subset(twentyfour_soupx, subset = nFeature_RNA > 200 & nFeature_RNA < 4000) 
twentyfour_soupx <- NormalizeData(twentyfour_soupx)
twentyfour_soupx <- FindVariableFeatures(twentyfour_soupx, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(twentyfour_soupx), 10)
plot1variable <- VariableFeaturePlot(twentyfour_soupx)
plot2variable <- LabelPoints(plot = plot1variable, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1variable, plot2variable))
all.genes5mo <- rownames(twentyfour_soupx)
twentyfour_soupx <- ScaleData(twentyfour_soupx, features = all.genes5mo)
twentyfour_soupx <- RunPCA(twentyfour_soupx, features = VariableFeatures(object = twentyfour_soupx))
VizDimLoadings(twentyfour_soupx, dims = 1:2, reduction = "pca")
DimPlot(twentyfour_soupx, reduction = "pca")
DimHeatmap(twentyfour_soupx, dims = 1:15, cells = 500, balanced = TRUE)
twentyfour_soupx <- FindNeighbors(twentyfour_soupx, dims = 1:15)
twentyfour_soupx <- FindClusters(twentyfour_soupx, resolution = 0.5)
twentyfour_soupx <- RunUMAP(twentyfour_soupx, dims = 1:15)
twentyfour_soupx@reductions[["umap"]] <- twentyfour_reductions
twentyfour_soupx <- RenameIdents(twentyfour_soupx, "0" = "Type IIb Myonuclei", "1" = "Type IIx Myonuclei", "2" = "Type IIb Myonuclei", "3" = "Type IIx Myonuclei", "4" = "FAPs", "5" = "Smooth Muscle", "6" = "Endothelial Cells", "7" = "Myotendinous Junction", "8" = "Satellite Cells", "9" = "Immune Cells", "10" = "Smooth Muscle", "11" = "Endothelial Cells", "12" = "Neuromuscular Junction", "13" = "Schwann Cells", "14" = "Unidentified")
twentyfour_soupx$CellType <- Idents(twentyfour_soupx)
DimPlot(twentyfour_soupx, reduction = "umap", label = TRUE)
twentyfour_muscle <- subset(twentyfour_soupx, idents = c("Type IIx Myonuclei", "Type IIb Myonuclei", "Neuromuscular Junction", "Myotendinous Junction"))

#thirty
library(Seurat)
library(dplyr)
thirty.data <- Read10X(data.dir = "E:/10X-datasets/10X-Millay-30-Month/filtered_feature_bc_matrix")
thirty_soupx <- CreateSeuratObject(counts = thirtycounts, project = "C", min.cells = 3, min.features = 200)
VlnPlot(thirty_soupx, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
featurescatterplot <- FeatureScatter(thirty_soupx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot(featurescatterplot)
thirty_soupx <- subset(thirty_soupx, subset = nFeature_RNA > 200 & nFeature_RNA < 4000)
thirty_soupx <- NormalizeData(thirty_soupx)
thirty_soupx <- FindVariableFeatures(thirty_soupx, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(thirty_soupx), 10)
plot1variable <- VariableFeaturePlot(thirty_soupx)
plot2variable <- LabelPoints(plot = plot1variable, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1variable, plot2variable))
all.genes5mo <- rownames(thirty_soupx)
thirty_soupx <- ScaleData(thirty_soupx, features = all.genes5mo)
thirty_soupx <- RunPCA(thirty_soupx, features = VariableFeatures(object = thirty_soupx))
VizDimLoadings(thirty_soupx, dims = 1:2, reduction = "pca")
DimPlot(thirty_soupx, reduction = "pca")
DimHeatmap(thirty_soupx, dims = 1:15, cells = 500, balanced = TRUE)
thirty_soupx <- FindNeighbors(thirty_soupx, dims = 1:20)
thirty_soupx <- FindClusters(thirty_soupx, resolution = 0.5)
thirty_soupx <- RunUMAP(thirty_soupx, dims = 1:20)
thirty_soupx@redctions[["umap"]] <- thirty_reductions
DimPlot(thirty_soupx, reduction = "umap", label = TRUE)
thirty_soupx <- RenameIdents(thirty_soupx, "0" = "Type IIb Myonuclei", "1" = "Type IIx Myonuclei", "2" = "Nr4a3+ Myonuclei", "3" = "FAPs", "4" = "Immune Cells", "5" = "Smooth Muscle", "6" = "Unnamed", "7" = "Enah+ Myonuclei", "8" = "Endothelial Cells", "9" = "Myotendinous Junction", "10" = "Ampd3+ Myonuclei", "11" = "Neuromuscular Junction", "12" = "Satellite Cells", "13" = "Tenocytes", "14" = "Smooth Muscle #2")
thirty_muscle <- subset(thirty_soupx, idents = c("Type IIx Myonuclei", "Type IIb Myonuclei", "Nr4a3+ Myonuclei", "Myotendinous Junction", "Ampd3+ Myonuclei", "Unnamed", "Enah+ Myonuclei", "Neuromuscular Junction"))


fivemonth_soupx@meta.data[["orig.ident"]] <- "A"
twentyfour_soupx@meta.data[["orig.ident"]] <- "B"
thirty_soupx@meta.data[["orig.ident"]] <- "C"
integration.list <- list(thirty_soupx, fivemonth_soupx, twentyfour_soupx)
for (i in 1:length(integration.list)) {
  integration.list[[i]] <- NormalizeData(integration.list[[i]], verbose = FALSE)
  integration.list[[i]] <- FindVariableFeatures(integration.list[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
}


integration.anchors <- FindIntegrationAnchors(object.list = integration.list, dims = 1:30, verbose = FALSE)
integration.integrated <- IntegrateData(anchorset = integration.anchors, dims = 1:30)
integration.integrated <- ScaleData(integration.integrated, verbose = FALSE)
integration.integrated <- RunPCA(integration.integrated, npcs = 30, verbose = FALSE)
integration.integrated <- RunUMAP(integration.integrated, reduction = "pca", dims = 1:30)
DimPlot(integration.integrated, reduction = "umap", label = TRUE)
DimPlot(integration.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE) + NoLegend()


thirty_muscle <- subset(thirty_soupx, idents = c("Type IIx Myonuclei", "Type IIb Myonuclei", "Nr4a3+ Myonuclei", "Myotendinous Junction", "Ampd3+ Myonuclei", "Unnamed", "Enah+ Myonuclei", "Neuromuscular Junction"))
twentyfour_muscle <- subset(twentyfour_soupx, idents = c("Type IIx Myonuclei", "Type IIb Myonuclei", "Neuromuscular Junction", "Musculotendinous Junction"))
fivemonth_muscle <- subset(fivemonth_soupx, idents = c("Type IIx Myonuclei", "Type IIb Myonuclei", "Neuromuscular Junction", "Myotendinous Junction", "Type IIb Myonuclei #2", "Type IIx Myonuclei #2"))

fivemonth_muscle@meta.data[["orig.ident"]] <- "A"
twentyfour_muscle@meta.data[["orig.ident"]] <- "B"
thirty_muscle@meta.data[["orig.ident"]] <- "C"

integration.list.muscle <- list(thirty_muscle, fivemonth_muscle, twentyfour_muscle)
for (i in 1:length(integration.list)) {
  integration.list.muscle[[i]] <- NormalizeData(integration.list.muscle[[i]], verbose = FALSE)
  integration.list.muscle[[i]] <- FindVariableFeatures(integration.list.muscle[[i]], selection.method = "vst", 
                                                       nfeatures = 2000, verbose = FALSE)
}


integration.anchors.muscle <- FindIntegrationAnchors(object.list = integration.list.muscle, dims = 1:30, verbose = FALSE)
integration.integrated.muscle <- IntegrateData(anchorset = integration.anchors.muscle, dims = 1:30)
DefaultAssay(integration.integrated.muscle) <- "Integrated"
integration.integrated.muscle <- ScaleData(integration.integrated.muscle, verbose = FALSE)
integration.integrated.muscle <- RunPCA(integration.integrated.muscle, npcs = 30, verbose = FALSE)
integration.integrated.muscle <- RunUMAP(integration.integrated.muscle, reduction = "pca", dims = 1:30)
DimPlot(integration.integrated.muscle, reduction = "umap", label = TRUE)
DimPlot(integration.integrated.muscle, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE) + NoLegend()
DefaultAssay(integration.integrated.muscle) <- "RNA"
integration.integrated.muscle <- ScaleData(integration.integrated.muscle, verbose = FALSE)

DoHeatmap(integration.integrated.muscle, features = c("Nr4a3", "Hk2", "Smad3", "Ankrd33b", "Tmem140", "Fam134b", "Irs2", "Map3k14", "Actr3b", "Gramd4", "Gm42418", "Mybpc1", "Ppp1r3a", "Zbtb20", "Hnrnpa2b1"), group.by = "orig.ident") + NoLegend()

A_markers <- FindMarkers(integration.integrated.muscle, ident.1 = "A", ident.2 = NULL, group.by = "orig.ident", only.pos = TRUE)
B_markers <- FindMarkers(integration.integrated.muscle, ident.1 = "B", ident.2 = "A", group.by = "orig.ident", only.pos = TRUE)
C_markers <- FindMarkers(integration.integrated.muscle, ident.1 = "C", ident.2 = "A", group.by = "orig.ident", only.pos = TRUE)


write.xlsx(A_markers, "Z:/Michael/snRNA-seq Revisions/5month myonuclei positive markers.xlsx")
write.xlsx(B_markers, "Z:/Michael/snRNA-seq Revisions/24month positive markers.xlsx")
write.xlsx(C_markers, "Z:/Michael/snRNA-seq Revisions/30month positive markers.xlsx")
