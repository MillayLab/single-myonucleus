p10_soupx <- CreateSeuratObject(counts = p10_counts, project = "p10", min.cells = 3, min.features = 200)
VlnPlot(p10_soupx, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
featurescatterplot <- FeatureScatter(p10_soupx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot(featurescatterplot)
p10_soupx <- subset(p10_soupx, subset = nFeature_RNA > 200 & nFeature_RNA < 4200)
p10_soupx <- NormalizeData(p10_soupx)
p10_soupx <- FindVariableFeatures(p10_soupx, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(p10_soupx), 10)
plot1variable <- VariableFeaturePlot(p10_soupx)
plot2variable <- LabelPoints(plot = plot1variable, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1variable, plot2variable))
all.genes5mo <- rownames(p10_soupx)
p10_soupx <- ScaleData(p10_soupx, features = all.genes5mo)
p10_soupx <- RunPCA(p10_soupx, features = VariableFeatures(object = p10_soupx))
VizDimLoadings(p10_soupx, dims = 1:2, reduction = "pca")
DimPlot(p10_soupx, reduction = "pca")
DimHeatmap(p10_soupx, dims = 1:15, cells = 500, balanced = TRUE)
p10_soupx <- FindNeighbors(p10_soupx, dims = 1:20)
p10_soupx <- FindClusters(p10_soupx, resolution = 0.5)
p10_soupx <- RunUMAP(p10_soupx, dims = 1:20)
p10_soupx@reductions[["umap"]] <- p10_reductions
DimPlot(p10_soupx, reduction = "umap", label = TRUE)
p10_soupx <- subset(p10_soupx, idents = c("6"), invert = TRUE)
set.seed(42)
p10_soupx <- FindNeighbors(p10_soupx, dims = 1:20)
p10_soupx <- FindClusters(p10_soupx, resolution = 0.5)
p10_soupx <- RunUMAP(p10_soupx, dims = 1:20) 
DimPlot(p10_soupx, reduction = "umap", label = TRUE)

p10_soupx <- RenameIdents(p10_soupx, "0" = "Type IIb Myonuclei", "1" = "Type IIx Myonuclei", "2" = "FAPs", "3" = "Endothelial Cells", "4" = "Tenocytes", "5" = "Smooth Muscle", "6" = "Satellite Cells", "7" = "Immune Cells", "8"= "Neuromuscular Junction", "9" = "Schwann Cells", "10" = "Myocytes", "11" = "Myotendinous Junction", "12" = "Type I Myonuclei")
DimPlot(p10_soupx, reduction = "umap", label = TRUE)
