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
#supplying dimensions for consistent graphs
fivemonth_soupx@reductions[["umap"]] <- fivemonth_reductions
DimPlot(fivemonth_soupx, reduction = "umap", label = TRUE)
fivemonth_soupx <- RenameIdents(fivemonth_soupx, "0" = "Type IIb Myonuclei", "1" = "Type IIx Myonuclei", "2" = "Type IIb Myonuclei #2", "3" = "Type IIx Myonuclei #2", "4" = "FAPs", "5" = "Endothelial Cells", "6" = "Myotendinous Junction", "7" = "Smooth Muscle", "8" = "Satellite Cells", "9" = "Immune Cells", "10" = "Smooth Muscle #2", "11" = "Neuromuscular Junction", "12" = "Tenocytes") 


#feature expression plots in figure 1
FeaturePlot(fivemonth_soupx, features = c("Myh4"), pt.size = 2.5, cols = c("lightgrey", "red")) + NoAxes() + NoLegend() + ggtitle("")
FeaturePlot(fivemonth_soupx, features = c("Myh1"), pt.size = 2.5, cols = c("lightgrey", "red")) + NoAxes() + NoLegend() + ggtitle("")
FeaturePlot(fivemonth_soupx, features = c("Chrne"), pt.size = 2.5, cols = c("lightgrey", "red")) + NoAxes() + NoLegend() + ggtitle("")
FeaturePlot(fivemonth_soupx, features = c("Col22a1"), pt.size = 2.5, cols = c("lightgrey", "red")) + NoAxes() + NoLegend() + ggtitle("")
FeaturePlot(fivemonth_soupx, features = c("Pax7"), pt.size = 2.5, cols = c("lightgrey", "red")) + NoAxes() + NoLegend() + ggtitle("")


fivemonth_soupx <- RenameIdents(fivemonth_soupx, "Type IIb Myonuclei #2" = "Type IIb Myonuclei", "Type IIx Myonuclei #2" = "Type IIx Myonuclei") 
VlnPlot(fivemonth_soupx, features = c("Myh4"), idents = c("Type IIb Myonuclei", "Type IIx Myonuclei", "Neuromuscular Junction", "Myotendinous Junction", "Satellite Cells"), pt.size = 0) + NoLegend()
VlnPlot(fivemonth_soupx, features = c("Myh1"), idents = c("Type IIb Myonuclei", "Type IIx Myonuclei", "Neuromuscular Junction", "Myotendinous Junction", "Satellite Cells"), pt.size = 0) + NoLegend()
VlnPlot(fivemonth_soupx, features = c("Chrne"), idents = c("Type IIb Myonuclei", "Type IIx Myonuclei", "Neuromuscular Junction", "Myotendinous Junction", "Satellite Cells"), pt.size = 0) + NoLegend()
VlnPlot(fivemonth_soupx, features = c("Col22a1"), idents = c("Type IIb Myonuclei", "Type IIx Myonuclei", "Neuromuscular Junction", "Myotendinous Junction", "Satellite Cells"), pt.size = 0) + NoLegend()
VlnPlot(fivemonth_soupx, features = c("Pax7"), idents = c("Type IIb Myonuclei", "Type IIx Myonuclei", "Neuromuscular Junction", "Myotendinous Junction", "Satellite Cells"), pt.size = 0) + NoLegend()
