sol_soupx <- CreateSeuratObject(counts = sol_counts, project = "sol", min.cells = 3, min.features = 200)
#Dimensionality Plots
VlnPlot(sol_soupx, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
featurescatterplot <- FeatureScatter(sol_soupx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot(featurescatterplot)
#Subset the data, used < 3200
sol_soupx <- subset(sol_soupx, subset = nFeature_RNA > 200 & nFeature_RNA < 3200)
#Normalize Data
sol_soupx <- NormalizeData(sol_soupx)
sol_soupx <- FindVariableFeatures(sol_soupx, selection.method = "vst", nfeatures = 2000)
#Top 10 variable features labelled on a variable feature plot
top10 <- head(VariableFeatures(sol_soupx), 10)
plot1variable <- VariableFeaturePlot(sol_soupx)
plot2variable <- LabelPoints(plot = plot1variable, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1variable, plot2variable))
#Scaling data
all.genes5mo <- rownames(sol_soupx)
sol_soupx <- ScaleData(sol_soupx, features = all.genes5mo)
#PCA
sol_soupx <- RunPCA(sol_soupx, features = VariableFeatures(object = sol_soupx))
#VizDimLoadings
VizDimLoadings(sol_soupx, dims = 1:2, reduction = "pca")
#DimPlot
DimPlot(sol_soupx, reduction = "pca")
#DimHeatMap
DimHeatmap(sol_soupx, dims = 1:15, cells = 500, balanced = TRUE)
#UMAP, dimensions to test, 12 looked best
sol_soupx <- FindNeighbors(sol_soupx, dims = 1:20)
sol_soupx <- FindClusters(sol_soupx, resolution = 0.5)
sol_soupx <- RunUMAP(sol_soupx, dims = 1:20)
sol_soupx@reductions[["umap"]] <- sol_reductions
DimPlot(sol_soupx, reduction = "umap", label = TRUE)
sol_soupx <- RenameIdents(sol_soupx, "0" = "Type IIa Myonuclei", "1" = "FAPs", "2" = "Type IIx Myonuclei", "3" = "Type I Myonuclei", "4" = "Endothelial Cells", "5" = "Smooth Muscle", "6" = "Tenocytes", "7" = "Myotendinous Junction", "8" = "Immune Cells", "9" = "Smooth Muscle #2", "10" = "Adipocytes", "11" = "Satellite Cells", "12" = "FAPs #2", "13" = "Unnamed", "14" = "Schwann Cells", "15" = "Neuronal")
DimPlot(sol_soupx, reduction = "umap", label = TRUE)
