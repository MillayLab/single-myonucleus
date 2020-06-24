set.seed(42)
p21 <- CreateSeuratObject(counts = p21_counts, project = "p21", min.cells = 3, min.features = 200)
p10 <- CreateSeuratObject(counts = p10_counts, project = "p10", min.cells = 3, min.features = 200)
fivemonth <- CreateSeuratObject(counts = fivemonthcounts, project = "fivemonth", min.cells = 3, min.features = 200)
twentyfourmonth <- CreateSeuratObject(counts = twentyfourcounts, project = "twentyfour", min.cells = 3, min.features = 200)
thirtymonth <- CreateSeuratObject(counts = thirtycounts, project = "thirty", min.cells = 3, min.features = 200)

#Integrated All Datasets
dataset.list <- c(p10, p21, fivemonth, twentyfourmonth, thirtymonth)
for (i in 1:length(dataset.list)) {
  dataset.list[[i]] <- NormalizeData(dataset.list[[i]], verbose = FALSE)
  dataset.list[[i]] <- FindVariableFeatures(dataset.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = FALSE)
}
dataset.list[[1]]@meta.data[["orig.ident"]] <- "p10"
dataset.list[[2]]@meta.data[["orig.ident"]] <- "p21"
dataset.list[[3]]@meta.data[["orig.ident"]] <- "fivemonth"
dataset.list[[4]]@meta.data[["orig.ident"]] <- "twentyfour"
dataset.list[[5]]@meta.data[["orig.ident"]] <- "thirty"

dataset.anchors <- FindIntegrationAnchors(object.list = dataset.list, dims = 1:30)
dataset.integrated <- IntegrateData(anchorset = dataset.anchors, dims = 1:30)
DefaultAssay(dataset.integrated) <- "integrated"




# Run the standard workflow for visualization and clustering
dataset.integrated <- ScaleData(dataset.integrated, verbose = FALSE)
dataset.integrated <- RunPCA(dataset.integrated, npcs = 30, verbose = FALSE)
dataset.integrated <- FindNeighbors(dataset.integrated, dims = 1:30)
dataset.integrated <- FindClusters(dataset.integrated, resolution = 0.5)
dataset.integrated <- RunUMAP(dataset.integrated, reduction = "pca", dims = 1:30)
DimPlot(dataset.integrated, reduction = "umap", label = TRUE, 
        repel = TRUE) + NoLegend()
DimPlot(dataset.integrated, reduction = "umap", label = TRUE, pt.size = 2.5) + NoAxes() + NoLegend()