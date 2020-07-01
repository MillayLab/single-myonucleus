sol_soupx_muscle <- subset(sol_soupx, idents = c("Type I Myonuclei", "Type IIx Myonuclei", "Type IIa Myonuclei", "Myotendinous Junction", "Satellite Cells"))
fivemonth_soupx_muscle <- subset(fivemonth_soupx, idents = c("Neuromuscular Junction", "Type IIx Myonuclei", "Type IIx Myonuclei #2", "Type IIb Myonuclei #2", "Type IIb Myonuclei", "Myotendinous Junction", "Satellite Cells"))

fivemonth_subset_counts<- GetAssayData(object = fivemonth_soupx_muscle, slot = "counts")
soleus_subset_counts <- GetAssayData(object = sol_soupx_muscle, slot = "counts")


fivemonth_muscle_integration <- CreateSeuratObject(counts = fivemonth_subset_counts, project = "5 Month", min.cells = 3, min.features = 200)
fivemonth_muscle_integration <- subset(fivemonth_muscle_integration, subset = nFeature_RNA > 200 & nFeature_RNA < 3200)
fivemonth_muscle_integration$sol_muscle_integration <- "TA"
fivemonth_muscle_integration <- NormalizeData(fivemonth_muscle_integration, verbose = FALSE)
fivemonth_muscle_integration <- FindVariableFeatures(fivemonth_muscle_integration, selection.method = "vst", nfeatures = 2000)
sol_muscle_integration <- CreateSeuratObject(counts = soleus_subset_counts, project = "24 Month", min.cells = 3, min.features = 200)
sol_muscle_integration<- subset(sol_muscle_integration, subset = nFeature_RNA > 200 & nFeature_RNA < 3200)
sol_muscle_integration$sol_muscle_integration <- "Soleus"
sol_muscle_integration <- NormalizeData(sol_muscle_integration, verbose = FALSE)
sol_muscle_integration <- FindVariableFeatures(sol_muscle_integration, selection.method = "vst", nfeatures = 2000)
#Performing integration
aged.anchors.muscle <- FindIntegrationAnchors(object.list = list(fivemonth_muscle_integration, sol_muscle_integration), dims = 1:20)
sol_ta_muscle_integration <- IntegrateData(anchorset = aged.anchors.muscle, dims = 1:20)
#Cluster, integrated assay should get rid of batch
DefaultAssay(sol_ta_muscle_integration) <- "integrated"
sol_ta_muscle_integration <- ScaleData(sol_ta_muscle_integration, verbose = FALSE)
sol_ta_muscle_integration <- RunPCA(sol_ta_muscle_integration, npcs = 30, verbose = FALSE)
sol_ta_muscle_integration <- RunUMAP(sol_ta_muscle_integration, reduction = "pca", dims = 1:30)
sol_ta_muscle_integration <- FindNeighbors(sol_ta_muscle_integration, reduction = "pca", dims = 1:30)
sol_ta_muscle_integration <- FindClusters(sol_ta_muscle_integration, resolution = 0.5)
sol_ta_muscle_integration@reductions[["umap"]] <- sol_ta_reductions
DimPlot(sol_ta_muscle_integration, reduction = "umap", group.by = "sol_muscle_integration")
DimPlot(sol_ta_muscle_integration, reduction = "umap", label = TRUE)
DefaultAssay(sol_ta_muscle_integration) <- "RNA"
#Make Heatmap
sol_ta_muscle_integration <- subset(sol_ta_muscle_integration, idents = c("8", "9"), invert = TRUE)
sol_ta_muscle_integration <- RenameIdents(sol_ta_muscle_integration, "0" = "Type IIb Myonuclei", "1" = "Type IIx Myonuclei", "2" = "Type IIb Myonuclei", "3" = "Type I Myonuclei", "4" = "Musculotendinous Junction", "5" = "Type IIa Myonuclei",  "6" = "Satellite Cells", "7" = "Neuromuscular Junction")
sol_ta_muscle_integration_markers <- FindAllMarkers(sol_ta_muscle_integration, only.pos = TRUE)
top5 <- sol_ta_muscle_integration_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
sol_ta_muscle_integration <- ScaleData(sol_ta_muscle_integration)
DoHeatmap(sol_ta_muscle_integration, features = top5$gene) + NoLegend()
gene <- c("Myh4", "Mybpc2", "Actn3", "Sox6", "Pde4d", "Myh1", "Vegfa", "Sorbs1", "Actn2", "Tead1", "Myh7", "Tnnc1", "Tpm3", "Tnnt1", "Myl2", "Myh2", "Csrp3", "Ankrd2", "Fam129a", "Myoz2", "Lama2", "Col22a1", "Ankrd1", "Slc24a2", "Adamts20", "Pax7", "Chodl", "Notch3", "Fgfr4", "Vcam1", "Chrne", "Vav3", "Musk", "Ufsp1", "Colq")
