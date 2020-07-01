p10.markers <- FindAllMarkers(p10_soupx)
Idents(p10_soupx) <- p10_soupx@meta.data[["seurat_clusters"]]
top10 <- p10.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(p10_soupx, features = top10$gene) + NoLegend()

p21.markers <- FindAllMarkers(p21_soupx)
Idents(p21_soupx) <- p21_soupx@meta.data[["seurat_clusters"]]
top10 <- p21.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(p21_soupx, features = top10$gene) + NoLegend()

fivemonth.markers <- FindAllMarkers(fivemonth_soupx)
Idents(fivemonth_soupx) <- fivemonth_soupx@meta.data[["seurat_clusters"]]
top10 <- fivemonth.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(fivemonth_soupx, features = top10$gene) + NoLegend()

sol.markers <- FindAllMarkers(sol_soupx)
Idents(sol_soupx) <- sol_soupx@meta.data[["seurat_clusters"]]
top10 <- sol.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(sol_soupx, features = top10$gene) + NoLegend()

twentyfour.markers <- FindAllMarkers(twentyfour_soupx)
Idents(twentyfour_soupx) <- twentyfour_soupx@meta.data[["seurat_clusters"]]
top10 <- twentyfour.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(twentyfour_soupx, features = top10$gene) + NoLegend()

thirty.markers <- FindAllMarkers(thirty_soupx)
Idents(thirty_soupx) <- thirty_soupx@meta.data[["seurat_clusters"]]
top10 <- thirty.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(thirty_soupx, features = top10$gene) + NoLegend()
