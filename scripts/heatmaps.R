Idents(p10_soupx) <- p10_soupx@meta.data[["seurat_clusters"]]
top10 <- p10.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(p10_soupx, features = top10$gene) + NoLegend()

Idents(p21_soupx) <- p21_soupx@meta.data[["seurat_clusters"]]
top10 <- p21.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(p21_soupx, features = top10$gene) + NoLegend()

Idents(fivemonth_soupx) <- fivemonth_soupx@meta.data[["seurat_clusters"]]
top10 <- fivemonth.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(fivemonth_soupx, features = top10$gene) + NoLegend()

Idents(sol_soupx) <- sol_soupx@meta.data[["seurat_clusters"]]
top10 <- sol.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(sol_soupx, features = top10$gene) + NoLegend()

Idents(twentyfour_soupx) <- twentyfour_soupx@meta.data[["seurat_clusters"]]
top10 <- twentyfour.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(twentyfour_soupx, features = top10$gene) + NoLegend()

Idents(thirty_soupx) <- thirty_soupx@meta.data[["seurat_clusters"]]
top10 <- thirty.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(thirty_soupx, features = top10$gene) + NoLegend()
