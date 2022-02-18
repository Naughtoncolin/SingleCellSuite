pa4.filtered<- RunUMAP(pa4.filtered, dims = 1:50)
pa4.filtered <- FindNeighbors(pa4.filtered, dims = 1:50)
pa4.filtered <- FindClusters(pa4.filtered, resolution = .1)

library(DoubletFinder)
sweep.res.list_pa4 <- paramSweep_v3(pa4.filtered, PCs = 1:50, sct = TRUE)
sweep.stats_pa4 <- summarizeSweep(sweep.res.list_pa4, GT = FALSE)
bcmvn_pa4 <- find.pK(sweep.stats_pa4)
annotations <- pa4.filtered$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(pa4.filtered@meta.data))
pa4.doub <- doubletFinder_v3(pa4.filtered, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
pa4.doub$DF.classifications_0.25_0.09_250
