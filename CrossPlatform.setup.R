library(ArchR)

# Load unmodified Seurat object
pbmc <- readRDS("pbmc.Seurat.rds")
pbmc

library(ArchR)
RNA.ATAC.Integration <- loadArchRProject("IBD.PBMC")

DefaultAssay(pbmc) <- "RNA"

#Un/constrained Integration l1
pbmcInt <- addGeneIntegrationMatrix(
    ArchRProj = RNA.ATAC.Integration, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix_l1",
    reducedDims = "IterativeLSI",
    seRNA = pbmc,
    addToArrow = TRUE, #Set to FALSE during Unconstrained integration. 
    groupRNA = "predicted.celltype.l1",
    nameCell = "predictedCell_l1_Un",
    nameGroup = "predictedGroup_l1_Un",
    nameScore = "predictedScore_l1_Un",
    force =TRUE
)

#Un/constrained Integration l2 (Why does matrixName already show up as taken?)
pbmcInt <- addGeneIntegrationMatrix(
    ArchRProj = RNA.ATAC.Integration, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix_l2",
    reducedDims = "IterativeLSI",
    seRNA = pbmc,
    addToArrow = TRUE, #Set to FALSE during Unconstrained integration. 
    groupRNA = "predicted.celltype.l2",
    nameCell = "predictedCell_l2_Un",
    nameGroup = "predictedGroup_l2_Un",
    nameScore = "predictedScore_l2_Un",
    force =TRUE
)


#Constrained Integration
cM <- as.matrix(confusionMatrix(pbmcInt$Clusters, pbmcInt$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments


unique(unique(pbmcInt$predictedGroup_Un))

#Label scATAC-seq clusters with scRNA-seq info.
cM <- as.matrix(confusionMatrix(pbmcInt$Clusters, pbmcInt$predictedGroup_l1_Un))
labelOld <- rownames(cM)
labelOld

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

pbmcInt$Clusters.l1 <- mapLabels(pbmcInt$Clusters, newLabels = labelNew, oldLabels = labelOld)

#Plot ATAC-sec UMAP with new cluster identities
p0 <- plotEmbedding(ArchRProj = RNA.ATAC.Integration, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1 <- plotEmbedding(pbmcInt, colorBy = "cellColData", name = "Clusters.l1")
p2 <- plotEmbedding(pbmcInt, colorBy = "cellColData", name = "Clusters.l2")
#p1

plotPDF(p0,p1,p2, name = "IBD.PBMC.ATAC.ClusterLabelComparison.pdf", ArchRProj = pbmcInt, addDOC = TRUE, width = 6, height = 6)


#Make pseudo-bulk replicates from clusters. NOTE: Multithreading didn't work with 12 threads.
pbmcInt <- addGroupCoverages(pbmcInt, groupBy = "Clusters.l1", threads = 1)

#Call peaks with MACS2
pathToMacs2 <- findMacs2()
pbmcInt <- addReproduciblePeakSet(
    ArchRProj = pbmcInt, 
    groupBy = "Clusters.l1", 
    pathToMacs2 = pathToMacs2
)

#Add Peak Matrix
pbmcInt <- addPeakMatrix(pbmcInt)

#Identify Marker Peaks in each cluster
markersPeaks <- getMarkerFeatures(
    ArchRProj = pbmcInt, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters.l1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#Get specific parts of the markersPeaks object
#markerList.Peaks <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
#markerList

#Get marker peaks for specific cell group.
#markerList$Erythroid

#Make & visualize marker peak heat map.
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
heatmapPeaks <- markerHeatmap( # markerHeatmap is depreciated.
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = pbmcInt, addDOC = TRUE)


#Add motif annotations
#NOTE: CIS-BP Database= Catalog of Inferred Sequence Binding Preferences
pbmcInt <- addMotifAnnotations(ArchRProj = pbmcInt, motifSet = "cisbp", name = "Motif")

#Find & visualize motifs enriched in marker peaks.
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = pbmcInt,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = pbmcInt, addDOC = TRUE)

#Add other enrichment annotations
pbmcInt <- addArchRAnnotations(ArchRProj = pbmcInt, collection = "EncodeTFBS")
enrichEncode <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = pbmcInt,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
#enrichEncode
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = pbmcInt, addDOC = TRUE)


#Identify Marker Genes in each cluster
markersGS <- getMarkerFeatures( #GS stands for Gene Score
    ArchRProj = pbmcInt, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters.l1",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & abs(Log2FC) >= 1.25")
#markerList$C6

markerGenes  <- c( #Optional
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
    "CD14", "CEBPB", "MPO", #Monocytes
    "IRF8", 
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells
  )

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes, #Include if you want specific genes labeled in the heatmap.
  transpose = TRUE
)

#ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = pbmcInt, addDOC = TRUE)


#ChromVAR motif deviation analysis



saveArchRProject(ArchRProj = pbmcInt, outputDirectory = "IBD.PBMC", load = FALSE)