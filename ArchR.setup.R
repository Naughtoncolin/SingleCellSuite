#Setup
library(ArchR)
set.seed(1991)

#Set Threads
#addArchRThreads(threads = 1) #What is optimum?

#Set input files
#inputFiles= c("atac_fragments.tsv.gz")
inputFiles <- getInputFiles(path="./raw")

#Add Genome
#addArchRGenome("hg38") #Already added

#Create Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = "PBMC.ATAC.1",
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "IBD.PBMC",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

#Find out what matrices are available with the project.
#getAvailableMatrices(proj)

#Filter Doublets
proj <- filterDoublets(ArchRProj = proj)


proj <- addIterativeLSI(ArchRProj = proj, 
	useMatrix = "TileMatrix", 
	name = "IterativeLSI",
	dimsToUse=1:10)

#Add Clusters using Seurat Method.
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

#Project scATAC -eq
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")


plotPDF(p1,p2, name = "PBMC.Clusters.UMAP.1.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)




proj <- addImputeWeights(proj)


#STOPPING POINT
markerGenes  <- c(
    "CD34",  #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "MME", #B-Cell Trajectory
    "CD14", "MPO", #Monocytes
    "CD3D", "CD8A"#TCells
  )

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

#Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(plotList = p, 
    name = "PBMC.MarkerGeneImpute.UMAP.1.pdf",
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

#Browser Track
p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000
)

grid::grid.newpage()
#grid::grid.draw(p$CD14)

plotPDF(plotList = p, 
    name = "PBMC.MarkerGene.Tracks.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)
