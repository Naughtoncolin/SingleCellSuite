library(Seurat)
set.seed(1991)
memory.limit(25000)# If you have memery problems you might need to run this

#Multithreading setup; NOTE: Doesn't work in RStudio!
library(future)
#options(future.globals.maxSize= 7000*1024^2)
#options(future.globals.maxSize= 1500*1024^2)
#plan(strategy = "multicore", workers = 8)
options(future.globals.maxSize= 5240*1024^2)
plan(strategy = "multicore")

# Set working directory, raw file directory, and initialize list to hold Seurat objects
setwd("C:\\Users\\Naugh\\Dropbox (GaTech)\\School\\BIOL-6150\\Project")
rawDirs <- Sys.glob("raw/PBMC-BAL.*")
pbmcList = list()

#---------------------Create Seurat Objects and add HTO assays-------------------------------------------

for (i in 6:length(rawDirs)){ #This loops through the raw data files & creates multiple Seurat objects
	pbmc.data <- Read10X_h5(rawDirs[i]) #Read raw data file
	try(pbmcList[[paste0("pbmc", i)]] <- CreateSeuratObject(pbmc.data, project = "BIOL-6150"), silent = TRUE) #Create Seurat Object for donor with only GEX data
	try(pbmcList[[paste0("pbmc", i)]] <- CreateSeuratObject(pbmc.data$"Gene Expression", project = "BIOL-6150"), silent = TRUE) #Create Seurat Object
	try(pbmcList[[paste0("pbmc", i)]][["HTO"]] <- mar(pbmc.data$`Antibody Capture`[8:9,]), silent = TRUE) # Add HTO data as a new assay to newly created Seurat object, determines PBMCvsBAL
	try(pbmcList[[paste0("pbmc", i)]][["TcellHTO"]] <- CreateAssayObject(pbmc.data$`Antibody Capture`[1:2,]), silent = TRUE) # Add HTO data as a new assay to newly created Seurat object, determines CD4+vsCD8+
	pbmcList[[paste0("pbmc", i)]][["Sample"]] <- paste0("Donor.",i)
}
rm(rawDirs, pbmc.data, i)
gc()

#----------------Main block for processing Seurat objects pre-integration---------------------------------

pbmcList <- lapply(pbmcList, function(x) {
	x$Disease.Status  <- "Diseased" #For some reason adding metadata only works if implemented at the beginning of the lapply block.
	
	# Filtering by number of mitochondrial reads
	x$percent.mt <-  PercentageFeatureSet(x, pattern = "^MT-") #Get mitochondrial counts for metadata
	x <- subset(x, percent.mt < 10)

	#Block dealing with donors that had HTO(antibodies with attached oligos) data. 
	try(x <- NormalizeData(x, assay = "HTO", normalization.method = "CLR"), silent=TRUE)
	try(x <- HTODemux(x, assay = "HTO", positive.quantile = 0.99), silent=TRUE)
	try(x <- NormalizeData(x, assay = "TcellHTO", normalization.method = "CLR"), silent=TRUE)
	try(x <- HTODemux(x, assay = "TcellHTO", positive.quantile = 0.7), silent=TRUE)
	try(x <- subset(x, HTO_classification == "Hashtag2TotalSeqC" & HTO_classification.global== "Singlet"), silent=TRUE)#Subset cells based on high confidence of PBMC surface markers. The paper's authors indicated "Hashtag1TotalSeqC" were PBMCs

	#Normalization Commands
	#x <- NormalizeData(x)
	x <- SCTransform(x, method = "glmGamPoi") #Requires installation of glmGamPoi, runs faster.
	#x <- SCTransform(x)
	})

pbmcList[['pbmc7']][['Disease.Status']]<- 'Healthy'
gc()

#-----------------------------------Integration Steps---------------------------------------------------------------

# select integration features and prep step
features <- SelectIntegrationFeatures(pbmcList)
pbmc.anchors <- PrepSCTIntegration(pbmcList, anchor.features = features) #https://rdrr.io/cran/Seurat/man/PrepSCTIntegration.html
#rm(pbmcList) # Is this necessary? To prevent Challenger from crashing
gc()

# downstream integration steps
#plan(strategy = "multicore", workers = 4)
pbmc.anchors <- FindIntegrationAnchors(pbmc.anchors, anchor.features = features )
gc()
pbmc.combined <- IntegrateData(pbmc.anchors) # this command creates an 'integrated' data assay

rm(pbmc.anchors, features)
gc()


#-------------------------------Post-Integration Steps----------------------------------------------------------
# specify that we will perform downstream analysis on the corrected "integrated" data; note that the
# original unmodified data still resides in the 'RNA' assay
pbmc.combined <- ScaleData(pbmc.combined) #SCTransform is supposed to ScaleData, but following commands won't work without it.

DefaultAssay(pbmc.combined) <- "integrated"

pbmc.combined <- ScaleData(pbmc.combined) #SCTransform is supposed to ScaleData, but following commands won't work without it.
pbmc.combined <- RunPCA(pbmc.combined, npcs = 30)
ElbowPlot(pbmc.combined) # Used to determine dimensionality.
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:10) #Should I use spca instead of pca?
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:10)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5) 

gc()
#SCTransform is only good for clustering!!!!

# Clustering Plots
DimPlot(pbmc.combined, reduction = "umap", group.by = "Disease.Status")
DimPlot(pbmc.combined, reduction = "umap", group.by = "seurat_clusters")
Idents(pbmc.combined) <- "Disease.Status"
DefaultAssay(pbmc.combined) <- "SCT"
#markersSCT <- FindMarkers(pbmc.combined, ident.1 = "Diseased", ident.2 = "Healthy")
#markersSCT

#write.csv(markersSCT, "VariableGenes.211019.csv")



# Disease vs Healthy Heat Maps
library(dplyr)
DefaultAssay(pbmc.combined) <- "SCT"
Idents(pbmc.combined) <- 'Disease.Status'
pbmc.diseased.pos.markers <- FindAllMarkers(pbmc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
pbmc.diseased.pos.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top10.pos.diseased
#DoHeatmap(pbmc.combined, features = top10.pos.diseased$gene)
DoHeatmap(subset(pbmc.combined, downsample=3000), features = top10.pos.diseased$gene)

pbmc.diseased.markers<- FindAllMarkers(pbmc.combined, min.pct = 0.25, logfc.threshold = 1)
pbmc.diseased.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top10.diseased
#DoHeatmap(pbmc.combined, features = top10.diseased$gene)
DoHeatmap(subset(pbmc.combined, downsample=3000), features = top10.diseased$gene)

pbmc.diseased.markers2<- FindAllMarkers(pbmc.combined, min.pct = 0.25, logfc.threshold = 1)
pbmc.diseased.markers2 %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top10.2.diseased
#DoHeatmap(pbmc.combined, features = top10.diseased$gene)
DoHeatmap(subset(pbmc.combined, downsample=3000), features = top10.2.diseased$gene)





# Cluster Heat Maps
Idents(pbmc.combined) <- 'seurat_clusters'
pbmc.cluster.markers <- FindAllMarkers(pbmc.combined, min.pct = 0.25, logfc.threshold = 1)
pbmc.cluster.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top10.cluster
#DoHeatmap(pbmc.combined, features = top10.cluster$gene)
DoHeatmap(subset(pbmc.combined, downsample=3000), features = top10.cluster$gene)

Idents(pbmc.combined) <- 'seurat_clusters'
pbmc.cluster.pos.markers <- FindAllMarkers(pbmc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
pbmc.cluster.pos.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top10.pos.cluster
#DoHeatmap(pbmc.combined, features = top10.pos.cluster$gene)
DoHeatmap(subset(pbmc.combined, downsample=3000), features = top10.pos.cluster$gene)

Idents(pbmc.combined) <- 'seurat_clusters'
pbmc.cluster.pos.markers2<- FindAllMarkers(pbmc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
pbmc.cluster.pos.markers2 %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top10.pos.cluster2
#DoHeatmap(pbmc.combined, features = top10.pos.cluster$gene)
DoHeatmap(subset(pbmc.combined, downsample=3000), features = top10.pos.cluster2$gene)

write.csv(pbmc.diseased.markers, "diseased.markers.211027.csv")
write.csv(pbmc.cluster.markers, "cluster.markers.211027.csv")

Idents(pbmc.combined) <- 'Disease.Status'
Diseased.ConMarkers <- FindConservedMarkers(pbmc.combined)

Diseased.ConMarkers <- FindConservedMarkers(pbmc.combined, grouping.var='Disease.Status')



par(las=2)
barplot(pbmc.combined@reductions$pca@stdev, ylab="Standard Deviation Explained", names.arg = PCnames)


#Altered Elbow Plot; plots variance explained by all PCs up to that point
foo <- pbmc.combined@reductions$pca@stdev
sum =0
for(i in foo){
    sum = sum + i
}

counter = 0
for(i in 1:length(foo)){
    counter = counter + foo[i]
    print(counter)
    var.sum[i] <- (counter/sum)*100
}
par(las=2)
barplot(var.sum, names.arg = 1:30)