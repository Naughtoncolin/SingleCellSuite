library(future)

options(future.globals.maxSize= 8000*1024^2)
plan(strategy = "multicore", workers = 10)
library(Seurat)

rawDirs <- Sys.glob("MultiomeData/M*/filtered_feature_bc_matrix.h5")
#length(rawDirs)
#---------------------Create Seurat Objects-------------------------------------------
for (i in 1:length(rawDirs)){ #This loops through the raw data files & creates multiple Seurat objects
	pbmc.data <- Read10X_h5(rawDirs[i]) #Read raw data file
	assign( paste0("pbmc.RNA.M", i) , CreateSeuratObject(counts = pbmc.data$"Gene Expression", project = paste0("IBD.PBMC.RNA.M", i))) #Create Seurat Object
	placeholder <- get(paste0("pbmc.RNA.M", i)) #Refer to new seurat object with another name
	placeholder[['percent.mt']] <-  PercentageFeatureSet( placeholder, pattern = "^MT-") #Get mitochondrial counts for metadata

	assign( paste0("pbmc.RNA.M", i), placeholder)
	rm(placeholder, pbmc.data)
}
gc()


#----------------------View QC Metrics------------------------------------------------
for (i in 1:length(rawDirs)){
	assign( paste0("p", i), VlnPlot(get(paste0("pbmc.RNA.M", i)), features =  "percent.mt", y.max=30))
}

p1+p2+p3+p4+p5+p6+p7+p8+p9+p10

VlnPlot(pbmc.RNA.M1, features =  "percent.mt")
ggAlignPlot(p1, p2)

FindMarkers on outlier cluster


#--------------------------Normalize Data----------------------------------------------------
for (i in 1:length(rawDirs)){
	placeholder <- get(paste0("pbmc.RNA.M", i)) #Refer to new seurat object with another name
	#IBD.PBMC.RNA.M1 <- NormalizeData(pbmc.RNA.M1)

	#placeholder <- SCTransform(placeholder, verbose = FALSE)
	placeholder <- SCTransform(placeholder)
	assign( pte0("pbmc.RNA.M", i), placeholder)
	rm(placeholder, pbmc.data)
}

testList = list()
for (i in 1:length(rawDirs)){
	#placeholder <- get(paste0("pbmc.RNA.M", i)) #Refer to new seurat object with another name
	testList[paste0("pbmc.RNA.M", i)] <- get(paste0("pbmc.RNA.M", i))
	

}

features <- SelectIntegrationFeatures(object.list = testList) ##Can different reduction methods be used?
features <- SelectIntegrationFeatures(object.list = pbmcList) ##Can different reduction methods be used?
pbmc.anchors <- PrepSCTIntegration(object.list = pbmcList, anchor.features = features)
pbmc.anchors <- PrepSCTIntegration(object.list = pbmcList)

pbmc.anchors <- FindIntegrationAnchors(object.list = testList, anchor.features = features)
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.anchors)
pbmc.anchors <- FindIntegrationAnchors(testList)

# this command creates an 'integrated' data assay
pbmc.combined <- IntegrateData(anchorset = pbmc.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(pbmc.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE) #SCTransform is supposed to ScaleData, but following commands won't work without it.
pbmc.combined <- RunPCA(pbmc.combined, npcs = 30, verbose = TRUE)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:30) #Should I use spca instead of pca?
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:30)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5)

p1 <- DimPlot(pbmc.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(pbmc.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1+p2

DimPlot(pbmc.combined, reduction = "umap", split.by = "orig.ident")