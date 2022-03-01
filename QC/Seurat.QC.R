set.seed(1991)
setwd("/home/colin/Dropbox (GaTech)/Gibson/Working")
#memory.limit(25000)# Run if having memory issues; Windows OS specific

#Multithreading setup; NOTE: Doesn't work in RStudio!
library(future)
#plan(strategy = "multicore", workers = 8)
plan(strategy = "multicore") # Use all cores available
options(future.globals.maxSize= 5240*1024^2) # Memory allocation per core (I think)

# Get list of paths to raw data (sorted)
library(gtools)
rawDirs <- mixedsort(Sys.glob("../Data/MultiomeData/M*/filtered_feature_bc_matrix.h5"))
detach("package:gtools", unload=TRUE)


library(Seurat)
#---------------------Create Seurat Objects and add HTO assays-------------------------------------------
pbmcList = list()
for (i in 1:length(rawDirs)){ #This loops through the raw data files & creates multiple Seurat objects
	pbmc.data <- Read10X_h5(rawDirs[i]) #Read raw data file
	seq.dir <- strsplit(rawDirs[i], '/')[[1]][4] # Extract specific name of directory sequencing files are in, e.g. M1 or M2
	pbmcList[[seq.dir]] <- CreateSeuratObject(pbmc.data$"Gene Expression", project = seq.dir) #Create Seurat Object
}
rm(rawDirs, pbmc.data, seq.dir, i) 


#----------------Main block for processing Seurat objects pre-integration---------------------------------
#Get mitochondrial counts for metadata
pbmcList <- lapply(pbmcList, function(x){x[['percent.mt']] <-  PercentageFeatureSet(x, pattern = "^MT-"); x})

# Plot QC metrics, each sample on different page.
pdf(file="QCplots_211029.pdf")
plot.list <- list()
for (i in 1:length(pbmcList)){
	plot.list[[i]] <- VlnPlot(pbmcList[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) #For viewing plots outside of PDF
	#plot(VlnPlot(pbmcList[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
}
dev.off()

p1<- VlnPlot(pbmcList[[1]], features = "percent.mt", pt.size = 0)

p2<- VlnPlot(pbmcList[[2]], features = "percent.mt", pt.size = 0)


pbmcList <- lapply(pbmcList, function(x){x <-  SCTransform(x)})
foo <-  merge(pbmcList[[1]], pbmcList[2:length(pbmcList)])
# Make Violin plots for merged seurat objects
VlnPlot(foo, features = "percent.mt", pt.size = 0, group.by = 'orig.ident') 
VlnPlot(foo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = 'orig.ident')

options(future.globals.maxSize= 8000*1024^2)
foo <- SCTransform(foo)
foo <- SCTransform(foo, conserve.memory=TRUE, verbose = TRUE)



foo <- RunPCA(foo, npcs = 30, verbose = TRUE)
foo <- RunUMAP(foo, reduction = "pca", dims = 1:30) #Should I use spca instead of pca?
foo <- FindNeighbors(foo, reduction = "pca", dims = 1:30)
foo <- FindClusters(foo, resolution = 0.5)

DimPlot(pbmc.combined, reduction = "umap", group.by = "Disease.Status")
DimPlot(foo, reduction = "umap", group.by = "seurat_clusters")