# This script should take our PBMC data and map to query, doing
# everything but merge the data with the reference.
setwd("C:\\Users\\Naugh\\ShitTest")

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

# Load Reference Dataset
reference <- LoadH5Seurat("./raw/pbmc_multimodal.h5seurat")
#DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

# Load Query Dataset
pbmc_Read10X <- Read10X_h5("./raw/filtered_feature_bc_matrix.h5")
pbmc <- CreateSeuratObject(counts = pbmc_Read10X$"Gene Expression")

# QC for doublets & dying cells
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2

#Subset query data based on QC metrics observed above.
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize Query Dataset
pbmc <- SCTransform(pbmc, verbose = FALSE)
gc()

# Find anchors between reference and query datasets and map in low-dim space.
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
gc()

# Transfer cell type labels & protein data from reference->query
# Project query onto reference's UMAP
pbmc <- MapQuery(
  anchorset = anchors,
  query = pbmc,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

#p1 = DimPlot(pbmc, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
#p2 = DimPlot(pbmc, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
#p1 + p2

#Idents(pbmc) <- 'predicted.celltype.l2'