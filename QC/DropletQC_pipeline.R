library(DropletQC)
library(patchwork)
library(Seurat)
#outDir = "~/Dropbox (GaTech)/Gibson/Data/cellranger_count_outputs/p21090-s002_P4/outs"
#barcodeDir = "~/Dropbox (GaTech)/Gibson/Data/cellranger_count_outputs/p21090-s002_P4/outs/filtered_feature_bc_matrix.h5"
barcodeDir = "~/Dropbox (GaTech)/Gibson/Data/cellranger_count_outputs/PA_7/outs/filtered_feature_bc_matrix.h5"
outDir = "~/Dropbox (GaTech)/Gibson/Data/cellranger_count_outputs/PA_7/outs"
nf7 <- nuclear_fraction_tags(
  outs = outDir,
  tiles = 1, cores = 11, verbose = FALSE)

# Create Seurat Object
#pa7 <- Read10X_h5("C:\\Users\\Naugh\\Dropbox (GaTech)\\Gibson\\Data\\cellranger_count_outputs\\p21090-s002_P4\\outs\\filtered_feature_bc_matrix.h5")
pa7.data <- Read10X_h5(barcodeDir)
pa7 <- CreateSeuratObject(pa7.data)
pa7$percent.mt <- PercentageFeatureSet(pa7, pattern ="^MT-")
#pa7.umi <- data.frame(pa7$nCount_RNA)


# Identify empty droplets
# Uses umi_rescue default, 5004 (local minima), and 89343 (covers everything)
test7 <- data.frame(nf7$nuclear_fraction,pa7$nCount_RNA)
ed7.1 <- identify_empty_drops(nf_umi=test7, include_plot = TRUE)
ed7.2 <- identify_empty_drops(nf_umi=test7, include_plot = TRUE, umi_rescue = 5004) # umi_rescue determined by local minima of nCount_RNA density distribution
ed7.3 <- identify_empty_drops(nf_umi=test7, include_plot = TRUE, umi_rescue = 89343) # umi_rescue determined by "max(test4$pa6.nCount_RNA)+1"

#table(ed4.3$cell_status)

# Identify damaged cells
# Cell Types set to unknown because no prior knowledge about cluster identities, but knowing would be prefered.
# Does this mean this pipeline should also be run after cell labeling? Is using "unknown" even appropriate?
ed7.1$cellType <- "unknown"
ed7.2$cellType <- "unknown"
ed7.3$cellType <- "unknown"
dc7.1 <- identify_damaged_cells(ed7.1, verbose = FALSE, output_plots = TRUE)
dc7.2 <- identify_damaged_cells(ed7.2, verbose = FALSE, output_plots = TRUE)
dc7.3 <- identify_damaged_cells(ed7.3, verbose = FALSE, output_plots = TRUE)
dc7.3.2 <- identify_damaged_cells(ed7.3, verbose = FALSE, output_plots = TRUE, nf_sep=14)
wrap_plots(dc7.1[[2]]) #For plotting damaged cell

# Assign cell status labels to cells in original Seurat object. Make violin plot of QC metrics grouped by cell status.

pa7$DropletQC1 <- dc7.1$df$cell_status
pa7$DropletQC2 <- dc7.2$df$cell_status
pa7$DropletQC3 <- dc7.3.2$df$cell_status
#pa7.filtered <- subset(pa7, subset = pa7$DropletQC = 'cell')
VlnPlot(pa7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "DropletQC1")
VlnPlot(pa7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "DropletQC2")
VlnPlot(pa7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "DropletQC3")

# Following steps are for creating violin plots that include the original distribution
pa7.1 <- Read10X_h5("~/Dropbox (GaTech)/Gibson/Data/cellranger_count_outputs/p21090-s002_P4/outs/filtered_feature_bc_matrix.h5")
pa7.1 <- CreateSeuratObject(pa7.1)
pa7.1$percent.mt <- PercentageFeatureSet(pa7.1, pattern ="^MT-")
pa7.1 <- pa7
pa7.1$DropletQC1 <- "original"
pa7.1$DropletQC2 <- "original"
pa7.1$DropletQC3 <- "original"
pa7.merge <- merge(pa7, pa7.1)
# The following two commands are to order the DropletQC metadata columns so they show up with "original" first
pa7.merge$DropletQC1 <- factor(pa7.merge$DropletQC1 , levels = c("original", "cell", "empty_droplet"))
pa7.merge$DropletQC2 <- factor(pa7.merge$DropletQC2 , levels = c("original", "cell", "damaged_cell", "empty_droplet"))
pa7.merge$DropletQC3 <- factor(pa7.merge$DropletQC3 , levels = c("original", "cell", "damaged_cell", "empty_droplet"))
o1<- VlnPlot(pa7.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "DropletQC1")
o2<- VlnPlot(pa6.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "DropletQC2")
o3<- VlnPlot(pa6.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "DropletQC3")

length(colnames(pa6)) - length(colnames(pa6.filtered))
length(colnames(pa6.filtered))
ggpubr::ggarrange(
  VlnPlot(pa6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0),
  VlnPlot(pa6.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0),
  nrow = 2)

#Check for apoptotic markers
CD45Pos <- WhichCells(pf, expression = PTPRC > 0)
pf$CommonImmuneCD45 <- ifelse(colnames(pf) %in% CD45Pos, "POS", "NEG")


pa6.filtered <- subset(pa6.filtered, subset = percent.mt < 15)
pa6.filtered <-SCTransform(pa6.filtered)
pa6.filtered <- RunPCA(pa6.filtered)
stDeviations <- pa6.filtered@reductions$pca@stdev
sum =0
for(i in stDeviations){
sum = sum + i
}
var.sum <-vector()
counter = 0
for(i in 1:length(stDeviations)){
counter = counter + stDeviations[i]
print(counter)
var.sum[i] <- (counter/sum)*100
}
par(las=2)
barplot(var.sum, names.arg = 1:50)
pa6.filtered<- RunUMAP(pa6.filtered, dims = 1:50)
pa6.filtered <- FindNeighbors(pa6.filtered, dims = 1:50)
pa6.filtered <- FindClusters(pa6.filtered, resolution = .1)

