library(Seurat)
library(DropletQC)
library(patchwork)
setwd("~/Dropbox (GaTech)/Gibson/Working/Perianal Fistula/")
#setwd("C:\\Users\\Naugh\\Dropbox (GaTech)\\Gibson\\Data\\cellranger_count_outputs")
dataDir<- "~/Dropbox (GaTech)/Gibson/Data/cellranger_count_outputs"
dataDirList <- list.dirs(dataDir, recursive = FALSE)
paList = list()
for (i in 1:length(list.dirs(recursive = FALSE))){
#for (i in 1:1){
  paSinglet<-list.dirs(recursive = FALSE)[i]
  Objname <- sub("./", '', paSinglet)
  paList[[Objname]] <- foo(paSinglet)
  paList[[Objname]]<- bar(paSinglet, paList[[Objname]])
}


foo <- function(x){
  name <- sub("./", '', x)
  name.data<- Read10X_h5(paste0(x,"/outs/filtered_feature_bc_matrix.h5"))
  name<- CreateSeuratObject(name.data, project = name)
  name$percent.mt <- PercentageFeatureSet(name, pattern ="^MT-")
  #nucFrac <- nuclear_fraction_tags(
   # outs = paste0(x,"/outs"),
    #tiles = 1, cores = 11, verbose = FALSE)
  #name$nuclear_fraction <- nucFrac$nuclear_fraction
  return(name)
}


bar <- function(path, pa){
  print('here')
  if(!("nuclear_fraction" %in% colnames(pa[[]]))){
    print("Adding nuclear fractions")
    nucFrac <- nuclear_fraction_tags(
      outs = paste0(path,"/outs"),
      tiles = 1, cores = 11, verbose = FALSE)
    pa$nuclear_fraction <- nucFrac$nuclear_fraction
  }
  nf.df<- data.frame(pa$nuclear_fraction,pa$nCount_RNA)
  print("nucFrac calculated")
  #ed <- identify_empty_drops(nf_umi=nf.df, include_plot = TRUE)
  ed <- identify_empty_drops(nf_umi=nf.df, include_plot = TRUE, umi_rescue = 0)
  print("Empty Drops calculated")
  ed$cellType <- "unknown"
  dc <- identify_damaged_cells(ed, verbose = FALSE, output_plots = TRUE, nf_sep = .1)
  print("Damaged Cells calculated")
  print(wrap_plots(dc[[2]])) #For plotting damaged cell
  pa$DropletQC <- dc$df$cell_status
  return(pa)
  
}


paList <- readRDS("/home/colin/Dropbox (GaTech)/Gibson/Working/Perianal Fistula/paSamples_210118.RDS")
bar2 <- function(path, pa){
  #metadata<-colnames(pa[[]])
  #print(metadata)
  ifelse(!("nuclear_fraction" %in% colnames(pa[[]])), print("Adding nuclear fractions"),print("It's there"))
}
#bar <- function(pa){print(colnames(pa[[]]))}

for (i in 1:length(paList)){
  #bar2(dataDirList[[i]], paList[[i]])
  paList[[i]]<- bar(dataDirList[[i]], paList[[i]])
}
#colnames(paList[[1]][[]])

ep5<-ggplot(paList$PA_5[[]], aes(x=nuclear_fraction))+
  geom_density()+
  ggtitle("PA_5")
