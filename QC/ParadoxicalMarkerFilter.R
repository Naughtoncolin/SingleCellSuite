#pa6.filtered <- subset(pa6.filtered, TFF3 >0 & PTPRC>0, invert = TRUE)
EpithelialMarkers <- c("EPCAM", "TFF3", "PHGR1")
othermarkers <- c("PTPRC", "CD3D")
TcellMarkers <- c("CD3D", "CD69", "CD8A") #First 2 from paper
markerLists <- list(EpithelialMarkers, othermarkers, TcellMarkers)
markerListCombos<- combn(markerLists, 2, simplify = FALSE)

#EpithelialMarkers <- c("EPCAM", "TFF3", "PHGR1")
EpithelialMarkers <- c("EPCAM", "TFF3", "PHGR1", "MUC2")
TcellMarkers <- c("CD3D", "CD69", "CD8A") #First 2 from paper
BcellMarkers <-  c('MS4A1', 'CD19', "CD79B")
NeuralMarkers <- c("ETV1", "BNC2") # Not good
#MesenchymalMarkers <- c("VIM", "S100A4", "COL4A2")
MesenchymalMarkers <- c("VIM", "HAND1", "HAND2", "PITX1", "ZEB2")
EndothelialMarkers <- c("PECAM1", "CDH5") 
MyeloidMarkers <- c("S100A4", "S100A6", "HLA-DRB1") #Seemingly shitty?
markerLists <- list(EpithelialMarkers, TcellMarkers, BcellMarkers, MesenchymalMarkers, EndothelialMarkers )
#markerLists <-list(EpithelialMarkers, "PTPRC", "VIM", EndothelialMarkers, MastMarker, MyeloidMarkers)
markerListCombos<- combn(markerLists, 2, simplify = FALSE)

#--------------------------------------------------------------------------------#
pa.merge.filtered <- pa.merge
for (i in 1:length(markerListCombos)){
  geneCombos <- expand.grid(markerListCombos[[i]][[1]], markerListCombos[[i]][[2]])
  for (i in 1:length(row.names(geneCombos))){
    expr <-  FetchData(pa.merge.filtered, geneCombos$Var1[i])
    expr2 <- FetchData(pa.merge.filtered, geneCombos$Var2[i])
    pa.merge.filtered<- pa.merge.filtered[,which(x = !(expr >0 & expr2>0) )]
  }
}
pa.merge.filtered
#--------------------------------------------------------------------------------#
pa.merge.filtered <- subset(pa.merge.filtered, EPCAM >0 & PTPRC>0, invert = TRUE)

ggpubr::ggarrange(
  VlnPlot(pa.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = 'orig.ident'),
  VlnPlot(pa.merge.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = 'orig.ident'),
  nrow = 2)


p1<-ggpubr::ggarrange(VlnPlot(pa.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = 'orig.ident'))
p2<-ggpubr::ggarrange(VlnPlot(pa.merge.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = 'orig.ident'))
library(ggpubr)
p1<-annotate_figure(p1, top = text_grob("Pre-filtering", face="bold", size=20))
p2<-annotate_figure(p2, top = text_grob("Post-filtering", face="bold", size=20))
outPlot<-ggpubr::ggarrange(
  p1,
  p2,
  nrow = 2)
