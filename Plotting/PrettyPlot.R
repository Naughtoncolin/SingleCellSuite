p.grid <- ggpubr::ggarrange(p1, p2, ggpubr::ggarrange(p3, p3.2, ncol = 2),nrow = 3)

ggpubr::ggarrange(o1, o2, o3,nrow = 3)


ggpubr::ggarrange(ggpubr::ggarrange(pd5, pd6, ncol=2), ggpubr::ggarrange(pd7, pd8, ncol=2),nrow = 2)

pd8<-wrap_plots(d8[[2]]) + plot_annotation(title = "PA_8") #For plotting damaged cell
pd8


ep8<-ggplot(paList$PA_8[[]], aes(x=nuclear_fraction))+
  geom_density()+
  ggtitle("PA_8")

ggplot(paList$PA_8[[]], aes(x=nuclear_fraction))+ 
  geom_density()+
  ggtitle("PA_8")

paMerge <- merge(paList)
metadata<- paMerge@meta.data
library(dplyr)
metadata %>%
  ggplot(aes(color=orig.ident, x=nuclear_fraction, fill=orig.ident, ..density..)) +
  geom_density(alpha = 0.2)


lapply(paList, function(x){
  print(median(x$percent.mt))
  print(median(x$nFeature_RNA))
  print(median(x$nCount_RNA))
  print(median(x$nFeature_RNA)-median(x$nFeature_RNA)*median(x$percent.mt)/100)
  print('')
})

# Plot Dot Plot with x-labels at a 45 degree angle
library(ggplot2)
p <- DotPlot(pa, features = c(Enterocyte.markers, Goblet.markers, StemCell.markers, Immune.CD45pos, Paneth_Like.markers, Enteriendocrine.markers ))
p + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust = 1))

