slimGMM.plot <-function(x, components){
  # Make data frame from Seurat object metadata
  df <- data.frame(x$nuclear_fraction_droplet_qc, x$log10_umi_count)
  #df <- data.frame(x$nuclear_fraction, log10(x$nCount_RNA))
  #df <- data.frame(x$nuclear_fraction, log10(x$nCount_RNA), x$percent.mt)
  #pamClass<-pam(df, 3)
  #df <- cbind(df, cluster=pamClass$clustering)
  #ggplot(df, aes(x=pa5.nuclear_fraction, y=log10.pa5.nCount_RNA., colour=cluster)) + geom_point()
  ##Gaussian Mixed Model Clustering
  dat <- center_scale(df, mean_center = T)
  gmm = GMM(dat, gaussian_comps = components, dist_mode = 'maha_dist')
  pr = predict_GMM(dat, gmm$centroids, gmm$covariance_matrices, gmm$weights)
  #pr$cluster_labels
  blah<-pr$cluster_labels
  df <- cbind(df, GMMclust=blah)
  
  #namechange <- c("1"="cell", "2"='damaged_cell', "3"='empty_droplet', '4'='unknown')
  #namechange <- c("3"="cell", "1"='damaged_cell', "2"='empty_droplet')
  #namechange <- c("1"='damaged_cell',"2"='empty_droplet',"3"="cell") # Using "'1'='damaged_cell'" doesn't actually set 1 to 'damaged_cell' if it's in an index other than 1
  #df$GMMclust <- as.character(namechange[df$GMMclust])
  x$GMMclust <- df$GMMclust
  xTemp <- x
  xTemp$GMMclust <- 'original'
  xmerge <-merge(x,xTemp)
  xmerge$GMMclust <- factor(xmerge$GMMclust , levels = c(0,1,2,3))
  print("GMM")
  print(table(df$GMMclust))
  print("DropletQC")
  print(table(x$DropletQC))
  #o1<- VlnPlot(xmerge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "GMMclust")
  #o2<- VlnPlot(xmerge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "GMMclust")
  #VlnPlot(xmerge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'nuclear_fraction'), ncol = 3, pt.size =0, group.by = "GMMclust")
  #ggplot(df, aes(x=x.nuclear_fraction, y=log10.x.nCount_RNA., colour=GMMclust)) + geom_point()
  df$GMMclust <- factor(df$GMMclust)
  return(ggplot(df, aes(x=x.nuclear_fraction_droplet_qc, y=x$log10_umi_count, colour=GMMclust)) + geom_point()) #+ggtitle(x$sample)
  #opt_gmm = Optimal_Clusters_GMM(dat, max_clusters = 10 )
}
p1 <- slimGMM.plot(dqc$PBMC, 2)
pList <- lapply(dqc, function(x){
  p1<-slimGMM.plot(x, 2) +ggtitle(x$sample)
  p2<-slimGMM.plot(x, 3)
  p.grid <- ggpubr::ggarrange(p1, p2, nrow = 2)
})
pList
p.grid <- ggpubr::ggarrange(pList$GBM, pList$HL, pList$MB, pList$PBMC,ncol = 4)
p.grid


#Scatterplot metadata
p2List <lapply(dqc, function(x){
  ggplot(x, aes(x=x.nuclear_fraction_droplet_qc, y=x$log10_umi_count, colour=GMMclust)) + geom_point()
})

library(dbscan)
denScan <-function(x, y, z){
  df <- data.frame(x$nuclear_fraction_droplet_qc, x$log10_umi_count)
  dat <- center_scale(df, mean_center = T)
  #print(head(df))
  #print(head(dat))
  dbtest <- dbscan(dat, eps = y, minPts = z)
  #print(dbtest)
  df <- cbind(df, dcscan=dbtest$cluster)
  p1<- ggplot(df, aes(x=x.nuclear_fraction_droplet_qc, y=x$log10_umi_count, colour=factor(dcscan))) + geom_point()
  p2 <- ggplot(x, aes(x=x$nuclear_fraction_droplet_qc, y=x$log10_umi_count, colour=x$flag)) + geom_point()
  print(table(df$dcscan))
  print(table(x$flag))
  ggpubr::ggarrange(p1, p2, nrow = 2)
  
}
denScan(dqc$GBM, 0.3, 100)
GBM <- split(dqc$GBM,dqc$GBM$cell_type)

lapply(GBM, function(x){
  denScan(x, 0.3, 100)
})

# Make scatter plot log10UMI-count and nuclear fraction colored by DropletQC cell classification
scatNuc <- function(x){
  ggplot(x, aes(x=nuclear_fraction_droplet_qc, y=log10_umi_count, colour=flag)) + geom_point() + ggtitle(x$sample)
}








