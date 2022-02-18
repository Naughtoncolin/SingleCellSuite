library(Seurat)
library(dplyr)
library(ggplot2)

options(bitmapType="cairo")
paList <- readRDS("~/Dropbox (GaTech)/Gibson/Working/Perianal Fistula/paSamples_210118.RDS")
pa5 <- paList$PA_5

df <- data.frame(pa5$nuclear_fraction, pa5$nCount_RNA)
fit <- kmeans(df, 4)
df <- data.frame(df, fit$cluster)
ggplot(df, aes(x=pa5.nuclear_fraction, y=log10(pa5.nCount_RNA), colour=fit.cluster)) + geom_point()

rm(fit)
df <- data.frame(pa5$nuclear_fraction, log10(pa5$nCount_RNA))
fit <- kmeans(df, 3)
aggregate(df,by=list(fit$cluster),FUN=mean)
df <- data.frame(df, fit$cluster)
ggplot(df, aes(x=pa5.nuclear_fraction, y=log10.pa5.nCount_RNA., colour=fit.cluster)) + geom_point()




df <- data.frame(gbm$nuclear_fraction, log10(gbm$nCount_RNA))
#pamClass<-pam(df, 3)
#df <- cbind(df, cluster=pamClass$clustering)
#ggplot(df, aes(x=pa5.nuclear_fraction, y=log10.pa5.nCount_RNA., colour=cluster)) + geom_point()
##Gaussian Mixed Model Clustering
dat <- center_scale(df, mean_center = T)
gmm = GMM(dat, gaussian_comps = 3)
pr = predict_GMM(dat, gmm$centroids, gmm$covariance_matrices, gmm$weights)
pr
pr$cluster_labels
blah<-pr$cluster_labels
df <- cbind(df, GMMclust=blah)
namechange <- c("1"="cell", "2"='damaged_cell', "3"='empty_droplet')
df$GMMclust <- as.character(namechange[df$GMMclust])
ggplot(df, aes(x=pa6.nuclear_fraction, y=log10.pa6.nCount_RNA., colour=GMMclust)) + geom_point()
opt_gmm = Optimal_Clusters_GMM(dat, max_clusters = 10 )

## GMM-based clustering 
library(ClusterR)
slimGMM <-function(x, components){
  # Make data frame from Seurat object metadata
  #df <- data.frame(x$nuclear_fraction_droplet_qc, x$log10_umi_count)
  df <- data.frame(x$nuclear_fraction, log10(x$nCount_RNA))
  #df <- data.frame(x$nuclear_fraction, log10(x$nCount_RNA), x$percent.mt)
  #pamClass<-pam(df, 3)
  #df <- cbind(df, cluster=pamClass$clustering)
  #ggplot(df, aes(x=pa5.nuclear_fraction, y=log10.pa5.nCount_RNA., colour=cluster)) + geom_point()
  ##Gaussian Mixed Model Clustering
  dat <- center_scale(df, mean_center = T)
  gmm = GMM(dat, gaussian_comps = components, dist_mode = 'maha_dist')
  pr = predict_GMM(dat, gmm$centroids, gmm$covariance_matrices, gmm$weights)
  pr$cluster_labels
  blah<-pr$cluster_labels
  df <- cbind(df, GMMclust=blah)
  #namechange <- c("1"="cell", "2"='damaged_cell', "3"='empty_droplet', '4'='unknown') #Commented out b/c different samples got different numbers for each category
  #namechange <- c("3"="cell", "1"='damaged_cell', "2"='empty_droplet')
  #namechange <- c("1"='damaged_cell',"2"='empty_droplet',"3"="cell") # Using "'1'='damaged_cell'" doesn't actually set 1 to 'damaged_cell' if it's in an index other than 1
  #df$GMMclust <- as.character(namechange[df$GMMclust])#Commented out b/c different samples got different numbers for each category
  x$GMMclust <- df$GMMclust
  ggplot(df, aes(x=x.nuclear_fraction, y=log10.x.nCount_RNA, colour=DropletQC)) + geom_point() +ggtitle(x$'orig.ident')
  return(factor(df$GMMclust)) #Uncomment if adding metadata to seurat object
  xTemp <- x
  xTemp$GMMclust <- 'original'
  xmerge <-merge(x,xTemp)
  xmerge$GMMclust <- factor(xmerge$GMMclust , levels = c("original", "cell", "damaged_cell", "empty_droplet"))
  print("GMM")
  print(table(df$GMMclust))
  print("DropletQC")
  print(table(x$DropletQC))
  #o1<- VlnPlot(xmerge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "DropletQC")
  #o2<- VlnPlot(xmerge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "GMMclust")
  #VlnPlot(xmerge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'nuclear_fraction'), ncol = 3, pt.size =0, group.by = "GMMclust")
  
  #ggplot(df, aes(x=x.nuclear_fraction, y=log10.x.nCount_RNA., colour=GMMclust)) + geom_point()
  ggplot(df, aes(x=x.nuclear_fraction, y=log10.x.nCount_RNA, colour=DropletQC)) + geom_point() +ggtitle(x$'orig.ident')
  #opt_gmm = Optimal_Clusters_GMM(dat, max_clusters = 10 )
  
}
slimGMM(dqc$PBMC, 2)
slimGMM(gbm, 3)
slimGMM(paList$PA_5,3)
length(colnames(paList$PA_5)[paList$PA_5$DropletQC=='cell'])

pDQC<- VlnPlot(pa5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'nuclear_fraction'), ncol = 3, pt.size =0, group.by = "DropletQC")
pGMM<- VlnPlot(pa5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'nuclear_fraction'), ncol = 3, pt.size =0, group.by = "GMMclust")
pDQC<- VlnPlot(pa5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "DropletQC")
pGMM<- VlnPlot(pa5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size =0, group.by = "GMMclust")
ggpubr::ggarrange(pGMM, pDQC,nrow = 2)

pa5$GMMclust

names(gbm)[names(gbm)=="umi_count"] <- "nCount_RNA"
names(gbm)[names(gbm)=="nuclear_fraction_droplet_qc"] <- "nuclear_fraction"
names(gbm)[names(gbm)=="percent_mt"] <- "percent.mt"
names(gbm)[names(gbm)=="flag"] <- "DropletQC"

#Elbow plot for setermining epsilon for dbscan
plot(rev(sort(kNNdist(df, 3))))
epsilonPar <- function(seurObj, k){
  plot(rev(sort(kNNdist(data.frame(seurObj$nuclear_fraction, log10(seurObj$nCount_RNA)), k))))
}
epsilonPar(paList$PA_7, 25)

#### Density based clustering
#x: Seurat object with nuclear_fraction metadata
#y: Radius parameter
#z: Minimum number of points within radius
# Returns plot
library(dbscan)
denScan <-function(x, y, z){
  #Add portion that incluses minpts and epsilon in plot
  df <- data.frame(x$nuclear_fraction, log10(x$nCount_RNA))
  dat <- center_scale(df, mean_center = T)
  #print(head(df))
  #print(head(dat))
  dbtest <- dbscan(dat, eps = y, minPts = z)
  print(dbtest)
  df <- cbind(df, dcscan=factor(dbtest$cluster))
  ggplot(df, aes(x=x.nuclear_fraction, y=log10.x.nCount_RNA., colour=dcscan)) + geom_point()
}

denScan(gbm,.2,25)
denScan(paList$PA_7,.15,24)

multiDenScan <- function(x, y, z){
  #Add parameter to return list of plots
  plotList <- lapply(x, function(x){
    denScan(x, y, z)
  })
  multiPlot <- ggpubr::ggarrange(ggpubr::ggarrange(plotList[[1]], plotList[[2]], ncol=2),
                                 ggpubr::ggarrange(plotList[[3]], plotList[[4]], ncol=2), nrow = 2)
  return(multiPlot)
}
multiDenScan(dqc, .2, 25)



#hdbscan
#seurObj: Seurat Object
#k: k neightbors
hd <- function(seurObj, k){
  df <- data.frame(nuclear_fraction=seurObj$nuclear_fraction, log10.nCount_RNA=log10(seurObj$nCount_RNA))
  print(hdbscan(df, k))
  cl<- hdbscan(df, k)
  print(colnames(df))
  return(ggplot(df, aes(x=nuclear_fraction, y=log10.nCount_RNA, colour=factor(cl$cluster))) + geom_point())
  #plot(df , col=cl$cluster+1, pch=20)
}
hd(dqc$MB,25)
plotList <- lapply(dqc, function(x){
  
})
library(DropletQC)
data('qc_examples')
# Parse DropletQC vignette data
dqc <- split(qc_examples, qc_examples$sample) #Could use filter if only wanting 1 element
dqc<- lapply(dqc, function(x){
  names(x)[names(x)=="umi_count"] <- "nCount_RNA"
  names(x)[names(x)=="nuclear_fraction_droplet_qc"] <- "nuclear_fraction"
  names(x)[names(x)=="percent_mt"] <- "percent.mt"
  names(x)[names(x)=="flag"] <- "DropletQC"
  return(x)
  #split(x, x$cell_type)
})

lapply(dqc, function{x}{
  lapply(x, function(y){
    slimGMM(y, 3)
  })
})

lapply(dqc$PBMC, function(y){
  slimGMM(y, 3)
})
