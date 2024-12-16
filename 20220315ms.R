
# 2022 MAR
# Yuying Rong

setwd('~/Labs/xinhua/massSpec/')
x <- read.csv('iFOT_20200913.csv')
x[1,]

# remove first two cols of x, which are Gene_Symbol and Gene_ID
# the remaining cols are all numericals, are assigned to a new var y
y <- x[3:24]

# colnames of y are now Sample_1, etc., change them to numbers, shorter for plot label
colnames(y) <- 1:22

# add gene name info from x to rownames of y
rownames(y) <- x[,1]
y[1,]

# because y is a data.frame, plotting looks weird
# coerce it into a matrix and assign to a new var z
z <- as.matrix(y)

plot(z)
# this plot is just plotting the values of the first two samples
# this mainly tests if z can be heatmapped


pdf(file='MS_plots.pdf',width=10,height=10)
### hierarchical clustering (dendrogram) and heatmap
heatmap(z)


### k-means clustering

# tutorial used: https://techvidvan.com/tutorials/cluster-analysis-in-r

# want the samples as rownames and parameters as colnames, take transpose: t(z)
# start with 3 centers, so will get 3 clusters
# start with 10 randomly sowed centers: nstart=10

clust <- kmeans(x=t(z),centers=3,nstart=10)
# same result when nstart=10 and =50

# try 5 clusters
clust2 <- kmeans(x=t(z),centers=5,nstart=10)
# same result when nstart=10 and =50, so kept 10

# ran other tests; found <3 and >5 clusters not informative
# centers<3: 2 clusters...
# centers>5: 6 clusters with roughly same num of items in each cluster

# view num of items in each cluster
table(clust2$cluster)
# 1 2 3 4 5 
# 1 9 2 5 5

# view detailed cluster-sample correlation
clust2$cluster
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 
# 1  1  4  3  4  1  1  1  4  3  3  2  4  2  4  3  3  1  1  5  1  1 

# view in order
sort(clust2$cluster)
# 1  2  6  7  8 18 19 21 22  3  5  9 13 15 12 14  4 10 11 16 17 20 
# 1  1  1  1  1  1  1  1  1  2  2  2  2  2  3  3  4  4  4  4  4  5 

# may use fviz_cluster to view nice plots, but it is clear enough for now


### PCA dimension reduction

# followed this tutorial for PCA: https://www.datacamp.com/community/tutorials/pca-analysis-r

# want samples as rownames, and parameters as colnames: take transpose t(z)
# set center point and scaling to TRUE
pca_z <- prcomp(t(z),center=TRUE,scale=TRUE)

# install new package to view nice plots
install.packages('devtools')
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

ggbiplot(pca_z, var.axes=FALSE, labels=colnames(z))

# sep by fast slow muscles
# 'IIb': 1,3,4,7,8,13,15,17,18,19
# 'I+IIx': 5,6,9,10,11,12,14,16
# 'IIx': 20,21,22
# 'I': 2

# pass categorical values (classifiers) to group option
z.type <- c('IIb','I','IIb','IIb','I+IIx','I+IIx','IIb','IIb','I+IIx','I+IIx','I+IIx','I+IIx','IIb','I+IIx','IIb','I+IIx','IIb','IIb','IIb','IIx','IIx','IIx')

ggbiplot(pca_z, ellipse=TRUE, labels=colnames(z), groups=z.type, var.axes=FALSE)

# find most highly variable features by stdev
sort(pca_z$scale,decreasing=TRUE)[1:10]

# however, these highly variable features also have the highest mean expression levels
sort(pca_z$center,decreasing=TRUE)[1:10]

# plot the stdev of the first 100 points with the highest stdev
plot(sort(pca_z$scale,decreasing=TRUE)[1:100])
# the stdev sharply decreases at ~20

summary(pca_z)
str(pca_z)

# view PC values of each sample
pca_z$x[1:5,1:3]


library(dplyr)

# plot PCA result
pca_z$x %>% as.data.frame %>% ggplot(aes(x=PC1,y=PC2))+geom_point(size=1)

# or use this
plot(pca_z$x)
# which is the same as
plot(x=pca_z$x[,1],y=pca_z$x[,2])
# but I will stick to the simpler for now
# plot with labels
plot(pca_z$x)
text(pca_z$x,colnames(z))
# plot labels with offsets
plot(pca_z$x)
text(pca_z$x,colnames(z),pos=1)

# find the most variable feature by rotation value

# find the most variable feature in each PC
# rotation, i think, is positively related to the contribution of a feature to the PC; has the largest eigenvalue
# inspired by: https://stats.stackexchange.com/questions/108148/pull-out-most-important-variables-from-pca
# and by: https://stackoverflow.com/questions/43811572/pca-in-r-how-to-determine-contribution-of-each-variable-to-a-pc-score

# look at the rotation score for each PC
length(pca_z$rotation[1,])
# 22 PCs
length(pca_z$rotation[,1])
# returns 3887; too long

# just view genes with the highest rotation val
sort(pca_z$rotation[,1],decreasing=TRUE)[1:10]
pc1 <- cbind(rownames(z), pca_z$rotation[,1])
# these are the rownames
sort(pc1[,1],decreasing=TRUE)[1:10]
# these are ordered rotation values, but they are coerced into strings
sort(as.numeric(pc1[,2]),decreasing=TRUE)[1:10]
# but this does not return rownames
# I'll figure out this method later, but use easier one for now

# view most negative rotations for PC1, PC2
sort(sort(pca_z$rotation[,1],decreasing=TRUE)[3878:3887],decreasing=FALSE)
# nah, that was stupid
sort(pca_z$rotation[,1],decreasing=FALSE)[1:10]
sort(pca_z$rotation[,2],decreasing=FALSE)[1:10]
# view most positive rotations for PC1, PC2
sort(pca_z$rotation[,1],decreasing=TRUE)[1:10]
sort(pca_z$rotation[,2],decreasing=TRUE)[1:10]

# getting the original indices from ordered array
# using index.return flag, which produces index list $ix in addition to content list $x
# retrieving only the index list with $ix
sort(pca_z$rotation[,1],decreasing=TRUE,index.return=TRUE)$ix[1:10]

# found that order() will just return indices
order(pca_z$rotation[,1],decreasing=TRUE)[1:10]
sort(pca_z$rotation[,1],decreasing=TRUE,index.return=TRUE)$ix[1:10]

# view heatmap of just top10 values

# using sort() with index.return
heatmap(z[c(
	sort(pca_z$rotation[,1],decreasing=TRUE,index.return=TRUE)$ix[1:10],
	sort(pca_z$rotation[,1],decreasing=FALSE,index.return=TRUE)$ix[1:10]
	),])

# using order()
heatmap(z[c(
	order(pca_z$rotation[,1],decreasing=TRUE)[1:10],
	order(pca_z$rotation[,1],decreasing=FALSE)[1:10]
	),])

heatmap(z[c(
	order(pca_z$rotation[,2],decreasing=TRUE)[1:10],
	order(pca_z$rotation[,2],decreasing=FALSE)[1:10]
	),])

heatmap(z[c(
	order(pca_z$rotation[,1],decreasing=TRUE)[1:10],
	order(pca_z$rotation[,1],decreasing=FALSE)[1:10],
	order(pca_z$rotation[,2],decreasing=TRUE)[1:10],
	order(pca_z$rotation[,2],decreasing=FALSE)[1:10]
	),])



# I also plotted them with plot() and looked at $ratation[,2]

# pc1 <- cbind(rownames(z), pca_z$rotation[,1])
# pc1[1,]
# wanted to add a col of gene names, but this worked poorly
# the numbers became strings

### kmeans clustering after PCA reduction

# inspired by this post: https://stats.stackexchange.com/questions/475883/dataset-transformation-for-clustering-after-pca

# view proportion of variance
summary(pca_z)$importance

# it seems each PC has quite some importance, not great
# say I want 90% importance accounted
# then I will make a new dataset with 1:14 PCs
z2=pca_z$x[,1:14]

k_z <- kmeans(x=z2, centers=3, nstart=20)
k_z$cluster
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 
# 2  2  2  2  2  2  2  2  2  2  2  3  2  1  2  2  2  2  2  2  2  2 
plot(z2, col=k_z$cluster, cex=0.1)
plot(z2, col=k_z$cluster)

k_z <- kmeans(x=z2, centers=5, nstart=20)
sort(k_z$cluster)
plot(sort(k_z$cluster))
plot(z2, col=k_z$cluster)


dev.off()

