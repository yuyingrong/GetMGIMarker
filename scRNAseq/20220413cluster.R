
# 20220413 clustering
# "as default as possible"

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
#BiocManager::install ("scater")
library(scater)

setwd('~/Labs/xinhua/scSeq/')

data <- readRDS(file='rds/merged.rds')

# script 1 

# check if nCells matches expectation based on nCells loaded and lib prep method used
data@meta.data %>% ggplot(aes(x=sample, fill=sample))+
	geom_bar()+
	theme_classic()+
	theme(axis.text.x=element_text(angle=45,hjust=1))+
	theme(plot.title=element_text(hjust=0.5,face='bold'))+
	ggtitle('nCells')

VlnPlot(data, features = c('nFeature_RNA','nCount_RNA'),ncol=2)
VlnPlot(data, features = c('percent_mt','log10GenesPerUMI'),ncol=2)

# another view of above matrix
data@meta.data %>% ggplot(aes(x=nCount_RNA,y=nFeature_RNA,color=percent_mt))+
	geom_point(alpha=0.5)+
	scale_color_gradient(low='gray',high='black')+
	stat_smooth(method=lm)+
	scale_x_log10()+
	scale_y_log10()+
	theme_classic()+
	geom_vline(xintercept=500)+
	geom_hline(yintercept=250)+
	facet_wrap(~sample)

# eyeballed QC cutoff
data <- subset(data, subset=(nFeature_RNA>200 & percent_mt<50))

data <- NormalizeData(data, normalization.method='LogNormalize', scale.factor=10000)

data <- FindVariableFeatures(data, selection.method='vst', nfeatures=2000)
# VariableFeatures(data) <- SelectIntegrationFeatures(SplitObject(data, split.by='sample'), nfeatures = 2000, verbose = TRUE, fvf.nfeatures = 2000, selection.method = "vst")

data <- ScaleData(data)

# linear dimensional reduction
data <- RunPCA(data, features=VariableFeatures(object=data))

DimPlot(data, reduction='pca')

ElbowPlot(data, ndim=50, reduction='pca')

DimHeatmap(data, dims=1:15, cells=500, balanced=TRUE)
DimHeatmap(data, dims=16:30, cells=500, balanced=TRUE)
DimHeatmap(data, dims=31:45, cells=500, balanced=TRUE)

# cluster cells
data <- FindNeighbors(data, dims=1:20)
# used 1:40 for merged

data <- FindClusters(data, resolution=0.1)

# non-linear dimensional reduction
data <- RunUMAP(data, dims=1:20)
# used 1:40 for merged

DimPlot(data, reduction='umap')

DimPlot(data,reduction='umap',label=TRUE,split.by='sample')

VlnPlot(data, features=c('PAX7'))
VlnPlot(data, features=c('SCX'))

FeaturePlot(data, features=c('PAX7', 'SCX'))

FeaturePlot(data, features=c('PAX7', 'SCX', 'MYOG', 'MYOD1', 'MYF5'), ncol=3)

FeaturePlot(data, features=c('PAX7'), split.by='sample')

length(which(GetAssayData(data,slot='data')['PAX7',]>0))
# 3834
length(GetAssayData(data,slot='data')['PAX7',])
# 6285


# script 2, with integration

data <- readRDS(file='rds/merged.rds')

data@meta.data[1,]
levels(factor(data@meta.data$cluster))
levels(factor(data@meta.data$sample))

# look at cluster 1-6
data <- subset(data, subset=(cluster!='clust7' & cluster!='unassigned'))
# look at sample 1-2, cluster 1-6
data <- subset(data, subset=(sample=='fx1' | sample=='fx2'))
data <- subset(data, subset=(cluster!='clust7' & cluster!='unassigned'))
# look at sample 1/3/4/5, cluster 1-6
data <- subset(data, subset=(sample=='fx1' | sample=='fx3' | sample=='fx4' | sample=='fx5'))
data <- subset(data, subset=(cluster!='clust7' & cluster!='unassigned'))
# look at sample 1/3/4/5/6, cluster 1-6
data <- subset(data, subset=(sample!='fx2'))
data <- subset(data, subset=(cluster!='clust7' & cluster!='unassigned'))

# pdf(file='sample_1to2_cluster_1to6.pdf',width=10,height=10)
# pdf(file='sample_1345_cluster_1to6.pdf',width=10,height=10)
# pdf(file='sample_13456_cluster_1to6.pdf',width=10,height=10)

data@meta.data %>% ggplot(aes(x=nCount_RNA,y=nFeature_RNA,color=percent_mt))+
	geom_point(alpha=0.5)+
	scale_color_gradient(low='gray',high='black')+
	stat_smooth(method=lm)+
	scale_x_log10()+
	scale_y_log10()+
	theme_classic()+
	geom_vline(xintercept=500)+
	geom_hline(yintercept=250)+
	facet_wrap(~sample)

# eyeballed QC cutoff
data <- subset(data, subset=(nFeature_RNA>300 & percent_mt<50))

# dimension reduction and clustering
data <- NormalizeData(data, normalization.method='LogNormalize', scale.factor=10000)
data <- FindVariableFeatures(data, selection.method='vst', nfeatures=2000)
data <- ScaleData(data)

data <- RunPCA(data, features=VariableFeatures(object=data))

DimPlot(data, reduction='pca', group.by='sample')
DimPlot(data, reduction='pca', group.by='cluster')
ElbowPlot(data, ndim=50, reduction='pca')

DimHeatmap(data, dims=1:15, cells=500, balanced=TRUE)
DimHeatmap(data, dims=16:30, cells=500, balanced=TRUE)
DimHeatmap(data, dims=31:45, cells=500, balanced=TRUE)
DimHeatmap(data, dims=36:50, cells=500, balanced=TRUE)

data <- FindNeighbors(data, dims=1:50)
# cluster1:6: 1:43
# sample1:2, cluster1:6: 1:30
# sample1345, cluster1:6: 1:50

data <- FindClusters(data, resolution=0.1)
data <- FindClusters(data, resolution=1)

data <- RunUMAP(data, dims=1:50)

DimPlot(data, reduction='umap', label=TRUE, group.by='seurat_clusters')
DimPlot(data,reduction='umap',label=TRUE,split.by='sample', ncol=3)
DimPlot(data,reduction='umap',label=TRUE,group.by='sample')
DimPlot(data,reduction='umap',label=TRUE,split.by='cluster', ncol=2)
DimPlot(data,reduction='umap',label=TRUE,group.by='cluster')

VlnPlot(data, features=c('PAX7'))
VlnPlot(data, features=c('MYF5'))
VlnPlot(data, features=c('MYOD1'))

FeaturePlot(data, features=c('PAX7', 'MYOG', 'MYOD1', 'MYF5'), ncol=2)
FeaturePlot(data, features=c('PAX7'), split.by='cluster', ncol=3)

dev.off()

saveRDS(data, file='rds/sample_all_cluster_1to6.rds')


###################################
###################################
###################################
# add integration

data <- NULL
data <- readRDS(file='rds/merged.rds')

data@meta.data[1,]
levels(factor(data@meta.data$cluster))
levels(factor(data@meta.data$sample))

# look at cluster 1-6
data <- subset(data, subset=(cluster!='clust7' & cluster!='unassigned'))
# look at sample 1-2, cluster 1-6
data <- subset(data, subset=(sample=='fx1' | sample=='fx2'))
data <- subset(data, subset=(cluster!='clust7' & cluster!='unassigned'))
# look at sample 1/3/4/5, cluster 1-6
data <- subset(data, subset=(sample!='fx2' & sample!='fx6'))
data <- subset(data, subset=(cluster!='clust7' & cluster!='unassigned'))
# look at sample 1/3/4/5/6, cluster 1-6
data <- subset(data, subset=(sample!='fx2'))
data <- subset(data, subset=(cluster!='clust7' & cluster!='unassigned'))
# look at sample 2, cluster 1-6
data <- subset(data, subset=(sample=='fx2'))
data <- subset(data, subset=(cluster!='clust7' & cluster!='unassigned'))
# look at sample 1-3, cluster 1-6
data <- subset(data, subset=(cluster=='clust1' | cluster=='clust2' | cluster=='clust3'))
# look at sample 1-5, cluster 1-6
data <- subset(data, subset=(sample!='fx6'))
data <- subset(data, subset=(cluster!='clust7' & cluster!='unassigned'))
# look at sample 1-6, cluster 4-6
data <- subset(data, subset=(cluster=='clust4' | cluster=='clust5' | cluster=='clust6'))
# look at sample 1-6, cluster 1-3
data <- subset(data, subset=(cluster=='clust1' | cluster=='clust2' | cluster=='clust3'))
# look at sample 1/3/4/5/6, cluster 1/3/4/6
data <- subset(data, subset=(sample!='fx2'))
data <- subset(data, subset=(cluster=='clust1' | cluster=='clust3' | cluster=='clust4' | cluster=='clust6'))



# QC
data <- subset(data, subset=(nFeature_RNA>300 & percent_mt<50))
data <- subset(data, subset=(nFeature_RNA>500 & percent_mt<8))
# caution! this will work incorrectly
# data <- subset(data, subset=(sample==paste0('fx', c(1:4))))

# split the dataset into a list of two seurat objects
data <- SplitObject(data, split.by='sample')

# normalize and identify variable features for each dataset independently
data <- lapply(X=data, FUN=function(x) {
	x <- NormalizeData(x)
	x <- FindVariableFeatures(x, selection.method='vst', nfeatures=2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list=data)

# identify anchors
anchors <- FindIntegrationAnchors(object.list=data, anchor.features=features)
rm(features)

# use anchors to integrate two datasets together
data <- IntegrateData(anchorset=anchors)
rm(anchors)

# specify to use anchored data for downstream analysis
DefaultAssay(object=data) <- 'integrated'

# standard clustering visualization workflow
data <- ScaleData(data)
data <- RunPCA(data)

DimPlot(data, reduction='pca', label=TRUE, group.by='sample')
DimPlot(data, reduction='pca', label=TRUE, group.by='cluster')
ElbowPlot(data, ndim=50, reduction='pca')

DimHeatmap(data, dims=1:15, cells=500, balanced=TRUE)
DimHeatmap(data, dims=16:30, cells=500, balanced=TRUE)
DimHeatmap(data, dims=31:45, cells=500, balanced=TRUE)
DimHeatmap(data, dims=36:50, cells=500, balanced=TRUE)

data <- RunUMAP(data, dims=1:50)

data <- FindNeighbors(data, dims=1:50)
# typically 0.4-1.2 for single-cell datasets of around 3K cells
data <- FindClusters(data, resolution=0.55)
data <- FindClusters(data, resolution=0.1)
data <- FindClusters(data, resolution=1)

# pdf(file='sample_1to6_cluster_1to6_integ.pdf',width=10,height=10)
# pdf(file='sample_1to2_cluster_1to6_integ.pdf',width=10,height=10)
# pdf(file='sample_1345_cluster_1to6_integ.pdf',width=10,height=10)
# pdf(file='sample_13456_cluster_1to6_integ.pdf',width=10,height=10)
# pdf(file='sample_2_cluster1to6.pdf',width=10,height=10)
# pdf(file='sample_all_cluster_1to3_integ.pdf',width=10,height=10)
# pdf(file='sample_1to5_cluster_1to6_integ.pdf',width=10,height=10)
# pdf(file='sample_1to6_cluster_4to6_integ.pdf',width=10,height=10)
pdf(file='sample_13456_cluster_1346_integ.pdf',width=10,height=10)
pdf(file='add_sample_13456_cluster_1346.pdf',width=10,height=10)

DimPlot(data, reduction='umap', label=TRUE, group.by='seurat_clusters')
DimPlot(data,reduction='umap',label=TRUE,split.by='sample', ncol=3)
DimPlot(data,reduction='umap',label=TRUE,group.by='sample')
DimPlot(data,reduction='umap',label=TRUE,split.by='cluster', ncol=2)
DimPlot(data,reduction='umap',label=TRUE,group.by='cluster')

DefaultAssay(object=data) <- 'RNA'

data@meta.data %>% ggplot(aes(x=nCount_RNA,y=nFeature_RNA,color=percent_mt))+
	geom_point(alpha=0.5)+
	scale_color_gradient(low='gray',high='black')+
	stat_smooth(method=lm)+
	scale_x_log10()+
	scale_y_log10()+
	theme_classic()+
	geom_vline(xintercept=500)+
	geom_hline(yintercept=250)+
	facet_wrap(~sample)

VlnPlot(data, features='percent_mt')
VlnPlot(data, features='nFeature_RNA')
VlnPlot(data, features=c('PAX7'))
VlnPlot(data, features=c('MYF5'))
VlnPlot(data, features=c('MYOD1'))

FeaturePlot(data, features=c('PAX7', 'MYOG', 'MYOD1', 'MYF5'), ncol=2)
FeaturePlot(data, features=c('PAX7'), split.by='cluster', ncol=3)

FeaturePlot(data, features=c('percent_mt','nFeature_RNA'))
#FeaturePlot(data, features='nCount_RNA')

dev.off()

str(data)

saveRDS(data, file='rds/sample_all_cluster_1to6_integ.rds')
saveRDS(data, file='rds/sample_1to2_cluster_1to6_integ.rds')

###################################
###################################
###################################
# marker visualization

data <- readRDS(file='rds/sample_all_cluster_1to6_integ.rds')

DefaultAssay(object=data) <- 'RNA'

# pdf(file='sample_all_cluster_1to6_integ_plots.pdf',width=10,height=10)
pdf(file='sample_13456_cluster_1to6_integ_plots.pdf',width=10,height=10)

# genes to plot
# use Bash Shell to paste genes copied from Excel; or just paste into a TXT
# use Python to format string:
#   f=open('temp2.txt','r')
#   s2=f.read()
#   re.sub('\n',"','",s2)
# use toupper() to capitalize
# did not use use levels(factor()) to remove duplicates, 
#   because the gene names are changed to alphabetical order
# to retain original order, use unique()
dot_plot_genes <- toupper(unique(c('PAX7','MYF5','MYOD1','NCAM1','ITGB1','CXCR4','VCAM1','ICAM1','CD82','DLK1','CHRDL2','NDRG2','DAG1','CHRNA1','SPRY1','HES1','HEY1','EGR1','CAV1','CHODL','CYCS','KLF4','MYC','CHRDL2','FOS','LPL','CALCR','FGFR4','APOD','APOE','LICAM','DLK1','CDH15','APOC1','MYOG','MYF6','MYH2','MEF2A','MEF2C')))

length(dot_plot_genes)
# 37

# plot DotPlot
# use ggplot2's coord_flip() to flip x- and y-axes
DotPlot(object=data, features=dot_plot_genes) + coord_flip()
# added some interesting genes
DotPlot(object=data, features=toupper(unique(c('CDKN1A','CDKN1B','CDKN1C','TP53','RBL2','CCND1','CCND3','MKI67','ANLN','BIRC5','CCNA2','CCNB1','CCNE2','CGOL1','CENPE','TOP2A','UBE2C','CDk1','CDK2','CDC6','CDC20')))) + coord_flip()

# good, but still too many genes on the y-axis
# use a loop


# use levels(factor()) to remove duplicates
feature_plot_genes <- toupper(unique(c('PAX7','MYF5','MYOD1','NCAM1','ITGB1','CXCR4','VCAM1','ICAM1','CD82','DLK1','CHRDL2','NDRG2','DAG1','CHRNA1','SPRY1','HES1','HEY1','EGR1','CAV1','CHODL','CYCS','KLF4','MYC','Sox8','MYOG','MYF6','MYH2','MEF2A','MEF2C','NDRG2','CHRDL2','FOS','LPL','CALCR','FGFR4','APOD','APOE','LICAM','FZD4','DLK1','GNAS','WNT5A','JAK1','JAK2','JAK3','STAT1','STAT2','STAT3','STAT6','HTRA1','HES4','CDH15','IRF1','SOD2','FOXO1','FOXO3','DKK3','CAV2','IGFBP5','IGFBP6','KLF4','IGF1','IGF2','APOC1','KLF14','KLF15','DKK2','IGTA1','JAG1','BMP','ACTA1','CKM','TNNC1','MYLPF','CDKN1A','CDKN1B','CDKN1C','TP53','RBL2','CCND1','CCND3','MKI67','ANLN','BIRC5','CCNA2','CCNB1','CCNE2','CGOL1','CENPE','TOP2A','UBE2C','CDk1','CDK2','CDC6','CDC20')))

#feature_plot_genes <- toupper(levels(factor(c('PAX7','MYF5','MYOD1','NCAM1','ITGB1','CXCR4','VCAM1'))))

# FeaturePlot(object=data, features=feature_plot_genes[1:4], ncol=2)
length(feature_plot_genes)
# 91

# use seq() to count every sep between i and j: seq(i, j, by=sep)
# in loops, auto printing does not work
# need to use print()
# "To drawn lattice plots on the device, one needs to print the object produced by a call to one of the lattice graphics functions. Normally, in interactive use, R auto prints objects if not assigned. In loops however, auto printing does not work, so one must arrange for the object to be printed, usually by wrapping it in print()."

for (i in seq(1, length(feature_plot_genes), by=4)) {
	if ((i+3)>length(feature_plot_genes)) {
		print(FeaturePlot(object=data, features=feature_plot_genes[i:length(feature_plot_genes)], ncol=2))
		break
	}
	print(FeaturePlot(object=data, features=feature_plot_genes[i:(i+3)], ncol=2))
}

dev.off()


# find all cluster 1 markers
# FindMarkers(data, ident.1=1, min.pct=0.25) %>% head(n=50)

# find all markers for each cluster against all other clusters
data.markers <- FindAllMarkers(data, logfc.threshold=0.25, min.pct=0.25, only.pos=FALSE, max.cells.per.ident=Inf)

data.markers[1,]
length(data.markers[,1])
# 15143
data.markers %>% top_n(n=15, wt=avg_log2FC)
data.markers %>% group_by(cluster) %>% top_n(n=8, wt=avg_log2FC)
data.markers %>% group_by(cluster) %>% top_n(n=-8, wt=avg_log2FC)
# the order looks strange
# prove with order()
sub_cluster <- subset(data.markers, cluster==0)
sub_cluster[order(sub_cluster$avg_log2FC, decreasing=FALSE),][1:8,]
# looks right
subset(data.markers, cluster==15) %>% top_n(n=50, wt=avg_log2FC)

# write markers to CSV
for (i in levels(factor(data.markers$cluster))) {
	print(i)
	subset(data.markers, cluster==i) %>% top_n(n=500, wt=avg_log2FC) %>% write.csv(file=paste0('markers/sample_all_clust_1to3/res_0.55/top_n/seurat_clust_',i,'_vs_all.csv'))
	subset(data.markers, cluster==i) %>% top_n(n=-500, wt=avg_log2FC) %>% write.csv(file=paste0('markers/sample_all_clust_1to3/res_0.55/bot_n/seurat_clust_',i,'_vs_all.csv'))
}
# top_n does not seem ranked by what I gave to wt, but by padj

# using order() solves the problem
for (i in levels(factor(data.markers$cluster))) {
	print(i)
	sub_cluster <- subset(data.markers, cluster==i)
	sub_cluster[order(sub_cluster$avg_log2FC, decreasing=TRUE),][1:500,] %>% write.csv(file=paste0('markers/sample_all_clust_1to6/res_1/order_by/seurat_clust_',i,'_vs_all.csv'))
}

saveRDS(data.markers, file='markers/sample_all_clust_1to3/data_markers.rds')
saveRDS(data.markers, file='markers/sample_13456_clust_1346/data_markers.rds')

# how many cells from Seurat clusters are in custom clusters 1-6?
data <- readRDS(file='rds/sample_all_cluster_1to6_integ.rds')

# DefaultAssay(object=data) <- 'integrated'

# levels(factor(data$idents))
# error

output <- matrix(0L, nrow=16, ncol=6)
rownames(output) <- paste0('ident',0:15)
colnames(output) <- paste0('fx',1:6)
for (i in 0:15) {
	print(paste0('Total number of cells in Seurat cluster ',i,':'))
	subcluster <- subset(data, idents=i)
	print(dim(subcluster)[2])
	for (j in 1:6) {
		# beware that i index starts from 0
		output[i+1,j] <- which(subcluster$sample==paste0('fx',j)) %>% length()
	}
}

write.csv(output, 'markers/sample_all_clust_1to3/sample_all_clust_1to3_cell_count.csv')

# cell count with stacked barplot
ggplot(data=data.frame(ident=factor(data$seurat_clusters), y=1:nrow(data@meta.data), sample=factor(data$sample)), aes(x=ident, fill=sample))+
	geom_bar()+
	ggtitle('Cell count per sample in cluster')
# adding count labels to bars:
# https://r-graphics.org/recipe-bar-graph-labels
# https://www.javaer101.com/en/article/15080744.html


# write all markers to CSV
markers <- readRDS(file='markers/sample_13456_clust_1346/data_markers.rds')

for (i in levels(factor(markers$cluster))) {
	print(i)
	sub_cluster <- subset(markers, cluster==i)
	sub_cluster[order(sub_cluster$avg_log2FC, decreasing=TRUE),] %>% write.csv(file=paste0('markers/sample_13456_clust_1346/res_0.55/seurat_clust_',i,'_vs_all.csv'))
}
