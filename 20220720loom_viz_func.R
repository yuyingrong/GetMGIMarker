
# 20220720
# pyscenic loom output viz


library(dplyr)
library(Seurat)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(pheatmap)
library(ggplot2)


setwd('Labs/xinhua/scSeq/regulon/')


get_AUC <- function (loom) {

	# pull info from loom
	regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
	regulons <- regulonsToGeneLists(regulons_incidMat)
	regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')

	# get regulon activity info	
	AUCmat <- getAUC(regulonsAUC[,])
	
	# get rid of '(+)' from regulon names, and add the number of genes the regulon controls	
	for (i in 1:length(regulons)) {
		rownames(AUCmat)[i] <- gsub(pattern='(+)', replacement=paste0('_',lengths(regulons)[i],'gene(s)'), x=rownames(AUCmat)[i], fixed=TRUE)
	}
	
	print('get_AUC: done!')
	return(AUCmat)
}



make_heatmap <- function (AUCmat, data) {
	
	# scale AUC matrix for heatmaps
	AUCmat <- t(scale(t(AUCmat)))
	
	# create annotation_col for plotting pheatmap
	# cellInfo <- get_cell_annotation(loom)
	# since sometimes get_cell_annotation() does not work
	cellInfo <- data.frame(sample=data@meta.data$sample)
	cellInfo$group <- factor(data@meta.data$group)
	cellInfo$seurat_cluster <- factor(data@meta.data$seurat_clusters)
	cellInfo$percent_mt <- data@meta.data$percent_mt
	cellInfo$nUMI <- data@meta.data$nCount_RNA
	cellInfo$nGene <- data@meta.data$nFeature_RNA
	rownames(cellInfo) <- rownames(data@meta.data)
	cellInfo[1:4,]
	cellInfo[(nrow(cellInfo)-3):nrow(cellInfo),]

	# binarize results
	dim(AUCmat)
	AUCmat[AUCmat>=2] <- 2
	AUCmat[AUCmat <=(-2)] <- (-2)
	AUCmat[1:4, 1:4]

	# print heatmap
	print(pheatmap(AUCmat, show_colnames=F, annotation_col=cellInfo, main='regulon activity heapmap'))
	print(pheatmap(AUCmat, show_colnames=F, annotation_col=cellInfo, color=c("white","black"), main='Regulon activity heapmap b/w'))
	
	print('make_heatmap: done!')
}



make_custom_heatmap <- function (AUCmat, data) {
	
	#AUCmat <- t(scale(t(AUCmat))) makes no difference

	# get mean regulon AUC per cell, grouped by seurat clusters
	AUC_by_seurat_clusters <- data.frame(
		cluster0=rowMeans(subset(data, seurat_clusters==0)), 
		row.names=rownames(AUCmat))
	print('custom_heatmap: calculating avg regulon activity for seurat cluster: 0')
	
	# create a df to be visualized
	for (i in levels(data$seurat_clusters)[-1]) {
		print(paste0('custom_heatmap: calculating avg regulon activity for seurat cluster: ', i))
		AUC_by_seurat_clusters[[paste0('cluster',i)]] <- rowMeans(subset(data, seurat_clusters==i))
	}
	
	write.csv(AUC_by_seurat_clusters, paste0(obj, '_avgAUC_by_seurat_cluster.csv'))
	
	pheatmap(mat=AUC_by_seurat_clusters, cluster_rows=T, cluster_cols=F, show_colnames=T, annotation_col=NA, main='Mean regulon activity per cell, grouped by seurat clusters')
	
	pheatmap(mat=AUC_by_seurat_clusters, cluster_rows=T, cluster_cols=F, show_colnames=T, annotation_col=NA, main='Mean regulon activity per cell, grouped by seurat clusters blue/red', color=c("blue","red"))


	# get mean regulon AUC per cell, grouped by custom groups
	AUC_by_custom_groups <- data.frame(row.names=rownames(AUCmat))
	
	AUC_by_custom_groups[[levels(factor(data$group))[1]]]=rowMeans(subset(data, group==levels(factor(data$group))[1]))

	for (i in levels(factor(data$group))) {
		print(paste0('custom_heatmap: calculating avg regulon activity for: ', i))
		AUC_by_custom_groups[[i]] <- rowMeans(subset(data, group==i))
	}
	
	write.csv(AUC_by_custom_groups, paste0(obj, '_avgAUC_by_custom_group.csv'))
	
	pheatmap(mat= AUC_by_custom_groups, cluster_rows=T, cluster_cols=F, show_colnames=T, annotation_col=NA, main='Mean regulon activity per cell, grouped by custom groups')

	pheatmap(mat= AUC_by_custom_groups, cluster_rows=T, cluster_cols=F, show_colnames=T, annotation_col=NA, main='Mean regulon activity per cell, grouped by custom groups blue/red', color=c("blue","red"))
	
	print('make_custom_heatmap: done!')
}
	
	

make_feature_plot_from_RNA <- function (regulons, data) {
	
	DefaultAssay(data) <- 'RNA'
	
	features <- gsub(pattern='[(+)]', replacement='', x=names(regulons))
	
	# print regulon feature plots from RNA assay
	for (i in seq(1, length(features), by=4)) {
		if ((i+3)>length(features)) {
			print(FeaturePlot(object=data, features=features[i:length(features)], reduction='pca', label=TRUE, ncol=2))
			print(FeaturePlot(object=data, features=features[i:length(features)], reduction='umap', label=TRUE, ncol=2))
			break
		}
		print(FeaturePlot(object=data, features=features[i:(i+3)], reduction='pca', label=TRUE, ncol=2))
		print(FeaturePlot(object=data, features=features[i:(i+3)], reduction='umap', label=TRUE, ncol=2))
	}
	
	print('make_feature_plot: regulons plotted from RNA clusters!')
}




make_feature_plot_from_AUC <- function (data) {
	
	DefaultAssay(data) <- 'AUC'
	
	# cluster based on AUC matrix
	#data <- NormalizeData(data)
	data <- FindVariableFeatures(data)
	data <- ScaleData(data)
	data <- RunPCA(data)
	data <- RunUMAP(data, dims=1:50)
	data <- FindNeighbors(data, dims=1:50)
	data <- FindClusters(data, resolution=0.9)
	print('viz_loom: clustering based on AUC matrix done!')
	
	# print AUC cluster plots
	print(DimPlot(data, reduction='pca', label=TRUE, group.by='seurat_clusters')+ggtitle('PCA: cells grouped by regulon activity colored by seurat clusters'))
	print(DimPlot(data, reduction='pca', label=TRUE, group.by='sample')+ggtitle('PCA: cells grouped by regulon activity colored by samples'))
	print(DimPlot(data, reduction='pca', label=TRUE, group.by='group')+ggtitle('PCA: cells grouped by regulon activity colored by custom groups'))

	print(DimPlot(data, reduction='umap', label=TRUE, group.by='seurat_clusters')+ggtitle('UMAP: cells grouped by regulon activity colored by seurat clusters'))
	print(DimPlot(data, reduction='umap', label=TRUE, group.by='sample')+ggtitle('UMAP: cells grouped by regulon activity colored by samples'))
	print(DimPlot(data, reduction='umap', label=TRUE, group.by='group')+ggtitle('UMAP: cells grouped by regulon activity colored by custom groups'))
	
	print('make_feature_plot: clustering based on AUC done!')
	
	features <- rownames(data)
	
	# print regulon feature plots from AUC assay
	for (i in seq(1, length(features), by=4)) {
		if ((i+3)>length(features)) {
			print(FeaturePlot(object=data, features=features[i:length(features)], reduction='pca', label=TRUE, ncol=2))
			print(FeaturePlot(object=data, features=features[i:length(features)], reduction='umap', label=TRUE, ncol=2))
			break
		}
		print(FeaturePlot(object=data, features=features[i:(i+3)], reduction='pca', label=TRUE, ncol=2))
		print(FeaturePlot(object=data, features=features[i:(i+3)], reduction='umap', label=TRUE, ncol=2))
	}
	
	print('make_feature_plot: regulons plotted from AUC clusters!')
	print('make_feature_plot: done!')
}




make_violin_plot_from_RNA <- function (data, regulons) {
	
	DefaultAssay(data) <- 'RNA'
	
	features <- gsub(pattern='[(+)]', replacement='', x=names(regulons))
	
	# print regulon feature plots from RNA assay
	for (i in seq(1, length(features), by=4)) {
		if ((i+3)>length(features)) {
			print(VlnPlot(object=data, features=features[i:length(features)], group.by='seurat_clusters', ncol=2))
			break
		}
		print(VlnPlot(object=data, features=features[i:(i+3)], group.by='seurat_clusters', ncol=2))
	}
	
	print('make_violin_plot: regulons plotted from RNA clusters!')
}




make_violin_plot_from_AUC <- function (data) {
	
	DefaultAssay(data) <- 'AUC'
		
	features <- rownames(data)
	
	# print regulon feature plots from AUC assay
	for (i in seq(1, length(features), by=4)) {
		if ((i+3)>length(features)) {
			print(VlnPlot(object=data, features=features[i:length(features)], group.by='seurat_clusters', ncol=2))
			break
		}
		print(VlnPlot(object=data, features=features[i:(i+3)], group.by='seurat_clusters', ncol=2))
	}
	
	print('make_violin_plot: regulons plotted from AUC clusters!')
	print('make_violin_plot: done!')
}




viz_loom <- function (obj) {
	
	# load loom file
	# must use relative path
	loom  <- open_loom(paste0('./out_loom/', obj, '_scenic.loom'), mode='r')
	regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
	regulons <- regulonsToGeneLists(regulons_incidMat)

	AUCmat <- get_AUC(loom)
	
	# load seurat object
	data <- readRDS(paste0('../output/', obj, '.rds'))
	data[['AUC']] <- CreateAssayObject(data=AUCmat)
	DefaultAssay(data) <- 'AUC'

	# check if cell names are identical between regulon mat and seurat obj
	if (!identical(colnames(AUCmat), colnames(data))) {
		print('viz_loom: cell names in AUCmat and data not identical!')
		colnames(AUCmat) <- colnames(GetAssayData(data))
	}
	
	# draw feature plots of regulons based on RNA assay
	pdf(file=paste0(obj, '_regulon_feature_plot_from_RNA_assay.pdf'), width=10, height=10)
	make_feature_plot_from_RNA(regulons, data)
	dev.off()
	print('viz_loom: feature plot RNA done!')
	
	# draw feature plots of regulons based on AUC assay
	pdf(file=paste0(obj, '_regulon_feature_plot_from_AUC_assay.pdf'), width=10, height=10)
	make_feature_plot_from_AUC(data)
	dev.off()
	print('viz_loom: feature plot AUC done!')


	pdf(file=paste0(obj, '_regulon_violin_plot_from_RNA_assay.pdf'), width=10, height=10)
	make_violin_plot_from_RNA(data, regulons)
	dev.off()
	
	pdf(file=paste0(obj, '_regulon_violin_plot_from_AUC_assay.pdf'), width=10, height=10)
	make_violin_plot_from_AUC(data)
	dev.off()

}

	# draw heatmaps
	pdf(file=paste0(obj, '_heatmaps.pdf'), width=40, height=50)
	make_heatmap(AUCmat, data)
	dev.off()
	
	# draw custom heatmaps
	pdf(file=paste0(obj, '_custom_heatmaps.pdf'), width=10, height=50)
	make_custom_heatmap(AUCmat, data)
	dev.off()
	
	


write_sheets <- function (obj) {
	
	# load loom file
	# must use relative path
	loom  <- open_loom(paste0('./out_loom/', obj, '_scenic.loom'), mode='r')
	
	# pull info from loom
	regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
	regulons <- regulonsToGeneLists(regulons_incidMat)
	regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')

	# write number/list of regulated genes by each regulon
	regulon_gene <- data.frame(regulon=names(regulons),
			number=lengths(regulons),
			genes=as.character(lapply(1:length(regulons), FUN = function(x) paste(regulons[[x]], collapse=';'))),
			row.names=NULL)
			
	write.csv(regulon_gene, paste0(obj, '_TF_Regulon_Genelist.csv'))

	# write regulon raw AUC
	AUCmat <- AUCell::getAUC(regulonsAUC)
	rownames(AUCmat) <- gsub('[(+)]', '', rownames(AUCmat))
	write.csv(AUCmat, file=paste0(obj, '_TF_Regulon_AUC_raw.csv'))

	# write regulon scaled AUC
	AUCmat <- t(scale(t(getAUC(regulonsAUC[,]))))
	rownames(AUCmat) <- gsub('[(+)]', '', rownames(AUCmat))
	write.csv(AUCmat, file=paste0(obj, '_TF_Regulon_AUC_scaled.csv'))
	
	print('write_sheets: done!')
}




for (filename in list.files('../output/', pattern='.rds', full.names=F)) {
	obj <- gsub('.rds', '', filename)
	print(paste0('Working on ', obj, '...'))
	
	viz_loom(obj)
	#write_sheets(obj)
}

# output several pdf instead of one
# output also violin plots


#obj <- 'sample_123456_group_123'
#loom <- open_loom('test_fx1/fx1_scenic.loom')


# discarded useful lines:
# this scales AUCmat:
# AUCmat <- t(scale(t(getAUC(regulonsAUC[,]))))
# this gets rid of only the '(+)' in the AUCat rownames: 
# rownames(AUCmat) <- gsub('[(+)]', '', rownames(AUCmat))
# however, cannot use gsub with this regex: '[(+)]'
# will grab each of '(', '+', and ')' and replace three times

DefaultAssay(data) <- 'AUC'
names <- rownames(data)[1:5]
names <- gsub('[(+)]', '', names(regulons)[1:5])
FeaturePlot(data, features=names, label=TRUE)
VlnPlot(data, features=names)



data <- data %>% 
	NormalizeData() %>% 
	FindVariableFeatures() %>% 
	ScaleData() %>% 
	RunPCA() %>% 
	RunUMAP(dims=1:25)

data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data)
data <- RunUMAP(data, dims=1:50)
data <- FindNeighbors(data, dims=1:50)
data <- FindClusters(data, resolution=0.5)


print(DimPlot(data, reduction='pca', label=TRUE, group.by='sample'))
print(DimPlot(data, reduction='umap', label=TRUE, group.by='sample'))
print(DimPlot(data, reduction='pca', label=TRUE, group.by='group'))
print(FeatureScatter(object=data, feature1=names[1], feature2=names[2]))
# scatter plot of two features (typically feature expression), across a set of single cells. Cells are colored by their identity class. Pearson correlation between the two features is displayed above the plot
FeaturePlot(data, feature='NR1D1-9gene(s)', reduction='umap', label=TRUE)
FeaturePlot(data, feature=names, reduction='pca', label=TRUE)
FeaturePlot(data, feature=names, reduction='umap', label=TRUE)

