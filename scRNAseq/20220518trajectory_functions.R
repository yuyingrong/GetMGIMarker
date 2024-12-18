
# 202205013
# Monocle 2 - "only the essentials"
# the trajectory script organized into functions for running on multiple samples on the server

library(Seurat)
library(dplyr)
library(ggplot2)
library(monocle)

setwd('~/fx/')

create_cds <- function(data) {
	# input: data of wanted samples/clusters
	
	# make a matrix from non-zero expression genes
	exp <- data@assays$RNA@counts[which(rowSums(data@assays$RNA@counts)!=0),]
	
	# create CellDataSet object from gene expression matrix
	cds <- newCellDataSet(as(as.matrix(exp),'sparseMatrix'),
					phenoData=new('AnnotatedDataFrame', data=data@meta.data),
					featureData=new('AnnotatedDataFrame', data.frame(gene_short_name=row.names(exp), row.names=row.names(exp))),
					lowerDetectionLimit=0.5,
					expressionFamily=negbinomial.size()
					)
	
	# rate determining steps
	cds <- estimateSizeFactors(cds)# necessary for calling dispersions
	cds <- estimateDispersions(cds)
	
	# now cds is ready for further calculation and plotting
	print('estimateDispersions done!')
	return(cds)
	
}

order_by_meanExpression <- function(cds) {
	# input: cds after size/dispersion estimation
	
	# find genes by mean_expression only
	disp_genes <- subset(dispersionTable(cds), subset=(mean_expression>=0.1))$gene_id
	cds <- setOrderingFilter(cds, ordering_genes=disp_genes)
	
	# rate determining step
	cds <- reduceDimension(cds, max_components=2, method='DDRTree')
	
	cds <- orderCells(cds)
	
	# now cds is ready for visualization
	print('ordering cells done!')
	return(cds)
}

order_by_SeuratVariableFeatures <- function(data, cds) {
	# input: cds after size/dispersion estimation
	
	# find Seurat variable genes
	data <- FindVariableFeatures(data)
	var_features <- VariableFeatures(data)
	
	cds <- setOrderingFilter(cds, ordering_genes=var_features)

	# rate determining step
	cds <- reduceDimension(cds, max_components=2, method='DDRTree')
	
	cds <- orderCells(cds)
	
	# now cds is ready for visualization
	return(cds)
}

viz <- function(cds) {
	# plot ordering genes
	print(plot_ordering_genes(cds))
	
	# color by pseudotime
	print(plot_cell_trajectory(cds, color_by='Pseudotime', cell_size=0.3, show_tree=TRUE))
	print(plot_cell_trajectory(cds, color_by='Pseudotime', cell_size=0.3, show_tree=TRUE)+ facet_wrap('~State', ncol=3))

	# color by state
	print(plot_cell_trajectory(cds, color_by='State', cell_size=0.3, show_tree=TRUE))
	print(plot_cell_trajectory(cds, color_by='State', cell_size=0.3, show_tree=TRUE) + facet_wrap('~State', ncol=3))

	# color by samples
	print(plot_cell_trajectory(cds, color_by='sample', cell_size=0.3, show_tree=TRUE))

	print(plot_cell_trajectory(cds, color_by='sample', cell_size=0.3, show_tree=TRUE) + facet_wrap('~sample', ncol=3))

	# color by clusters
	print(plot_cell_trajectory(cds, color_by='ident', cell_size=0.3, show_tree=TRUE))

	print(plot_cell_trajectory(cds, color_by='ident', cell_size=0.3, show_tree=TRUE) + facet_wrap('~ident', ncol=3))

	print(plot_cell_trajectory(cds, color_by='cluster', cell_size=0.3, show_tree=TRUE))

	print(plot_cell_trajectory(cds, color_by='cluster', cell_size=0.3, show_tree=TRUE) + facet_wrap('~cluster', nrow=2))

}

trajectory <- function(samples, clusters) {
	
	# unload/reload dataset
	data <- NULL
	data <- readRDS('rds/data.rds')
	
	# subset by wanted samples/clusters	
	data <- subset(data, cells=rownames(data@meta.data)[which(data$sample %in% paste0('fx', samples))])
	data <- subset(data, cells=rownames(data@meta.data)[which(data$cluster %in% paste0('clust', clusters))])
	
	print(levels(factor(data@meta.data$sample)))
	print(levels(factor(data@meta.data$cluster)))
	
	# create cds from data
	cds <- create_cds(data)

	# order cells by dispersionTable: mean expression only
	cds <- order_by_meanExpression(cds)
	saveRDS(object=cds, file=paste0('~/fx/rds/sample_', paste(samples, collapse=''), '_cluster_', paste(clusters, collapse=''), '_ordered_cells_by_expr.rds'))
	#cds <- order_by_SeuratVariableFeatures(data, cds)
	#saveRDS(object=cds, file=paste0('~/fx/rds/sample_', paste(samples, collapse=''), '_cluster_', paste(clusters, collapse=''), '_ordered_cells_by_SeuratVariableFeatures.rds'))

	
	# plot trajectories
	pdf(file=paste0('~/fx/output/sample_', paste(samples, collapse=''), '_cluster_', paste(clusters, collapse=''), '_trajectory.pdf'), width=20, height=20)
	viz(cds)
	dev.off()

}


trajectory(c(1,3,4,5,6), 1:6)
#trajectory(c(1,3,4,5,6), c(1,3,4,6))
#trajectory(1:6, 1:6)
#trajectory(1:6, 1:3)
