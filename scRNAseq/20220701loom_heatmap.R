
# 20220701
# pyscenic loom output viz
# AUCell binarize


library(dplyr)
library(Seurat)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(pheatmap)


setwd('Labs/xinhua/scSeq/regulon/')


# retrieve AUC scores per cell
loom  <- open_loom('fx1_scenic.loom')

# pull info from loom
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')

#AUCmat <- AUCell::getAUC(regulonsAUC)
#AucThresholdsMat <- AUCell::getAUC(regulonsAucThresholds)
# format AUCmat rownames
#rownames(AUCmat) <- gsub('[(+)]', '', rownames(AUCmat))
#pheatmap(AUCmat, show_colnames=F, annotation_col=cellInfo$seurat_cluster)

data <- readRDS('../rds/fx1_for_scenic_viz.rds')
data


# Basic heatmap

AUCmat <- t(scale(t(getAUC(regulonsAUC[,]))))
dim(AUCmat)
# 348 6248
length(AUCmat)
# 2174304
length(AUCmat[AUCmat>=2])
# 81037
length(AUCmat[AUCmat <=(-2)])
# 9886
AUCmat[AUCmat>=2] <- 2
AUCmat[AUCmat <=(-2)] <- (-2)
AUCmat[1:4, 1:4]
ac <- data.frame(group=as.character(Idents(data)))

#AUCmat <- AUCmat[,colnames(AUCmat) %in% colnames(data)]
rownames(ac) <- colnames(AUCmat)

# if you need to use a subset of tf
#cg <- read.table('my_chosen_tf.txt')[,1]
#cg
#cg <- AUCmat[rownames(AUCmat) %in% cg,]

#pheatmap(cg, show_colnames=F, show_rownames=T, annotation_col=ac)
pheatmap(AUCmat, show_colnames=F, show_rownames=T, annotation_col=ac)
# Error: 'gpar' element 'fill' must not be length 0
# this error arises when you do not have EXACT SAME rownames between test and annotation_col, so pheatmap cannot correlate values

identical(rownames(ac), colnames(AUCmat))
length(rownames(ac))
# 6248
length(colnames(AUCmat))
# 6248
head(rownames(ac))
# "1" "2" "3" "4" "5" "6"
head(colnames(AUCmat))
# "AAACCCAAGATGGTAT-1" "AAACCCAAGTGGTGAC-1" "AAACCCACACAAGGTG-1"

rownames(ac) <- colnames(AUCmat)
pheatmap(AUCmat, show_colnames=F, show_rownames=T, annotation_col=ac)
# just random first 10 genes
pheatmap(AUCmat[1:10,], show_colnames=F, show_rownames=T, annotation_col=ac)


# Basic heatmap with more info

# swap clust info as Idents(), rather than seurat_cluster
Idents(data) <- factor(data$cluster, levels=c('clust1','clust2','clust3','clust4','clust5','clust6','clust7','unassigned'))

# AUC scores
AUCmat <- AUCell::getAUC(regulonsAUC)
AUCmat <- t(scale(t(getAUC(regulonsAUC[,]))))
# format AUCmat rownames
rownames(AUCmat) <- gsub('[(+)]', '', rownames(AUCmat))

length(colnames(AUCmat))==length(colnames(data))
identical(colnames(AUCmat), colnames(data))
# because right now colnames(AUCmat) has no 'fx1_' prefix, while colnames(data) does
# colnames don't match, but subsequent steps require them to match
# make colnames(AUCmat) identical with colnames(data)
colnames(AUCmat) <- colnames(GetAssayData(data))


# plot the first four regulons
FeaturePlot(data, features=rownames(AUCmat)[1:4])

# create annotation_col for plotting pheatmap
cellInfo <- get_cell_annotation(loom)
# now nGene and nUMI are added
cellInfo$cluster <- data@meta.data$cluster
# now self-defined cluster info is added
cellInfo$seurat_cluster <- data@meta.data$seurat_cluster
cellInfo$percent_mt <- data@meta.data$percent_mt
cellInfo[1:4,]

# unify colnames, and put the cellInfo in
rownames(cellInfo) <-  colnames(GetAssayData(data))
pheatmap(AUCmat[1:10,], show_colnames=F, annotation_col=cellInfo)

# oh no, forget to format the data
summary(AUCmat)
AUCmat[AUCmat>=2] <- 2
AUCmat[AUCmat <=(-2)] <- (-2)
AUCmat[1:4, 1:4]

pheatmap(AUCmat[1:10,], show_colnames=F, annotation_col=cellInfo)

# black and white heatmap; but not sure if this is meaningful
pheatmap(AUCmat[1:10,], show_colnames=F, annotation_col=cellInfo, color=c("black", "white"))





# produce sheets

# use lengths() to find length of every element in a list
regulon_gene <- data.frame(regulon=gsub('[(+)]', '', names(regulons)),
			number=lengths(regulons),
			genes=as.character(lapply(1:length(regulons), FUN = function(x) regulons[[x]])),
			row.names=NULL)

write.csv(regulon_gene, 'TF_Regulon_Genelist.csv')
# not perfect, but works

# write.csv() writes the entire file content; use write.table() instead
for (i in 1:length(regulons)) {
	X <- data.frame(Regulon=gsub('[(+)]', '', names(regulons)[i]), 
				Gene_number=length(regulons[[i]]), 
				Genes=paste(regulons[[i]], collapse=';')
				)
	write.table(X, file='TF_Regulon_Genelist.csv', append=T, quote=F, sep=',', eol='\n', row.names=F, col.names=T)
}
# prints headers per line; undesirable

# header
paste(c('Regulon', 'Gene_number', 'Genes'), collapse=',') %>%
	write.table(file='TF_Regulon_Genelist.csv', append=T, quote=F, sep=',', eol='\n', row.names=F, col.names=F)
# content
for (i in 1:length(regulons)) {
	c(gsub('[(+)]', '', names(regulons)[i]), 
		length(regulons[[i]]), 
		paste(regulons[[i]], collapse=';')
		) %>%
	paste(collapse=',') %>%
	write.table(file='TF_Regulon_Genelist.csv', append=T, quote=F, sep=',', eol='\n', row.names=F, col.names=F)
}
# good



AUCmat <- AUCell::getAUC(regulonsAUC)
rownames(AUCmat) <- gsub('[(+)]', '', rownames(AUCmat))
write.csv(AUCmat, file='TF_Regulon_AUC_raw.csv')


AUCmat <- t(scale(t(getAUC(regulonsAUC[,]))))
rownames(AUCmat) <- gsub('[(+)]', '', rownames(AUCmat))
write.csv(AUCmat, file='TF_Regulon_AUC_scale.csv')





# custom heatmaps

heatmap(AUCmat, Rowv=NA, Colv=NA)

# unify colnames, or cell names
colnames(AUCmat) <- colnames(GetAssayData(data))

data[['AUC']] <- CreateAssayObject(data=AUCmat)
DefaultAssay(data) <- 'AUC'

#data <- ScaleData(data, assay='AUC', features=rownames(AUCmat))

rowSums(subset(data, seurat_clusters==0))

nrow(data[['AUC']])
# 348


AUC_by_seurat_clusters <- data.frame(
	cluster0=rowMeans(subset(data, seurat_clusters==0)),
	cluster1=rowMeans(subset(data, seurat_clusters==1)),
	cluster2=rowMeans(subset(data, seurat_clusters==2)),
	cluster3=rowMeans(subset(data, seurat_clusters==3)),
	cluster4=rowMeans(subset(data, seurat_clusters==4)),
	cluster5=rowMeans(subset(data, seurat_clusters==5)),
	cluster6=rowMeans(subset(data, seurat_clusters==6)),
	cluster7=rowMeans(subset(data, seurat_clusters==7)),
	row.names=rownames(data[['AUC']])
)

?pheatmap

pdf(file='20220706pheatmap_avg_AUC_by_seurat_clusters.pdf', width=6, height=42)
pheatmap(mat=AUC_by_seurat_clusters, cluster_rows=T, cluster_cols=F, show_colnames=T, annotation_col=NA)

dev.off()


levels(factor(data$cluster))
AUC_by_custom_groups <- data.frame(
	group1=rowMeans(subset(data, cluster=='clust1')),
	group2=rowMeans(subset(data, cluster=='clust2')),
	group3=rowMeans(subset(data, cluster=='clust3')),
	group4=rowMeans(subset(data, cluster=='clust4')),
	group5=rowMeans(subset(data, cluster=='clust5')),
	group6=rowMeans(subset(data, cluster=='clust6')),
	group7=rowMeans(subset(data, cluster=='clust7')),
#	unassigned=rowMeans(subset(data, cluster=='unassigned')),
	row.names=rownames(data[['AUC']])
)

pdf(file='20220706pheatmap_avg_AUC_by_custom_groups.pdf', width=6, height=42)
pheatmap(mat= AUC_by_custom_groups, cluster_rows=T, cluster_cols=F, show_colnames=T, annotation_col=NA)

dev.off()


# thinking of a way to type less
calc_row <- function(x) {
	assign(paste0('cluster', x), rowSums(subset(data, seurat_clusters==x)))
}
length(levels(data$seurat_clusters))
lapply(0:length(levels(data$seurat_clusters)), calc_row)
# not working 

# another attempt
mean_AUC <- function(type) {
	categories <- levels(factor(data[[type]]))
	out <- data.frame()
	for (i in 1:length(categories)) {
		out$i <- rowMeans(subset(data, type==categories[i]))
	}
}

# succeed!
AUC_by_seurat_clusters <- data.frame(
	cluster0=rowMeans(subset(data, seurat_clusters==0)), 
	row.names=rownames(data[['AUC']]))


for (i in levels(data$seurat_clusters)[-1]) {
	print(i)
	AUC_by_seurat_clusters[[paste0('cluster',i)]] <- rowMeans(subset(data, seurat_clusters==i))
	
}

colnames(AUC_by_seurat_clusters)
dim(AUC_by_seurat_clusters)
AUC_by_seurat_clusters[1:4,1:3]
