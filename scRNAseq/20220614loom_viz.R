
# 20220614
# gene regulatory network inference and clustering using pyscenic .loom output visualization

library(dplyr)
library(Seurat)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(pheatmap)


setwd('Labs/xinhua/scSeq/regulon/')


# retrieve AUC scores per cell
loom  <- open_loom('Labs/xinhua/scSeq/regulon/fx1_scenic.loom')

regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')

#rm(loom)

regulons[3]

# binary AUC scores (skipped)


# load seurat obj unfiltered merged
data <- readRDS('Labs/xinhua/scSeq/rds/merged.rds')

print(data@meta.data[1,])
# rm col not needed
data@meta.data <- data@meta.data[,-6]
print(levels(factor(data@meta.data$cluster)))
# all clusters, 1-8
length(data@meta.data[,1])
# 38164

data <- subset(data, subset=percent_mt<15)
data <- subset(data, subset=nFeature_RNA>200)

length(data@meta.data[,1])
# 36766

# only look at certain samples/clusters
data <- subset(data, subset=(sample=='fx1'))

length(data@meta.data[,1])
# 6248

# cluster cells
DefaultAssay(data) <- 'RNA'

data <- NormalizeData(data, normalization.method='LogNormalize', scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method='vst', nfeatures=2000)
data <- ScaleData(data, features=rownames(data))
# scale data without mt
data <- ScaleData(data, vars.to.regress='percent_mt')
data <- RunPCA(data, features=VariableFeatures(object=data), verbose=FALSE)
data <- FindNeighbors(data, dims=1:20)
data <- FindClusters(data, resolution = 0.5)
data <- RunUMAP(data, dims = 1:20)
data <- RunTSNE(data, dims = 1:20) 

saveRDS(data, 'rds/fx1_for_scenic_viz.rds')

pdf(file='~/Labs/xinhua/20220614umap.pdf', width=10, height=10)
DimPlot(data, reduction='umap', label=TRUE)
DimPlot(data, reduction='tsne', label=TRUE)
VlnPlot(data, features=c('PAX7', 'MYF5', 'MYOD1'), group.by='seurat_clusters')
dev.off()

Idents(data)[1:10]
colnames(data@meta.data)

# change default output from Idents(), which is originally sample
# which then becomes seurat_clusters
Idents(data) <- factor(data$cluster, levels=c('clust1','clust2','clust3','clust4','clust5','clust6','clust7','unassigned'))

# AUC scores
AUCmat <- AUCell::getAUC(regulonAUC)
# format AUCmat rownames
rownames(AUCmat) <- gsub('[(+)]', '', rownames(AUCmat))

dim(AUCmat)
# 348 6248
dim(GetAssayData(data))
# 36601 6248

data[['AUC']] <- CreateAssayObject(data=AUCmat)
# Error: All cells in the object being added must match the cells in this object

length(colnames(AUCmat))==length(colnames(data))
# TRUE
identical(colnames(AUCmat), colnames(data))
# FALSE
# cell names do not match

AUCmat[1:4, 1:4]
# regulons as rownames, cells as colnames
GetAssayData(data)[1:4, 1:4]
# genes as rownames, cells as colnames
# however, all cell names have prefix: fx1_

co <- colnames(GetAssayData(data))
co[1:4]
identical(sub('fx1_', '', co), colnames(AUCmat))
# TRUE
# subbing 'fx1_' with '' has rm the unwanted prefix

RenameCells(data[['RNA']], new.names=sub('fx1_', '', colnames(GetAssayData(data))))

identical(colnames(AUCmat), colnames(data))
# FALSE
# seems like RenameCells() had no use
str(data)
str(data$RNA@data)
str(data$RNA@data@Dimnames)
data$RNA@data@Dimnames[[2]] <- sub('fx1_', '', colnames(GetAssayData(data)))
# Error: All cells in the object being added must match the cells in this object

# okay, I think it is too difficult to modify cellnames of a seurat obj
# next time, when creating merged.rds, do not add the sample prefixes

# I'll make the AUCmat colnames identical with the AssayObject colnames
colnames(AUCmat) <- colnames(GetAssayData(data))
data[['AUC']] <- CreateAssayObject(data=AUCmat)

DefaultAssay(data) <- 'AUC'

data <- ScaleData(data, assay='AUC', features=rownames(AUCmat))

# plot the first four regulons
FeaturePlot(data, features=rownames(AUCmat)[1:4])

# heatmap
pheatmap(AUCmat, show_colnames=F, annotation_col=data@meta.data$cluster)
# Error in annotation_col[colnames(mat), , drop = F] : 
#  incorrect number of dimensions

cellInfo <- get_cell_annotation(loom)
colnames(cellInfo)
cellInfo[1:4,]
cellInfo$cluster <- data@meta.data$cluster
cellInfo[1:4,]
cellType <- subset(cellInfo, select='cluster')
pheatmap(AUCmat, show_colnames=F, annotation_col=cellType)
# Error in check.length("fill") : 
#  'gpar' element 'fill' must not be length 0
# this error arose due to not have EXACT SAME rownames between :test: and :annotation_col:, or pheatmap cannot correlate values

identical(rownames(cellType), colnames(AUCmat))
# FALSE
rownames(cellType)[1:4]
colnames(AUCmat)[1:4]
rownames(cellType) <- colnames(GetAssayData(data))

length(data@meta.data$seurat_clusters)
cellInfo$seurat_cluster <- data@meta.data$seurat_clusters
colnames(cellInfo)
cellSeuratCluster <- subset(cellInfo, select='seurat_cluster')
rownames(cellSeuratCluster) <- colnames(GetAssayData(data))


#
#
#
# also add mito%

pheatmap(AUCmat[1:10,], show_colnames=F, annotation_col=cellSeuratCluster)
pheatmap(AUCmat[1:10,], show_colnames=F, annotation_col=cellType)
# this works
pheatmap(AUCmat[1:10,], show_colnames=F, annotation_col=cellInfo$cluster)
# this cannot work

dim(cellType)# 6248 1
dim(cellInfo$cluster)# NULL
dim(cellInfo[,1])# NULL
dim(as.matrix(cellInfo[,1]))# 6248 1, but still error: no 'dimnames' attribute for array

# how about cbind()?
pheatmap(AUCmat[1:10,], show_colnames=F, annotation_col=cbind(cellSeuratCluster, cellType))
# worked! 

# how about the entire cellInfo matrix
rownames(cellInfo) <-  colnames(GetAssayData(data))
pheatmap(AUCmat[1:10,], show_colnames=F, annotation_col=cellInfo)
# worked!

pheatmap(AUCmat, show_colnames=F, annotation_col=cellInfo, filename='~/Labs/xinhua/pheatmap.png', width=20, height=80)


pheatmap(AUCmat, show_colnames=F, annotation_col=cellInfo$seurat_cluster)





pheatmap(AUCmat[rownames(AUCmat) %in% c('BCLAF1','BCL6','ATF6B','ATF5','ABL1'),], show_colnames=F)

#######

following: https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/importing_pySCENIC.html

# import regulons, AUC, embeddings from loom file
# regulons_incidMat <- get_regulons(loom)
# regulons <- regulonsToGeneLists(regulons_incidMat)
# regulonAUC <- get_regulonsAuc(loom)
regulonAucThresholds <- get_regulonThresholds(loom)
embeddings <- get_embeddings(loom)

# import saved info about expression matrix
# cellInfo <- get_cell_annotation(loom)
exprMat <- get_dgem(loom)
clusterings <- get_clusterings_withName(loom)

close_loom(loom)


# motif enrichment analysis
# read table ###### look in the terminal, see if you find any tsv
motifsDf <- data.table::fread('/Users/shell/Labs/xinhua/scSeq/regulon/resource/motifs-v9-nr.hgnc-m0.001-o0.0.tbl', header=T, sep='\t')
# visualize: not working
tableSubset <- motifsDf[motif_name=='Abd-B']##
# only view a few lines
#tableSubset <- tableSubset[1:20,]
colsToShow <- colnames(motifsDf)[-c(2, 9:11)]#???
viewMotifs(tableSubset, colsToShow=colsToShow)
# Error
# docpage: https://rdrr.io/github/aertslab/SCENIC/man/viewMotifs.html


