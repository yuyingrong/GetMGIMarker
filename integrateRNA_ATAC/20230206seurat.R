
# same-cell single-cell RNA- and ATAC-seq multiomic data integration
# merging two datasets in one seurat object
# integrating same-cells and clustering

# install.packages('Seurat')
library(Seurat)
library(ggplot2)
library(dplyr)
#library(stringr)
#install.packages('scater')
#library(scater) not available in this R version


setwd('~/Labs/wijst/')

# script copied from "20220428clustering.R"

viz_QC <- function(data) {
	
	# cell count
	print(
		data@meta.data %>% ggplot(aes(x=orig.ident, fill=orig.ident))+
		geom_bar()+
		theme_classic()+
		theme(axis.text.x=element_text(angle=45,hjust=1))+
		theme(plot.title=element_text(hjust=0.5,face='bold'))+
		ggtitle('Cell count by sample')
	)
	
	# overall
	print(
		data@meta.data %>% ggplot(aes(x=nCount_RNA,y=nFeature_RNA,color=percent_mt))+
		geom_point(alpha=0.5)+
		scale_color_gradient(low='gray',high='black')+
		stat_smooth(method=lm)+
		scale_x_log10()+
		scale_y_log10()+
		theme_classic()+
		geom_vline(xintercept=500)+
		geom_hline(yintercept=250)+
		facet_wrap(~orig.ident)
	)
	
	# mt% and gene count	
	print(VlnPlot(data, features='percent_mt', group.by='orig.ident') + ggtitle('Mito% by sample'))
	print(VlnPlot(data, features='nFeature_RNA', group.by='orig.ident') + ggtitle('Gene count by sample'))
		
	# distribution of mt% and gene count

	print(summary(data@meta.data$percent_mt))
	#mitoDrop <- isOutlier(data@meta.data$percent_mt, nmads=3, type='higher')
	print(
		ggplot(data@meta.data, aes(x=percent_mt))+
		geom_histogram(binwidth=0.5)+
		#geom_vline(xintercept=as.numeric(attributes(mitoDrop)$thresholds[2]), linetype='dashed', color='red')+
		#geom_text(aes(x=as.numeric(attributes(mitoDrop)$thresholds[2], digits=3), y=-20, label=as.character(round(attributes(mitoDrop)$thresholds[2], digits=3))))+
		ggtitle('Mito% of merged dataset')
	)
		
	print(summary(data@meta.data$nFeature_RNA))
	#geneDrop <- isOutlier(data@meta.data$nFeature_RNA,nmads=3,type='both')
	print(
		ggplot(data@meta.data,aes(x=nFeature_RNA))+
		geom_histogram(binwidth=50)+
		#geom_vline(xintercept=as.numeric(attributes(geneDrop)$thresholds[1]),linetype='dashed',color='red')+
		#geom_vline(xintercept=as.numeric(attributes(geneDrop)$thresholds[2]),linetype='dashed',color='red')+
		#geom_text(aes(x=attributes(geneDrop)$thresholds[[1]]+2,y=-10,label=as.character(attributes(geneDrop)$thresholds[[1]])))+
		#geom_text(aes(x=attributes(geneDrop)$thresholds[[2]]+2,y=-10,label=as.character(attributes(geneDrop)$thresholds[[2]])))+
		ggtitle('Gene count of merged dataset')
	)

}

# read in 10X data
data <- Read10X(data.dir = "raw/filtered_feature_bc_matrix")

# reported containing two types of data

str(data)
# $ Gene Expression and $ Peaks

# retrieve only Gene Expression from data
pbmc <- CreateSeuratObject(counts = data$`Gene Expression`, project = "pbmc", min.cells = 1, min.features = 0)

pbmc$percent_mt <- PercentageFeatureSet(pbmc, pattern = '^MT-')

pdf(file='/out/filtered_qc.pdf', width=10, height=10)

viz_QC(pbmc)

#data <- subset(data, subset=(nFeature_RNA>500 & percent_mt<8))

dev.off()

###
# atac data
pbmc.atac <- CreateSeuratObject(counts = data$Peaks, project = "pbmc", min.cells = 1, min.features = 0)

# pbmc[['ATAC']] <- pbmc.atac
# not possible

# 20230228
setwd('~/Labs/wijst/')
library(Seurat)
library(ggplot2)
library(dplyr)

data <- Read10X(data.dir = "raw/filtered_feature_bc_matrix_gz/")
# path must contain use .tsv.gz files

data <- CreateSeuratObject(counts = data$`Gene Expression`, project = "pbmc", min.cells = 1, min.features = 0)
data
# 29717 genes across 11898 cells

data$percent_mt <- PercentageFeatureSet(data, pattern = '^MT-')

#
# QC parameters based on mito%

summary(data@meta.data$percent_mt)
# median 9.747, mean 10.265
mad(data@meta.data$percent_mt)
# mad 3.414726
mitoDrop <- median(data@meta.data$percent_mt)+2*mad(data@meta.data$percent_mt)
mitoDrop
# 16.5768

ggplot(data@meta.data, aes(x=percent_mt))+
  geom_histogram(binwidth=0.5)+
  geom_vline(xintercept=mitoDrop, linetype='dashed', color='blue')+
  geom_text(aes(x=as.numeric(mitoDrop), y=-20, label=as.character(round(mitoDrop, digits=2))))+
  ggtitle('Mito%')

# QC parameters based on number of genes
# not using
summary(data@meta.data$nFeature_RNA)
geneDrop <- c(median(data@meta.data$nFeature_RNA)-2*mad(data@meta.data$nFeature_RNA),
              median(data@meta.data$nFeature_RNA)+2*mad(data@meta.data$nFeature_RNA))
geneDrop
# 735.3064 2917.6936

ggplot(data@meta.data,aes(x=nFeature_RNA))+
  geom_histogram(binwidth=50)+
  geom_vline(xintercept=geneDrop[1],linetype='dashed',color='blue')+
  geom_vline(xintercept=geneDrop[2],linetype='dashed',color='blue')+
  geom_text(aes(x=geneDrop[1],y=-10,label=as.character(round(geneDrop[1], digits=2))))+
  geom_text(aes(x=geneDrop[2],y=-10,label=as.character(round(geneDrop[2], digits=2))))+
  ggtitle('Gene count')

# QC parameters based on number of mRNA molecules

summary(data@meta.data$nCount_RNA)
umiDrop <- c(median(data@meta.data$nCount_RNA)-2*mad(data@meta.data$nCount_RNA),
              median(data@meta.data$nCount_RNA)+2*mad(data@meta.data$nCount_RNA))
umiDrop
# 705.5528 6849.4472

ggplot(data@meta.data,aes(x=nCount_RNA))+
  geom_histogram(binwidth=50)+
  geom_vline(xintercept=umiDrop[1],linetype='dashed',color='blue')+
  geom_vline(xintercept=umiDrop[2],linetype='dashed',color='blue')+
  geom_text(aes(x=umiDrop[1],y=-10,label=as.character(round(umiDrop[1], digits=2))))+
  geom_text(aes(x=umiDrop[2],y=-10,label=as.character(round(umiDrop[2], digits=2))))+
  ggtitle('mRNA count')


#
# filter cells based on QC parameters
# data <- subset(data, subset=(nFeature_RNA>500 & percent_mt<8))
# instead of subset(), I will make a table with cell barcodes and QC parameters

nrow(data@meta.data)
# 11898, number of cells

# first, extract all atac-seq cell barcodes
# since I already have this file from doublet labeling
rna <- read.csv('out/scrublet_predictions.csv', header=FALSE)
nrow(rna)
# 11898
rna[1:4,]

# add colnames
colnames(rna) <- c('barcode', 'doublet')

# replace (True with 1) and (False with 0) for doublet
table(rna$doublet)
# False  True 
# 10994   904
rna$doublet[rna$doublet == 'True'] <- 1
rna$doublet[rna$doublet == 'False'] <- 0
table(rna$doublet)

# add cols
colnames(data@meta.data)
# "orig.ident"   "nCount_RNA"   "nFeature_RNA" "percent_mt"

# import clusters
pbmc <- readRDS('rds/signac_pbmc_multiome_unfiltered.rds')

colnames(pbmc@meta.data)
DefaultAssay(pbmc) <- 'RNA'

# cell clusters and cell types anchored to Azimuth ref clusters
DimPlot(pbmc, reduction = 'ref.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

# what about seurat defined cell clusters?
pbmc
# 4 dimensional reductions calculated: pca, integrated_dr, ref.umap, lsi
# no UMAP, ...doesn't look like clustering has been done!

# clustering cells
library(dplyr)
library(Seurat)

# should I remove low quality cells before clustering?
# create a new obj first

# to remove doublets, add doublet into to obj data
identical(rownames(data@meta.data), rna$barcode)
# TRUE, then I don't need to worry about order
data@meta.data$doublet <- rna$doublet

# remove doublets by scrublet predictions
data1 <- subset(data, doublet==0)
data1
# 10994 cells, used to be 11898 cells; removed 904 doublets

# remove cells by seurat QC metrics
# we don't filter based on gene counts bc cells like RBC can have many UMIs but few genes
data1 <- subset(data1, subset=(
#  nFeature_RNA > geneDrop[1] & 
  nCount_RNA > umiDrop[1] & 
  percent_mt < mitoDrop
))

data1
# if removed by geneDrop: 10288 cells; further removed 706 low quality cells
# if removed by umiDrop: 10338 cells; further removed 656 cells

# run normalize, scale, and find variable features in one command
data1 <- SCTransform(data1, method = "glmGamPoi", vars.to.regress = "percent_mt", verbose = FALSE)

# dimension reduction
data1 <- RunPCA(data1, verbose = FALSE)
DimPlot(data1, reduction = 'pca', label = TRUE)

data1 <- RunUMAP(data1, dims = 1:50, verbose = FALSE)
DimPlot(data1, reduction = 'umap', label = TRUE)

# cluster
data1 <- FindNeighbors(data1, dims = 1:50, verbose = FALSE)
data1 <- FindClusters(data1, verbose = FALSE)
data1
DimPlot(data1, reduction = 'umap', label = TRUE)

# save filtered clustered rds
saveRDS(data1, file = 'rds/seurat_pbmc_rna_filtered_clustered.rds')

colnames(data1@meta.data)

rna$nUMI <- data@meta.data$nCount_RNA[match(rna$barcode, rownames(data@meta.data))]
rna$nGene <- data@meta.data$nFeature_RNA[match(rna$barcode, rownames(data@meta.data))]
rna$mt <- data@meta.data$percent_mt[match(rna$barcode, rownames(data@meta.data))]
rna$seurat_cluster <- data1@meta.data$seurat_clusters[match(rna$barcode, rownames(data1@meta.data))]
rna$predicted.celltype.l1 <- pbmc@meta.data$predicted.celltype.l1[match(rna$barcode, rownames(pbmc@meta.data))]
rna$predicted.celltype.l2 <- pbmc@meta.data$predicted.celltype.l2[match(rna$barcode, rownames(pbmc@meta.data))]

rna$seurat_cluster
# contains <NA>, good
rna$predicted.celltype.l1
# every cell has been assigned a type

rna[1:5,]
# without a seurat_cluster means the cell does not pass QC

write.csv(rna, 'out/rna_metadata.csv', row.names=FALSE)

# NOTE!!!
# the file can be improved by: re-running azimuth with filtered data1

# retrive table
rna <- read.csv('out/rna/rna_metadata.csv', header=TRUE)
rna[1:10,]

# add info to data1 metadata
data1@meta.data[1:4,]
identical(rownames(data1@meta.data), rna$barcode)
# FALSE
# oh, I subset() the obj; rna has 11898 rows, whereas data1 metadata has 10338
# now I cannot directly transfer azimuth labels
data1@meta.data$predicted.celltype.l1 <- rna$predicted.celltype.l1[match(rownames(data1@meta.data), rna$barcode)]
data1@meta.data$predicted.celltype.l2 <- rna$predicted.celltype.l2[match(rownames(data1@meta.data), rna$barcode)]
data1@meta.data[1:4,]

# view azimuth labels on umap
DimPlot(data1, reduction = 'umap', group.by = 'predicted.celltype.l1', label = TRUE)
DimPlot(data1, reduction = 'umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE)

# save updated metadata
saveRDS(data1, file = 'rds/seurat_pbmc_rna_filtered_clustered.rds')


#
# export barcodes, features, and count matrix

#data1 <- readRDS(file = 'rds/seurat_pbmc_rna_filtered_clustered.rds')
# write the files out from unfiltered data, not filtered data1!
# data1 is only for seurat clustering

# found wanted info!
str(data@assays$RNA@counts)
str(data@assays$RNA@data)

identical(data@assays$RNA@counts,data@assays$RNA@data)
# TRUE

str(data@assays$RNA@counts)

# create a list of barcodes
# write.csv(data@assays$RNA@counts@Dimnames[2][[1]], file = 'out/rna/rna_barcodes.tsv', row.names = FALSE, col.names = NULL, quote = FALSE)
# did not use write.csv() because it kept attaching a header
write.table(data@assays$RNA@counts@Dimnames[2][[1]], file = 'out/rna/rna_barcodes.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)
# check if write out are correct
data@assays$RNA@counts@Dimnames[2][[1]][1:4]
# good

# create a list of features
write.table(data@assays$RNA@counts@Dimnames[1][[1]], file = 'out/rna/rna_features.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)
data@assays$RNA@counts@Dimnames[1][[1]][1:4]

# write a matrix market file
install.packages('Matrix')
library(Matrix)
#writeMM(data1@assays$RNA@counts, 'out/rna/rna_matrix.mtx')
writeMM(data@assays$RNA@counts, 'out/rna/rna_matrix.mtx')

# export filtered data
data1 <- readRDS(file = 'rds/seurat_pbmc_rna_filtered_clustered.rds')
str(data1@assays$RNA@counts)
# 29717 10338
# that's the correct number of cells after doublets and mt% and nUMI removed
# sparse matrix
writeMM(data1@assays$RNA@counts, 'out/rna/rna_matrix.mtx')
# barcodes
write.table(data1@assays$RNA@counts@Dimnames[2][[1]], file = 'out/rna/rna_barcodes.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)
# features
write.table(data1@assays$RNA@counts@Dimnames[1][[1]], file = 'out/rna/rna_features.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)
# metadata
str(data1@meta.data)
# nCount_RNA nFeature_RNA percent_mt doublet seurat_clusters predicted.celltype.l1 predicted.celltype.l2
# I don't want everything here
str(data1@meta.data[,-c(1,5,6,7,8)])
# that's better, I don't even need doublet info, since doublets have been removed
# write.csv(data1@meta.data[,-c(1,5,6,7,8)], 'out/rna/rna_metadata.csv', row.names=FALSE)
# file written out this way did not have barcodes!
data1@meta.data[1,]
# barcodes are in rownames
colnames(data1@meta.data)[1] <- 'barcode'
data1@meta.data[,1] <- rownames(data1@meta.data)
# also round the mt% calculation to 2 digits
data1@meta.data$percent_mt <- round(data1@meta.data$percent_mt, 2)
write.csv(data1@meta.data[,-c(5,6,7,8)], 'out/rna/rna_metadata.csv', row.names=FALSE)

# check if SCTransform affects the raw count matrix
# data obj here has bad cells removed, but no normalization performed
identical(data@assays$RNA@counts, data1@assays$RNA@counts)
# these matrices are identical, so my output is the un-normalized count matrix

# how many cells intersect between RNA and ATAC?
rna <- read.csv('out/rna/rna_metadata.csv', header=TRUE)
rna$barcode[1:4]
atac <- read.csv('out/atac/atac_metadata.csv', header=TRUE)
atac$barcode[1:4]

length(rna$barcode)
# 10338, from 11898 (if not rm by qc)
length(atac$barcode)
# 9939
length(intersect(rna$barcode, atac$barcode))
# 9213, from 9902 (if rna does not rm by qc)

# slice two sparse matrices by intersecting barcodes
library(Matrix)
atacMM <- readMM(file='out/atac/atac_matrix.mtx')
rnaMM <- readMM(file='out/rna/rna_matrix.mtx')
str(atacMM)
str(rnaMM)
# so all the col/row names are gone

# import the barcodes
atac <- read.table(file = 'out/atac/atac_barcodes.tsv', header = FALSE)
atac[1:4,]
nrow(atac)
# 9939
rna <- read.table(file = 'out/rna/rna_barcodes.tsv', header = FALSE)
rna[1:4,]
nrow(rna)
# 10338

joint <- intersect(atac, rna)
nrow(joint)
joint[1:4,]
# 9213

# match
match(c(4,1),c(1,2,3,4))
# 4,1
match(joint, atac)
# NA, why?
match(joint[,1], atac[,1])
length(match(joint[,1], atac[,1]))
# 9213, good
# now I have the indices, slice atacMM
dim(atacMM)
# 24919 9939, slice by col
joint_atacMM <- atacMM[,match(joint[,1], atac[,1])]
dim(joint_atacMM)
# 24919 9213
joint_rnaMM <- rnaMM[,match(joint[,1], rna[,1])]
dim(joint_rnaMM)
# 29717 9213

# slice metadata
atac_meta <- read.table(file = 'out/scMVP_input/atac_metadata.csv', header = TRUE, sep = ',')
rna_meta <- read.table(file = 'out/scMVP_input/rna_metadata.csv', header = TRUE, sep = ',')
atac_meta[1:2,]
rna_meta[1:2,]
joint_meta <- atac_meta[match(joint[,1], atac[,1]), c(1,4,5,6)]
joint_meta[1:2,]
# check if the l1 col in joint_meta is the same as what I would retrive from rna_meta
identical(joint_meta$predicted.celltype.l1, rna_meta[match(joint[,1], rna[,1]), 6])
# TRUE, great
joint_meta$seurat_clusters <- rna_meta[match(joint[,1], rna[,1]), 5]
joint_meta[1:2,]

# write out
writeMM(joint_atacMM, 'out/scMVP_input/atac_matrix.mtx')
writeMM(joint_rnaMM, 'out/scMVP_input/rna_matrix.mtx')
write.table(joint, file = 'out/scMVP_input/joint_barcodes.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(joint_meta, file = 'out/scMVP_input/joint_metadata.tsv', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)


##
# to color cells that are in RNA but not in ATAC on nFeature vs nCount colored by mt% gradient plot

data <- readRDS(file = 'rds/seurat_pbmc_rna_filtered_clustered.rds')

joint_barcodes <- read.table(file = 'out/scMVP_input/joint_barcodes.tsv', header = FALSE)

data$joint <- 1.0# these are missing cells
data$joint[which(joint_barcodes[,1] %in% colnames(data))] <- 0.3
table(data$joint)
# 9213: 0.3; 1125: 1

data@meta.data %>% ggplot(aes(x=nCount_RNA,y=nFeature_RNA,color=joint))+
  geom_point(alpha=0.5)+
  scale_color_gradient(low='grey',high='blue')+
  stat_smooth(method=lm)+
  scale_x_log10()+
  scale_y_log10()+
  theme_classic()+
  geom_vline(xintercept=500)+
  geom_hline(yintercept=250)+
  facet_wrap(~orig.ident)

# nah, I think I need to keep the bad QC cells

data <- Read10X(data.dir = "raw/filtered_feature_bc_matrix_gz/")
data <- CreateSeuratObject(counts = data$`Gene Expression`, project = "pbmc", min.cells = 1, min.features = 0)
data$percent_mt <- PercentageFeatureSet(data, pattern = '^MT-')

data$joint <- 1.0# these are missing cells
data$joint[which(joint_barcodes[,1] %in% colnames(data))] <- 0.3
table(data$joint)
# 9213: 0.3; 1125: 1

data@meta.data %>% ggplot(aes(x=nCount_RNA,y=nFeature_RNA,color=joint))+
  geom_point(alpha=0.5)+
  scale_color_gradient(low='grey',high='blue')+
  stat_smooth(method=lm)+
  scale_x_log10()+
  scale_y_log10()+
  theme_classic()+
  geom_vline(xintercept=500)+
  geom_hline(yintercept=250)+
  facet_wrap(~orig.ident)


##
# I want to find out how to add umap coordinates from Portal to a seurat obj

# look at seurat data structure
#data <- readRDS(file = 'rds/seurat_pbmc_rna_filtered_clustered.rds')
str(data)
# cannot find where I can insert my umap coordinate data for 20 dims
# also I have nCell(RNA) + nCell(ATAC); the barcodes will not match
# can I create my own seurat obj?
library(Seurat)
library(dplyr)
library(Matrix)
atacMM <- readMM(file='out/atac/atac_matrix.mtx')
rnaMM <- readMM(file='out/rna/rna_matrix.mtx')
dim(atacMM)# 24919 9939
dim(rnaMM)# 29717 10338

atac_barcode <- read.table(file = 'out/atac/atac_barcodes.tsv', header = FALSE)
atac_barcode <- atac_barcode$V1
rna_barcode <- read.table(file = 'out/rna/rna_barcodes.tsv', header = FALSE)
rna_barcode <- rna_barcode$V1

atac_feature <- read.table(file = 'out/atac/atac_features.tsv', header = FALSE)
atac_feature <- atac_feature$V1
rna_feature <- read.table(file = 'out/rna/rna_features.tsv', header = FALSE)
rna_feature <- rna_feature$V1

atac_meta <- read.table(file = 'out/atac/atac_metadata.csv', header = TRUE, sep = ',')
atac_meta[1:4,]
atac_meta <- atac_meta[,c(1,4,5,6)]
rna_meta <- read.table(file = 'out/rna/rna_metadata.csv', header = TRUE, sep = ',')
rna_meta[1:4,]
rna_meta <- rna_meta[,c(1,5,6,7)]
# rename seurat_clusters to Clusters, so seurat_clusters will be added to ArchR Clusters in one col
colnames(rna_meta)
colnames(rna_meta) <- c('barcode','Clusters','predicted.celltype.l1','predicted.celltype.l2')

atac <- CreateSeuratObject(
  counts = atacMM,
  project = 'pbmc',
  assay = 'atac',
  #row.names = atac_feature,#feature names
  meta.data = atac_meta
)
# No cell names (colnames) names present in the input matrix
# add cell names to count matrix
colnames(atacMM) <- atac_barcode
# No feature names (rownames) names present in the input matrix
rownames(atacMM) <- atac_feature
# Some cells in meta.data not present in provided counts matrix
dim(atac_meta)# 9939 8
# that's weird, there are 9939 cells in meta and mtx
identical(atac_meta$barcode, colnames(atac))
# TRUE
# then they are identical and should be all present

# same for creating seurat obj for rna
colnames(rnaMM) <- rna_barcode
rownames(rnaMM) <- rna_feature
rna <- CreateSeuratObject(
  counts = rnaMM,
  project = 'pbmc',
  assay = 'rna',
  meta.data = rna_meta
)

# merging two seurat obj? they have diff number of genes
# Merging Seurat objects will add zero counts 
# for features present in one object but not another
# https://github.com/satijalab/seurat/issues/1676
merged <- merge(
  x = atac,
  y = rna,
  add.cell.ids = c('atac', 'rna'),
  project = 'pbmc'
)

merged
merged@meta.data[1:4,]
# but the barcodes have prefix now
# what if I don't add.cell.ids?
# now barcodes have suffix for repetitive cells
# besides, cells are in different assays

# but for metadata, they are merged
# missing metadata cols in either assay is left empty
# which means I can still add umap coordinates
# I will also remove unwanted metadata cols in the prev lines, to keep things simple
merged@meta.data[1:4,]
nrow(merged@meta.data)
merged@meta.data[20274:20277,]
# seems CreateSeuratObject() auto added nCount and nFeature, okay

# import umap coordinates
coord <- read.table(file = 'out/portal_output/test_latent.csv', header = TRUE, sep = ',')
coord[1:2,]
dim(coord)
coord <- coord[,-1]
coord[1:2,]
colnames(coord) <- paste0('dim',1:20)
coord[1:2,]

# add to metadata
merged@meta.data <- cbind(merged@meta.data, coord)

# wait, I can't actually viz this from metadata
# remove
merged@meta.data <- merged@meta.data[,1:9]
colnames(merged@meta.data)

# I need to follow this: https://github.com/satijalab/seurat/issues/5113
# convert to matrix
coord_mat <- as(coord, 'matrix')
coord_mat[1:4,]
#colnames(UMAP_coordinates) <- gsub(pattern = "\\.", replacement = "_", x = colnames(UMAP_coordinates))

merged[['UMAP']] <- CreateDimReducObject(
  embeddings = coord_mat[1:10338,],
  key = 'dim',
  global = T,
  assay = 'rna'# since I have to specify an assay, I'd trim the atac assay cells
)

merged[['UMAP']] <- CreateDimReducObject(
  embeddings = coord_mat[10339:20277,],
  key = 'dim',
  global = T,
  assay = 'atac'
)

str(merged)

DimPlot(merged, dims = c('dim_1','dim_2'), reduction = 'UMAP', label = TRUE)
# Error: subscript out of bounds

# stuck
merged_umap <- RunUMAP(merged, dims = 1:20, reduction = 'UMAP')
merged_umap <- FindNeighbors(object = merged_umap, reduction = 'UMAP', dims = 1:20, verbose = FALSE)
# Please provide rownames (cell names) with the input object
# but row.names(merged_umap) has values!
merged_umap <- FindClusters(merged_umap, verbose = FALSE)
merged_umap
DimPlot(merged_umap, reduction = 'umap', label = TRUE)
DimPlot(merged, dims = c('dim_1','dim_2'), reduction = 'UMAP', label = TRUE)
# Error: subscript out of bounds
DimPlot(merged, dims = 1:2, reduction = 'UMAP', label = TRUE)
# Error: subscript out of bounds

ggplot(data=coord[,1:2], mapping = aes(x=dim1, y=dim2)) + geom_point()
ggplot(data=coord[,3:4], mapping = aes(x=dim3, y=dim4)) + geom_point()

library(Azimuth)
merged <- RunAzimuth(merged, assay = 'rna', reference = 'pbmcref')

DefaultAssay(merged) <- 'rna'
merged <- SCTransform(merged, assay = 'rna', verbose = FALSE)


# check azimuth clusters

library(Seurat)
library(Azimuth)
library(patchwork)
library(dplyr)

setwd('~/Labs/wijst/')

data <- Read10X(data.dir = "raw/filtered_feature_bc_matrix_gz/")
pbmc <- CreateSeuratObject(counts = data$`Gene Expression`, project = "pbmc", min.cells = 1, min.features = 0)
pbmc$percent_mt <- PercentageFeatureSet(pbmc, pattern = '^MT-')

rna_barcode <- read.table(file = 'out/rna/rna_barcodes.tsv', header = FALSE)
rna_barcode <- rna_barcode$V1

rna_barcode[1:4]
pbmc@assays$RNA@counts@Dimnames[[2]][1:4]
colnames(pbmc)[1:4]
# subset(pbmc, subset=which(pbmc@assays$RNA@counts@Dimnames[[2]] %in% rna_barcode))
# Error in FetchData.Seurat(object = object, vars = unique(x = expr.char[vars.use]),
# None of the requested variables were found: rna_barcode
# the var name must be in the seurat obj
colnames(pbmc@meta.data)
colnames(pbmc@meta.data)[1] <- 'barcode'
colnames(pbmc@meta.data)
pbmc$barcode <- colnames(pbmc)
pbmc$barcode[1:4]
subset(pbmc, subset=which(barcode %in% rna_barcode))



DefaultAssay(pbmc) <- 'RNA'
pbmc <- RunAzimuth(pbmc, reference = 'pbmcref')

DimPlot(pbmc, reduction = 'ref.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
colnames(pbmc@meta.data)
DimPlot(pbmc, reduction = 'ref.umap', group.by = 'predicted.celltype.l1', label = TRUE, repel = TRUE, label.size = 3) + NoLegend()



