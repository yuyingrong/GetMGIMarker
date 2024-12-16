
#
#
#

# signac multiome

setwd('~/Labs/wijst/')

# $ wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
# $ wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
# $ wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi

#install.packages('Signac')
#BiocManager::install('EnsDb.Hsapiens.v86')
# this is an ensembl based annotation package
#BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
# kept failing to install "ade4", so ran the two lines below
# but without installing "ade4", all packages can load
#install.packages("remotes")
#remotes::install_github("sdray/ade4")
#

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1)

# load RNA and ATAC data
counts <- Read10X_h5('raw/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5')

#install.packages('hdf5r')
#library(hdf5r)
fragpath <- 'raw/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
# Please install biovizBase
# BiocManager::install("biovizBase")
# library(biovizBase)
seqlevelsStyle(annotation) <- 'UCSC'

# create a Seurat object containing the RNA data
pbmc <- CreateSeuratObject(
	counts = counts$`Gene Expression`,
	assay = 'RNA'
)

colnames(pbmc@meta.data)
# "orig.ident" "nCount_RNA" "nFeature_RNA"

# create ATAC assay and add it to the object
pbmc[['ATAC']] <- CreateChromatinAssay(
	counts = counts$Peaks, 
	sep = c(':', '-'), 
	fragments = fragpath, 
	annotation = annotation
)
colnames(pbmc@meta.data)
# "orig.ident" "nCount_RNA" "nFeature_RNA" "nCount_ATAC" "nFeature_ATAC"


pbmc
# 36601 genes, 11898 cells



#
#
### Quality control

# compute per-cell QC metrics using DNA accessibility data
# removes cells that rae outliers for these metrics
# also remove cells with low or unusually high counts for RNA/ATAC assay
DefaultAssay(pbmc) <- 'ATAC'

pbmc <- NucleosomeSignal(pbmc)
# Calculate the strength of the nucleosome signal per cell. Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to fragments < 147 bp (nucleosome-free)
pbmc <- TSSEnrichment(pbmc)


colnames(pbmc@meta.data)
# "orig.ident" "nCount_RNA" "nFeature_RNA" "nCount_ATAC" "nFeature_ATAC"
# "nucleosome_signal" "nucleosome_percentile" "TSS.enrichment" "TSS.percentile"


VlnPlot(
	object = pbmc, 
	features = c('nCount_RNA', 'nCount_ATAC', 'nucleosome_signal', 'TSS.enrichment'), 
	ncol = 4, 
	pt.size = 0
)

# instead of subset(), I will add QC scores to a col in metadata
'
pbmc <- subset(
	x = pbmc, 
	subset = nCount_RNA < 25000 & 
	nCount_RNA > 1000 & 
	nCount_ATAC < 10000 & 
	nCount_ATAC > 1000 & 
	nucleosome_signal < 2 & 
	TSS.enrichment > 1
)
'

# why should nucleosome signal be <2?
# nucleosome signal is (mono-nucleosome frag) to (nucleosome-free frag) ratio ...
# why would TSS.enrichment be nearly 5?
# this TSS.e is (signal from TSS+/-1000bp region) to (signal from background region) ratio
# why take only those >1?
# if not >1, then the chromatin openness is random and the cell is bad quality
# assumption is that TSS regions should always be open

pbmc


#
#
### Peak Calling

# here I call peaks on all cells together
# if I identified cell types, I can use "group.by" to id peaks for each cell type
# then I can id peaks specific to rare cell populations

# call peaks using MACS2
# peaks <- CallPeaks(pbmc)
# Error in CallPeaks.default(object = allfragpaths, macs2.path = macs2.path,  : 
# MACS2 not found. Please install MACS:https://macs3-project.github.io/MACS/
# but MACS2 is python; how can I install on cmdline for python and make it talk to R?

# enter a path for macs2 (I actually installed MACS3)
peaks <- CallPeaks(object = pbmc, macs2.path = '/Users/shell/miniforge3/bin/MACS3')

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc), 
  features = peaks, 
  cells = colnames(pbmc)
)

# create a new assay using the MACS peak set and add it to the Seurat obj
pbmc[['peaks']] <- CreateChromatinAssay(
  counts = macs2_counts, 
  fragments = fragpath, 
  annotation = annotation
)

#
# 
## Gene expression data processing
DefaultAssay(pbmc) <- 'RNA'
pbmc <- SCTransform(pbmc)
# outliers are found and ignored for later analysis
pbmc <- RunPCA(pbmc)

# making my own UMAP clusters
pbmc
# 4 dimensional reductions calculated: pca, integrated_dr, ref.umap, lsi
DimPlot(pbmc, reduction = 'pca', label = TRUE)
DimPlot(pbmc, reduction = 'ref.umap', label = TRUE)
# following the SCTransform tutorial: https://satijalab.org/seurat/articles/sctransform_vignette.html
DefaultAssay(pbmc) <- 'SCT'
pbmc <- RunUMAP(pbmc, dims = 1:50)
pbmc <- FindNeighbors(pbmc, dims = 1:50)
pbmc <- FindClusters(pbmc)
pbmc
# 5 dimensional reductions calculated: pca, integrated_dr, ref.umap, lsi, umap
DimPlot(pbmc, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
DimPlot(pbmc, reduction = 'ref.umap', group.by = 'seurat_clusters', label = TRUE)
DimPlot(pbmc, reduction = 'umap', group.by = 'predicted.celltype.l1', label = TRUE)
DimPlot(pbmc, reduction = 'umap', group.by = 'predicted.celltype.l2', label = TRUE)

# DNA accessibility data processing
DefaultAssay(pbmc) <- 'peaks'
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc) #term-frequency inverse-document frequency, normalization
pbmc <- RunSVD(pbmc)

#saveRDS(pbmc, file='rds/signac_pbmc_multiome_unfiltered.rds')
dim(pbmc)
# 36601 1404
# wait, that looks weird
# how can there be only 1404 cells??

# after re-run without subset(): 132291 11898, good


#
#
## Annotating cell types

# annotating using Azimuth
# following: https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html

devtools::install_github('satijalab/seurat-data')
devtools::install_github('satijalab/azimuth')

library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

# this is trial dataset, but I don't have to use!!
# InstallData('pbmcsca')
# pbmcsca <- LoadData('pbmcsca')
# rm(pbmcsca)
# detach("package:pbmcsca.SeuratData", unload=TRUE)

DefaultAssay(pbmc) <- 'RNA'
pbmc <- RunAzimuth(pbmc, reference = 'pbmcref')

DimPlot(pbmc, reduction = 'ref.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
colnames(pbmc@meta.data)
DimPlot(pbmc, reduction = 'ref.umap', group.by = 'predicted.celltype.l1', label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
DimPlot(pbmc, reduction = 'ref.umap', group.by = 'predicted.celltype.l3', label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

# 
#
## mapping using Seurat v4
# following: https://stuartlab.org/signac/articles/pbmc_multiomic.html#annotating-cell-types

library(SeuratDisk)

reference <- LoadH5Seurat('ref/pbmc_multimodal.h5seurat')
# this is a CITE-seq dataset

DimPlot(
  object = reference, 
  reduction = 'wnn.umap', 
  group.by = 'celltype.l2', 
  label = TRUE, 
  label.size = 3, 
  repel = TRUE
) + NoLegend()

# find anchors between ref and query
# decided to not pursue this path

#DefaultAssay(pbmc) <- 'SCT'

#anchors <- FindTransferAnchors(
#  reference = reference, 
#  query = pbmc, 
#  normalization.method = 'SCT', 
#  reference.reduction = 'spca', # recommended to use supervised-PCA for CITE-seq datasets
#  recompute.residuls = FALSE, 
#  dims = 1:50
#)

#predictions <- TransferData(
#  anchorset = anchors, 
#  refdata = reference$celltype.l2, 
#  weight.reduction = pbmc[['pca']], 
#  dims = 1:50
#)

# make sure you don't over-write the celltype.l2 predictions from RNA and Azimuth

#pbmc <- AddMetaData(
#  object = pbmc, 
#  metadata = predictions
#)

#Idents(pbmc) <- 'predicted.id'

#
#
## visualize cell quality in DimPlot
pbmc <- readRDS('rds/signac_pbmc_multiome_unfiltered.rds')
DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  pbmc, 
  features = c('percent.mt', 'nFeature_RNA', 'TSS.enrichment', 'nucleosome_signal'), 
  label = FALSE
)

#
#
## import from ArchR clusters based on ATAC analysis

# check rna metadata
dim(pbmc@meta.data)
# 11898 23
rownames(pbmc@meta.data)[1:3]
colnames(pbmc@meta.data)
pbmc@meta.data[1,]
pbmc@meta.data$X <- rownames(pbmc@meta.data)

# import atac clusters from saved csv from ArchR
metadata_atac <- read.csv(file = 'out/metadata_ATAC.csv', header = TRUE)
dim(metadata_atac)
# 10846 16
colnames(metadata_atac)
rownames(metadata_atac)[1:3]
metadata_atac[1,]
rownames(metadata_atac) <- metadata_atac$X

library(dplyr)
dim(anti_join(metadata_atac, pbmc@meta.data, by='X'))
# 574 16

length(intersect(metadata_atac$X, pbmc@meta.data$X))
# 10272
# 574 + 10272 = 10846, # of barcodes in atac metadata

# are metadata barcodes a subset of seurat obj barcodes?
all(metadata_atac$X %in% pbmc@meta.data$X)
# FALSE, i.e. not all values are TRUE for if each item in atac is also in pbmc
all(pbmc@meta.data$X %in% metadata_atac$X)
# FALSE

dim(anti_join(pbmc@meta.data, metadata_atac, by='X'))
# 1626 24

metadata_rna <- pbmc@meta.data

metadata_merged <- merge(metadata_rna, metadata_atac, by = 'X', all = TRUE)
dim(metadata_merged)
# 12472 39
# but originally the seurat obj only has 11898 cells!

colnames(metadata_merged)

metadata_rna <- merge(metadata_rna, metadata_merged, by = 'X', all = FALSE)
dim(metadata_rna)
colnames(metadata_rna)
metadata_rna[1,]

# import merged metadata back to seurat obj
pbmc@meta.data <- metadata_rna

# add rownames back to pbmc@meta.data
rownames(pbmc@meta.data) <- pbmc@meta.data$X

# visualize
FeaturePlot(
  pbmc, 
  features = c('TSSEnrichment', 'NucleosomeRatio', 'DoubletScore', 'Clusters'), 
  label = FALSE
)

DimPlot(
  object = pbmc, 
  group.by = 'Clusters', 
  label = TRUE, 
  label.size = 3, 
  repel = TRUE
) + NoLegend()

DimPlot(
  object = pbmc, 
  group.by = 'DoubletScore', 
  label = TRUE, 
  label.size = 3, 
  repel = TRUE
) + NoLegend()

DimPlot(
  object = pbmc, 
  group.by = 'predicted.celltype.l1', 
  label = TRUE, 
  label.size = 3, 
  repel = TRUE
) + NoLegend()


#
# now I want to project rna clusters/celltypes onto the atac clusters

# read in metadata_atac
pbmc <- readRDS('rds/signac_pbmc_multiome_unfiltered.rds')
metadata_atac <- read.csv(file = 'out/metadata_ATAC.csv', header = TRUE)
colnames(metadata_atac)[1] <- 'barcode'
metadata_atac[1,]

# duplicate pbmc@meta.data to a new var
metadata_rna <- pbmc@meta.data
colnames(metadata_rna)[1] <- 'barcode'
metadata_rna$barcode <- rownames(metadata_rna)
metadata_rna[1,]

# (optional: add relevant cols to var rna_clusters)
rna_clusters <- metadata_rna[,c(1:3, 10, 12, 18, 20)]
colnames(rna_clusters)[1] <- 'barcode'
rna_clusters$barcode <- rownames(rna_clusters)
rna_clusters[1,]

# add only rna cluster col to var atac_clusters
# option 1: new var atac_clusters, use merge(), but somehow all cells become NA
atac_clusters <- merge(metadata_atac, rna_clusters, by = 'barcode', all = TRUE)
atac_clusters[1,]
dim(atac_clusters)
# 12472 22
metadata_atac <- merge(metadata_atac, atac_clusters, all = FALSE)
metadata_atac[1,]
dim(metadata_atac)
# randomly check a few cells
metadata_atac$predicted.celltype.l1[c(1:5, 22, 999)]
# write to file
write.csv(metadata_atac, 'out/metadata_merged.csv', row.names=FALSE)

# option 2: match(), map rna barcode index to atac index - DID NOT WORK
# match(c(1,2,3,4), c(1,3,2,4))
# 1 3 2 4
dim(metadata_atac)
# 10846 16
dim(metadata_rna)
# 11898 23
metadata_atac$predicted.celltype.l1 <- metadata_rna$predicted.celltype.l1[match(metadata_rna$barcode, metadata_atac$barcode, nomatch = '0')]
# Error: replacement has 10272 rows, data has 10846






