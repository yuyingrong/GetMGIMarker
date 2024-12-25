
# same-cell single-cell RNA- and ATAC-seq multiomic data integration
# ArchR for ATAC-seq data analysis and manipulation


setwd('~/Labs/wijst')

# install.packages("devtools")
library(ArchR)
# ArchR::installExtraPackages()

# set a known seed to replicate operations requiring randomization
set.seed(1)

# set default n of threads for parallelized operations
# n should match the specifications of your local machine
# you can only set cores to ncores minus 2: so 8-2=6
addArchRThreads(threads=6)

# set input files
# tutorial filename looks like so: “HemeFragments/scATAC_BMMC_R1.fragments.tsv.gz”
# so I would use my 'raw/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'
inputFiles <- getInputFiles(paths = '../raw')
inputFiles

# add reference geome annotation, to load chromosome and gene info
# tutorial uses 'hg19'
# mine is from 'hg38'
addArchRGenome('hg38')

# create arrow files
# this means to read accessible fragments from provided input files
# calc QC info for each cell (TSSE scores and nucleosome info)
# filter cells based on QC parameters
# create a genome-wide TileMatrix using 500-bp bins
# create a GeneScoreMatrix

# QC metric 1: # of unique nuclear fragments (not mitoDNA, not PCR duplicates)
# QC metric 2: signal-to-background ratio: low ratio suggests dead/dying cells with de-chromatinized DNA
# QC metric 3: frag size distribution: nucleosomal periodicity

ArrowFiles <- createArrowFiles(
	inputFiles = inputFiles, 
	sampleNames = names(inputFiles), 
	minTSS = 4, # you can always increase this later
	minFrags = 1000, 
	addTileMat = TRUE, 
	addGeneScoreMat = TRUE, 
	#subThreading = FALSE
)

# Number of Cells Pass Filter = 12162
# Median Frags = 17122
# Median TSS Enrichment = 18.3305

ArrowFiles

# infer doublets
# this step is done to the ArrowFiles before creating proj obj from the ArrowFiles
doubScores <- addDoubletScores(
	input = ArrowFiles, 
	k = 10, # refers to how many cells near a 'pseudo-doublet' to count
	knnMethod = 'UMAP', # refers to the embedding to use for nn search
	LSIMethod = 1
)

# create an ArchRProjet
proj <- ArchRProject(
	ArrowFiles = ArrowFiles, 
	outputDirectory = 'pbmc10x_archr', 
	copyArrows = TRUE
)

proj

colnames(proj@cellColData)
length(which(proj@cellColData$DoubletScore!=0))
# 928
length(which(proj@cellColData$DoubletEnrichment==0))
# 45


# save ArchRProject
saveArchRProject(
	ArchRProj = proj, 
	outputDirectory = getOutputDirectory(proj), 
	overwrite = TRUE, 
	load = TRUE, 
	dropCells = FALSE, 
	logFile = createLogFile('saveArchRProject'), 
	threads = getArchRThreads()
)

# load ArchRProject
proj <- loadArchRProject(path = './pbmc10x_archr', force = FALSE, showLogo = FALSE)







# ask which data matrices are available within the ArchRProject
# later we will start adding to this project
getAvailableMatrices(proj)

# filter putative doublets based on prev determined doublet scores
# cells won't be physically removed from the ArrowFiles but 
# ArchRProject will ignore these cells for later analyses
proj <- filterDoublets(ArchRProj = proj)
# Filtering 1316 cells from ArchRProject!
# pbmc_granulocyte_sorted_10k_atac : 1316 of 12162 (10.8%)
length(which(proj@cellColData$DoubletScore==0))
# 10846 = 12162 - 1316

# implement an iterative LSI dimensionality reduction
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = 'TileMatrix', name = 'IterativeLSI')
# read more about iterative LSI: https://www.archrproject.com/bookdown/iterative-latent-semantic-indexing-lsi.html

# call clusters in this reduced dimension sub-space
proj <- addClusters(input = proj, reducedDims = 'IterativeLSI')

# visualize data using 2D representation such as UMAP
# add a UMAP embedding to ArchRProject
proj <- addUMAP(ArchRProj = proj, reducedDims = 'IterativeLSI')


# visualize various attributes of cells stored in `cellColData` in ArchRProject
# specify the var to use for coloration via `colorBy` and `name` parameters

# color by 'sample'
p1 <- plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'Sample', embedding = 'UMAP')

# color by 'cluster'
p2 <- plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'Clusters', embedding = 'UMAP')

# view plots
ggAlignPlots(p1, p2, type = 'h')
# since I have one sample: 
ggAlignPlots(p2, type = 'h')

# save to PDF
plotPDF(p1, p2, name = 'pbmc10x_umap_clusters.pdf', ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
plotPDF(p2, name = 'pbmc10x_umap_clusters.pdf', ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


#
# 
#
## custom analysis

# access metadata
dim(proj@cellColData)
# 10846 16

# store metadata in a new var
metadata <- proj@cellColData

typeof(metadata)
# S4

# I want just the data.frame
metadata <- as.data.frame(proj@cellColData)

typeof(metadata)
# list

dim(metadata)
# 10846 16

rownames(metadata)[1:3]
# "pbmc_granulocyte_sorted_10k_atac#GGAGTCTGTGAGCAAG-1"
# "pbmc_granulocyte_sorted_10k_atac#GTTAAGTGTCACTCGC-1"
# "pbmc_granulocyte_sorted_10k_atac#GGTTGCGGTAAACAAG-1" 

# remove the string preceding and including #
gsub("^.*?#", "", rownames(metadata)[1])
# see ?regex
# ^: matches beginning of str
# .: matches any char
# *: repeated zero or more times
# #: matches #
# ?: "lazy" match, only matches to the first #
# "": replace matched with ""
gsub("^.*?#", "", rownames(metadata))[1:3]
rownames(metadata) <- gsub("^.*?#", "", rownames(metadata))

#metadata[1,]
#metadata$barcode <- rownames(metadata)
#metadata <- metadata[,2:17]
#metadata[1,]

colnames(metadata)
# "Sample" "TSSEnrichment" "ReadsInTSS" "ReadsInPromoter"
# "ReadsInBlacklist" "PromoterRatio" "PassQC" "NucleosomeRatio"
# "nMultiFrags" "nMonoFrags" "nFrags" "nDiFrags" "DoubletScore"
# "DoubletEnrichment" "BlacklistRatio" "Clusters"         

# remove Sample col
metadata <- metadata[,2:16]
metadata[1,]

# write data.frame to csv
write.csv(metadata, 'out/metadata_ATAC.csv', row.names=TRUE)


#
# now I want to project rna clusters/celltypes onto the atac clusters
proj <- loadArchRProject(path = 'pbmc10x_archr', force = FALSE, showLogo = TRUE)
proj
# nCells: 10846

metadata_merged <- read.csv(file = 'out/metadata_merged.csv', header = TRUE)
metadata_merged[1,]

# shorten rownames of ArchRProject proj
# rownames(proj@cellColData) <- gsub("^.*?#", "", rownames(proj@cellColData))
# proj@cellColData[1,]

# save ArchRProject with shorter rownames - LATER I REVERSED THE STEP

# add prefix back to barcode
metadata_merged$barcode <- paste0("pbmc_granulocyte_sorted_10k_atac#", metadata_merged$barcode)
metadata_merged$barcode[1:4]

# check if prefixed barcode names match with ArchRProject cellnames
identical(rownames(proj@cellColData), metadata_merged$barcode)
# FALSE
all(rownames(proj@cellColData) %in% metadata_merged$barcode)
# TRUE
# so the order is different, but names are identical

proj@cellColData$predicted.celltype.l1 <- metadata_merged$predicted.celltype.l1[match(rownames(proj@cellColData), metadata_merged$barcode)]
proj@cellColData$predicted.celltype.l2 <- metadata_merged$predicted.celltype.l2[match(rownames(proj@cellColData), metadata_merged$barcode)]

proj@cellColData[1,]

plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'predicted.celltype.l1', embedding = 'UMAP')
# Error: Not all cells in embedding are present in ArchRProject!

# why error? can I still plot sth that I could plot earlier?
plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'Clusters', embedding = 'UMAP')
# Error: Not all cells in embedding are present in ArchRProject!
# okay, so I guess this is due to I changed cellnames?

str(proj)
# I still see several places where cellnames are in this format
# "pbmc_granulocyte_sorted_10k_atac#GGAGTCTGTGAGCAAG-1"

# for expediency, I would just add the string back to cellnames in CellColData
rownames(proj@cellColData) <- paste0("pbmc_granulocyte_sorted_10k_atac#", rownames(proj@cellColData))
proj@cellColData[1,]

# now plot
plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'predicted.celltype.l1', embedding = 'UMAP')
plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'predicted.celltype.l2', embedding = 'UMAP')
plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'TSSEnrichment', embedding = 'UMAP')
plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'Clusters', embedding = 'UMAP')

# what are the NA's in the plot?
#proj@cellColData$predicted.celltype.l1 <- metadata_merged$predicted.celltype.l1[match(rownames(proj@cellColData), metadata_merged$barcode)]
#c(10,20,40,50,NA) == c(10,20,30,40,50)[match(c(1,2,3,4,5), c(1,2,2,3,4))]


##
# updated 20230301
setwd('~/Labs/wijst')
library(ArchR)
set.seed(1)

# load proj
proj <- loadArchRProject(path = './pbmc10x_archr', force = FALSE, showLogo = FALSE)

# goal is to filter cells also by FRIP (fraction of reads in peaks)
# perhaps I can add this FeatureMatrix to get FRIP info
# right now I don't have this matrix

getAvailableMatrices(proj)
# "GeneScoreMatrix" "TileMatrix"

# following: https://github.com/GreenleafLab/ArchR/blob/f6c0388bd37023400794c9ae8562ad69e3ba9fd7/R/MatrixFeatures.R
# addFeatureMatrix(input = proj, threads = 6, force = FALSE)
# however, for the 2nd arg, features, Idk what to input for a GRanges obj

# trying this: https://www.archrproject.com/reference/getGroupSE.html
getGroupSE(ArchRProj = proj, useMatrix = 'GeneScoreMatrix')
# nah, I think that would be over complicating

# calling peaks
findMacs2()
# could not find; I'll supply the path
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = 'Clusters', 
  pathToMacs2 = '/Users/shell/miniforge3/bin/MACS3'
)
# Error: No Coverage Metadata found for: Clusters. Please run addGroupCoverages!

proj <- addGroupCoverages(ArchRProj = proj, groupBy = 'Clusters')

# then re-run addReproduciblePeakSet()

# to retrieve peak set as a GRanges obj
getPeakSet(proj)

# now I actually know what to input for adding a feature matrix
proj <- addFeatureMatrix(
  input = proj, 
  features = getPeakSet(proj), 
  threads = 6, 
  force = FALSE
)
# Error: Matrix Group Already Exists! Set force = TRUE to overwrite!
# When did it come to exist?
getAvailableMatrices(proj)
# "FeatureMatrix"   "GeneScoreMatrix" "TileMatrix"
# ok, somehow it does exist now

getGroupSE(
  ArchRProj = proj, 
  useMatrix = 'GeneScoreMatrix', 
)
# colData names(15): TSSEnrichment ReadsInTSS ... BlacklistRatio nCells
# there's still no FRIP in colData names!!!
# according to: https://www.archrproject.com/bookdown/identification-of-positive-tf-regulators.html
# FRIP should have appeared, why?

# what if I force TRUE
proj <- addFeatureMatrix(
  input = proj, 
  features = getPeakSet(proj), 
  threads = 6, 
  force = TRUE
)

# oh, I noticed that I looked at GeneScoreMatrix earlier, changed that
# also grouped by clusters, otherwise it's grouped by sample
seGroupMotif <- getGroupSE(
  ArchRProj = proj, 
  useMatrix = 'MotifMatrix',
  groupBy = 'Clusters'
)
# Error: Object 'MotifMatrix/Info/FeatureDF' does not exist in this HDF5 file.
# according to: https://github.com/GreenleafLab/ArchR/issues/823
# addDeviationsMatrix() and even more steps needs to be run

# did not create new proj suggested by https://www.archrproject.com/bookdown/add-peak-matrix.html
proj <- addPeakMatrix(proj)

getAvailableMatrices(proj)

seGroupMotif <- getGroupSE(
  ArchRProj = proj, 
  useMatrix = 'PeakMatrix',
  groupBy = 'Clusters'
)

seGroupMotif

# FRIP appeared
# now I don't need to add DeviationsMatrix anymore!

seGroupMotif$FRIP[1:5]

# save
saveArchRProject(
  ArchRProj = proj, 
  outputDirectory = getOutputDirectory(proj), 
  overwrite = TRUE, 
  load = TRUE, 
  dropCells = FALSE, 
  logFile = createLogFile('saveArchRProject'), 
  threads = getArchRThreads()
)


#
# Now performing cell filtering

# remove doublets
all(proj@cellColData$DoubletScore==0)
# TRUE
# doublets are already removed

# remove too small TSS Enrichment
all(proj@cellColData$TSSEnrichment>=4)
# TRUE
all(proj@cellColData$TSSEnrichment>=5)
# so signal to bg ratio is at least >= 4

# remove FRiP
# I saw in MAESTRO tutorial that 0.2 was used as FRiP cutoff for PBMC cells
all(proj@cellColData$FRIP>=0.2)
# FALSE
summary(proj@cellColData$FRIP)

ggplot(proj@cellColData, aes(x=FRIP))+
  geom_histogram(binwidth=0.01)+
  geom_vline(xintercept=0.2, linetype='dashed', color='blue')+
  geom_text(aes(x=0.2, y=-20, label='cutoff = 0.2'))+
  ggtitle('FRiP')
# Error: data must be a data.frame
typeof(proj@cellColData)
# S4
typeof(data@meta.data)# if you have seurat open
# list

# add as.data.frame()
ggplot(as.data.frame(proj@cellColData), aes(x=FRIP))+
  geom_histogram(binwidth=0.01)+
  #geom_vline(xintercept=0.35, linetype='dashed', color='blue')+
  #geom_text(aes(x=0.35, y=-20, label='cutoff = 0.35'))+
  # prof suggested to cutoff at 0.5
  geom_vline(xintercept=0.5, linetype='dashed', color='orange')+
  geom_text(aes(x=0.5, y=-20, label='cutoff = 0.5'))+
  ggtitle('FRiP')

# 0.35 looks like a better cutoff than 0.2
length(which(proj@cellColData$FRIP<0.35))
# 613, cutoff 613 cells
length(which(proj@cellColData$FRIP<0.5))
# 907
# >=0.5 is our FRiP cutoff

# instead of using subset() which might cause "not all cells in embedding" error later
# I would just label these cells as "NO" in PassQC or somewhere

# what is in PassQC anyway?
summary(proj@cellColData$PassQC)
# all 1's!

# actually I don't need to label anything based on FRiP, the info is in FRiP already
# at most I can make a list of barcodes without these cells with FRIP < 0.35


#
# create cell barcode metadata files
# first, extract all atac-seq cell barcodes
# since I already have this file from doublet labeling
#atac <- read.csv('out/archr_predictions.csv', header=FALSE)
atac <- read.csv('out/archr_doublet.csv', header=TRUE)
nrow(atac)
# 12163? should be 12162
atac[1:4,]
# found empty 1st row; remove
atac <- atac[-1,]

# remove cell barcode prefixes
atac[,1] <- gsub("^.*?#", "", atac[,1])
colnames(atac) <- c('barcode', 'doublet')

# replace (True with 1) and (False with 0) for doublet
table(atac$doublet)
# False  True 
# 10846  1316
atac$doublet[atac$doublet == 'True'] <- 1
atac$doublet[atac$doublet == 'False'] <- 0
table(atac$doublet)

# update doublet csv
#write.csv(atac, 'out/archr_doublet.csv', row.names=FALSE)
#atac <- read.csv('out/archr_doublet.csv', header=TRUE)

# add cols
colnames(proj@cellColData)

# unnecessary: atac$TSSEnrichment <- rep(-1, nrow(atac))

atac$TSSEnrichment <- proj@cellColData$TSSEnrichment[match(paste0("pbmc_granulocyte_sorted_10k_atac#", atac$barcode), rownames(proj@cellColData))]
atac$Clusters <- proj@cellColData$Clusters[match(paste0("pbmc_granulocyte_sorted_10k_atac#", atac$barcode), rownames(proj@cellColData))]
atac$predicted.celltype.l1 <- proj@cellColData$predicted.celltype.l1[match(paste0("pbmc_granulocyte_sorted_10k_atac#", atac$barcode), rownames(proj@cellColData))]
atac$predicted.celltype.l2 <- proj@cellColData$predicted.celltype.l2[match(paste0("pbmc_granulocyte_sorted_10k_atac#", atac$barcode), rownames(proj@cellColData))]
atac$FRIP <- proj@cellColData$FRIP[match(paste0("pbmc_granulocyte_sorted_10k_atac#", atac$barcode), rownames(proj@cellColData))]
atac$nFrags <- proj@cellColData$nFrags[match(paste0("pbmc_granulocyte_sorted_10k_atac#", atac$barcode), rownames(proj@cellColData))]

summary(atac$nFrags)
# min is 1000 unique number of fragments mapped to ref genome

typeof(atac$TSSEnrichment)
# typeof(atac$TSSEnrichment) looked different from the other attributes
# the others are character, whereas typeof(atac$TSSEnrichment) is double
# so I just changed it
atac$TSSEnrichment <- as.character(atac$TSSEnrichment)
atac$TSSEnrichment[1:5]
atac$predicted.celltype.l1[1:5]
# but FRIP is a double, too
# change it back
atac$TSSEnrichment <- as.double(atac$TSSEnrichment)

atac[1:4,]

# write to csv
write.csv(atac, 'out/atac_metadata.csv', row.names=FALSE)
# note: later changed path

# NOTE!!!
# the file can be improved by: not setting a cutoff bar for nFrags and TSSE, so I know how many cells are filtered in this step
# also by: before filtering doublets, output the TSSE and nFrag info so the doublets have them

# retrive table
atac <- read.csv('out/atac/atac_metadata_unfiltered.csv', header=TRUE)


#
# plotting the TSS by unique_frags plot

# these are not what I wanted
#plotTSSEnrichment(ArchRProj = proj)
#plotFragmentSizes(ArchRProj = proj)

# modified from https://github.com/GreenleafLab/ArchR/blob/master/R/CreateArrow.R

# after doublet removal

ggtitle <- sprintf("%s\n%s\n%s",
  paste0("\nnCells Pass Filter = ", nrow(proj@cellColData)),
  paste0("Median Frags = ", median(proj@cellColData$nFrags)),
  paste0("Median TSS Enrichment = ", median(proj@cellColData$TSSEnrichment))
)

filterTSS <- 4
filterFrags <- 1000

ggPoint(
  x = pmin(log10(proj@cellColData$nFrags), 5) + rnorm(length(proj@cellColData$nFrags), sd = 0.00001),
  y = proj@cellColData$TSSEnrichment + rnorm(length(proj@cellColData$nFrags), sd = 0.00001), 
  colorDensity = FALSE,
  xlim = c(2.5, 5),
  ylim = c(0, max(proj@cellColData$TSSEnrichment) * 1.05),
  baseSize = 12,
  size = 0.2, 
  continuousSet = "sambaNight",
  xlabel = "Log 10 (Unique Fragments)",
  ylabel = "TSS Enrichment",
  title = ggtitle,
  rastr = TRUE) + 
  geom_hline(yintercept=filterTSS, lty = "dashed", size = 0.25) +
  geom_vline(xintercept=log10(filterFrags), lty = "dashed", size = 0.25)

# after FRiP filtering
ggtitle <- sprintf("%s\n%s\n%s",
                   paste0("\nnCells Pass Filter = ", length(which(proj@cellColData$FRIP>=0.5))),
                   paste0("Median Frags = ", median(proj@cellColData$nFrags[proj@cellColData$FRIP>=0.5])),
                   paste0("Median TSS Enrichment = ", median(proj@cellColData$TSSEnrichment[proj@cellColData$FRIP>=0.5]))
)

ggPoint(
  x = pmin(log10(proj@cellColData$nFrags[proj@cellColData$FRIP>=0.5]), 5) + rnorm(length(proj@cellColData$nFrags[proj@cellColData$FRIP>=0.5]), sd = 0.00001),
  y = proj@cellColData$TSSEnrichment[proj@cellColData$FRIP>=0.5] + rnorm(length(proj@cellColData$nFrags[proj@cellColData$FRIP>=0.5]), sd = 0.00001), 
  colorDensity = TRUE,
  xlim = c(2.5, 5),
  ylim = c(0, max(proj@cellColData$TSSEnrichment[proj@cellColData$FRIP>=0.5]) * 1.05),
  baseSize = 12,
  size = 0.2, 
  continuousSet = "sambaNight",
  xlabel = "Log 10 (Unique Fragments)",
  ylabel = "TSS Enrichment",
  title = ggtitle,
  rastr = TRUE) + 
  geom_hline(yintercept=filterTSS, lty = "dashed", size = 0.25) +
  geom_vline(xintercept=log10(filterFrags), lty = "dashed", size = 0.25)

# Wow! after filtering by FRiP, the smaller cluster completely disappeared
# almost as if a TSSEnrichment >10 filter is applied
# but TSSEnrichment is somewhat related to fraction in peaks, so that's kinda expected

# but why can I still see TSSEnrichment below 4?
summary(proj@cellColData$TSSEnrichment)
# Min = 4.024

summary(proj@cellColData$TSSEnrichment[proj@cellColData$FRIP>=0.5])
# Min = 10.33, other params don't change much
# basically with the small hill in FRiP barplot is cutoff, the small hill in TSSEnrichment is also cutoff

###### needs to discuss with prof whether to apply log10(nFrags) >= 3.5

# also quickly plot celltypes over the TSSEnrichment vs log10(nFrags) plot
ggPoint(
  x = pmin(log10(proj@cellColData$nFrags[proj@cellColData$FRIP>=0.5]), 5) + rnorm(length(proj@cellColData$nFrags[proj@cellColData$FRIP>=0.5]), sd = 0.00001),
  y = proj@cellColData$TSSEnrichment[proj@cellColData$FRIP>=0.5] + rnorm(length(proj@cellColData$nFrags[proj@cellColData$FRIP>=0.5]), sd = 0.00001), 
  colorDensity = FALSE,
  color = proj@cellColData$predicted.celltype.l2[proj@cellColData$FRIP>=0.5],
  xlim = c(2.5, 5),
  ylim = c(0, max(proj@cellColData$TSSEnrichment[proj@cellColData$FRIP>=0.5]) * 1.05),
  baseSize = 12,
  size = 0.2, 
  continuousSet = "sambaNight",
  xlabel = "Log 10 (Unique Fragments)",
  ylabel = "TSS Enrichment",
  rastr = TRUE)


#
# export barcodes, features, and count matrix

str(proj)

getAvailableMatrices(proj)

featureSE <- getGroupSE(
  ArchRProj = proj, 
  useMatrix = 'FeatureMatrix', 
)

# counts for each feature per cell in the sample
feature_mtx <- getMatrixFromProject(
  ArchRProj = proj, 
  useMatrix = 'FeatureMatrix'
)

expSE <- getGroupSE(
  ArchRProj = proj, 
  useMatrix = 'GeneScoreMatrix', 
)

expSE_byClusters <- getGroupSE(
  ArchRProj = proj, 
  useMatrix = 'GeneScoreMatrix', 
  groupBy = 'Clusters'
)

exp_mtx <- getMatrixFromProject(
  ArchRProj = proj, 
  useMatrix = 'GeneScoreMatrix'
)

peakSE <- getGroupSE(
  ArchRProj = proj, 
  useMatrix = 'PeakMatrix', 
)

# we decided feature_mtx is the input for Portal

# finally found the data I need
dim(feature_mtx@assays@data$FeatureMatrix)
# 164675  10846
# 164675 is the number of genomic intervals
# 10846 is nCell: all cells minus doublets
# A GRanges object contains a set of genomic intervals
# These will form the rows of the matrix, with each entry
# recording the number of unique reads falling in the
# genomic region for each cell
feature_mtx@assays@data$FeatureMatrix

typeof(feature_mtx@assays@data$FeatureMatrix)
# S4

feature_mtx@assays@data$FeatureMatrix[1:4, 1:4]
colnames(feature_mtx@assays@data$FeatureMatrix) <- gsub("^.*?#", "", colnames(feature_mtx@assays@data$FeatureMatrix))
colnames(feature_mtx@assays@data$FeatureMatrix)[1:2]

write.csv(as.matrix(feature_mtx@assays@data$FeatureMatrix), 'out/atac_feature_matrix.csv', row.names=FALSE)
# cannot coerce class ‘structure("dgCMatrix", package = "Matrix")’ to a data.frame
# coerce into matrix succeeded

# however, the file is getting too large!
# need to export to a sparse matrix instead
# see more in "20230206seurat.R"
library(Matrix)

# after loading the ArchR proj
feature_mtx <- getMatrixFromProject(
  ArchRProj = proj, 
  useMatrix = 'FeatureMatrix'
)

copy <- feature_mtx

# check what is inside
str(feature_mtx@assays@data$FeatureMatrix)
# 164675 10846

# remove cell barcode prefix
# feature_mtx@assays@data$FeatureMatrix@Dimnames[2] <- gsub("^.*?#", "", feature_mtx@assays@data$FeatureMatrix@Dimnames[2])
# this didn't work
colnames(feature_mtx@assays@data$FeatureMatrix) <- gsub("^.*?#", "", colnames(feature_mtx@assays@data$FeatureMatrix))
# colnames(feature_mtx@assays@data$FeatureMatrix)[1:4]
# don't check colnames this way!
# use str()

# write a matrix market file
writeMM(feature_mtx@assays@data$FeatureMatrix, 'out/atac/atac_matrix.mtx')

# create a list of barcodes
write.table(feature_mtx@assays@data$FeatureMatrix@Dimnames[2][[1]], file = 'out/atac/atac_barcodes.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)
# check if write out are correct
feature_mtx@assays@data$FeatureMatrix@Dimnames[2][[1]][1:4]
# good

# since FeatureMatrix has no rownames for genomic intervals
# I won't create a features list
# when in Python, I can just make a list of range(0, len()+1)

# also, I want to remove unwanted cells from the count matrix
# how do I remove entries from a sparse matrix?

# perhaps I can use signac, since its data structure is easy to manipulate
# but first try manipulating a sparse matrix
# see: https://stackoverflow.com/questions/74854550/efficient-way-to-remove-sparse-matrix-row-having-large-data-in-r

str(feature_mtx@assays@data$FeatureMatrix)
# rows are genomic intervals
# cols are cell barcodes
# want to filter by cols to select only cells passed QC


# load metadata table
atac <- read.csv('out/atac/atac_metadata_unfiltered.csv', header=TRUE)
atac[1:2,]

identical(colnames(feature_mtx@assays@data$FeatureMatrix), atac$barcode)
# FALSE, ?
length(colnames(feature_mtx@assays@data$FeatureMatrix))
# 10846
length(atac$barcode)
# 12162, this is before doublet removal

# then I can't use metadata's indices on feature matrix
colnames(feature_mtx@assays@data$FeatureMatrix)[1:2]
# ok, at least feature matrix barcodes do not have prefixes
# I can remove by barcodes, instead of indices

# apply doublet==0, FRIP>=0.5
atac[which(atac$doublet==0 & atac$FRIP>=0.5),][1:2,]
# looks right
nrow(atac[which(atac$doublet==0 & atac$FRIP>=0.5),])
# 9939, correct

atac[which(atac$doublet==0 & atac$FRIP>=0.5),1][1:2]
# view first two barcodes

# create a filtered metadata
filtered_atac <- atac[which(atac$doublet==0 & atac$FRIP>=0.5),]
dim(filtered_atac)
# 9939 8
write.csv(filtered_atac, 'out/atac/atac_metadata.csv', row.names=FALSE)

filtered_barcodes <- atac[which(atac$doublet==0 & atac$FRIP>=0.5),1]
length(filtered_barcodes)
# 9939

filtered_mtx <- feature_mtx@assays@data$FeatureMatrix[,which(colnames(feature_mtx@assays@data$FeatureMatrix) %in% filtered_barcodes)]

dim(feature_mtx@assays@data$FeatureMatrix)
# 164675  10846
dim(filtered_mtx)
# 164675   9939
str(filtered_mtx)
# looks right

# then apply the strategy to extract the GeneScoreMatrix
# because it turned out in python that FeatureMatrix wasn't what we wanted

rm(feature_mtx)
rm(filtered_mtx)
rm(atac)

exp_mtx <- getMatrixFromProject(
  ArchRProj = proj, 
  useMatrix = 'GeneScoreMatrix'
)

# look for the wanted sparse matrix
str(exp_mtx)
# this is the wanted matrix
dim(exp_mtx@assays@data@listData$GeneScoreMatrix)
# 24919 10846
# extract the wanted matrix
mtx <- exp_mtx@assays@data@listData$GeneScoreMatrix

rm(proj)

# follow what I did for FeatureMatrix

# remove barcode prefix
colnames(mtx)[1:2]
colnames(mtx) <- gsub("^.*?#", "", colnames(mtx))

# use filtered_barcodes to slice GeneScoreMatrix
filtered_mtx <- mtx[,which(colnames(mtx) %in% filtered_barcodes)]

# check if I have all needed info from this matrix
str(filtered_mtx)
# No!! it has the right cell barcodes, but for feature names: it is NULL
# where can I find the gene names?
# see: https://github.com/GreenleafLab/ArchR/issues/374#issuecomment-718335424

# rowData(mtx)$name
# Error: unable to find an inherited method for function ‘rowData’ for signature ‘"dgCMatrix"’
# there's no rowData in this extracted sparse matrix

# look into the original GeneScoreMatrix extracted from the ArchR proj obj
rowData(exp_mtx)$name[1:100]
# these are the gene names!
length(rowData(exp_mtx)$name)
# 24919, this is also correct
dim(mtx)
# 24919 10846
dim(filtered_mtx)
# 24919  9939

# add gene names to filtered_mtx
rownames(filtered_mtx) <- rowData(exp_mtx)$name

# does the count mtx actually have values?
rowSums(mtx)[1:100]
# only the first three rowSums are zero
rowSums(exp_mtx@assays@data@listData$GeneScoreMatrix)[1:100]
# confirmed with original GSMatrix, so there wasn't accidental modification on my part
# what are the first three rows anyway?
rowData(exp_mtx)$name[1:4]
# "OR4F5" "LOC729737" "LOC101928626" "FAM87B"      
# how many of the rowSums are zero?
length(which(rowSums(exp_mtx@assays@data@listData$GeneScoreMatrix)==0))
# 159
length(which(rowSums(mtx)==0))
# 159
identical(mtx@x,exp_mtx@assays@data@listData$GeneScoreMatrix@x)
# TRUE
# what about the mtx after filtering by unwanted cell barcodes?
length(which(rowSums(filtered_mtx)==0))
# 162, removing unwanted cells increased the number of genes with zero scores
# sort of expected

# okay, since the filtered count matrix has values, start extracting
writeMM(filtered_mtx, 'out/atac/atac_matrix.mtx')

# write a list of barcodes
write.table(filtered_mtx@Dimnames[2][[1]], file = 'out/atac/atac_barcodes.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)
# write a list of gene names
write.table(filtered_mtx@Dimnames[1][[1]], file = 'out/atac/atac_features.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)


##
# we want to check how marker genes visualize in umap clusters
# see ArchR handbook chapter 7.4: https://www.archrproject.com/bookdown/visualizing-marker-genes-on-an-embedding.html

markerGenes <- c(
  "CD34",  #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "MME", #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A"#TCells
)

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = 'GeneScoreMatrix', 
  name = markerGenes, 
  embedding = 'UMAP', 
  quantCut = c(0.01, 0.95), 
  imputeWeights = NULL
)

# view CD14
p$CD14

# note: right now the umap clusters contain unwanted cells
dim(proj@cellColData)
# 10846 20
# 10846 means these cells don't include doublets (before doublet removel: 12162)
# after filtering by FRiP, we should only retain 9939 cells
colnames(proj@cellColData)
summary(proj@cellColData$PassQC)
# so all these are 1
# maybe I can change some to 0
# but since this is manipulating metadata anyway, why not directly subset by FRiP?

# see: https://www.archrproject.com/bookdown/manipulating-an-archrproject.html

# subset ArchR proj by metadata col
filtered_proj <- proj[proj$FRIP >= 0.5,]

saveArchRProject(
  ArchRProj = filtered_proj, 
  outputDirectory = paste0(getOutputDirectory(proj), '_filtered'), 
  overwrite = TRUE, 
  load = TRUE, 
  dropCells = FALSE, 
  logFile = createLogFile('saveArchRProject'), 
  threads = getArchRThreads()
)

filtered_p <- plotEmbedding(
  ArchRProj = filtered_proj, 
  colorBy = 'GeneScoreMatrix', 
  name = markerGenes, 
  embedding = 'UMAP', 
  quantCut = c(0.01, 0.95), 
  imputeWeights = NULL
)

filtered_p$CD14
# the umap embedding has not changed

# viz the markers on the clusters from unfiltered cells
umap <- plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'Clusters', embedding = 'UMAP')
umap_azimuth1 <- plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'predicted.celltype.l1', embedding = 'UMAP')
umap_azimuth2 <- plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'predicted.celltype.l2', embedding = 'UMAP')
plotPDF(umap, umap_azimuth, umap_azimuth1, name = 'filtered_pbmc10x_umap_clusters.pdf', ArchRProj = filtered_proj, addDOC = FALSE, width = 5, height = 5)

# regenerate umap embedding with filtered cells

filtered_proj <- addIterativeLSI(force = TRUE, ArchRProj = filtered_proj, useMatrix = 'TileMatrix', name = 'IterativeLSI')
# read more about iterative LSI: https://www.archrproject.com/bookdown/iterative-latent-semantic-indexing-lsi.html

# call clusters in this reduced dimension sub-space
filtered_proj <- addClusters(force = TRUE, input = filtered_proj, reducedDims = 'IterativeLSI')

# visualize data using 2D representation such as UMAP
# add a UMAP embedding to ArchRProject
filtered_proj <- addUMAP(force = TRUE, ArchRProj = filtered_proj, reducedDims = 'IterativeLSI')

# view umap
filtered_umap <- plotEmbedding(ArchRProj = filtered_proj, colorBy = 'cellColData', name = 'Clusters', embedding = 'UMAP')
filtered_umap_azimuth1 <- plotEmbedding(ArchRProj = filtered_proj, colorBy = 'cellColData', name = 'predicted.celltype.l1', embedding = 'UMAP')
filtered_umap_azimuth2 <- plotEmbedding(ArchRProj = filtered_proj, colorBy = 'cellColData', name = 'predicted.celltype.l2', embedding = 'UMAP')

# save umap
plotPDF(filtered_umap, filtered_umap_azimuth1, filtered_umap_azimuth2, name = 'filtered_pbmc10x_umap_clusters.pdf', ArchRProj = filtered_proj, addDOC = FALSE, width = 5, height = 5)

# viz markers on clusters from filtered cells
filtered_p <- plotEmbedding(
  ArchRProj = filtered_proj, 
  colorBy = 'GeneScoreMatrix', 
  name = markerGenes, 
  embedding = 'UMAP', 
  quantCut = c(0.01, 0.95), 
  imputeWeights = NULL
)

filtered_p2 <- lapply(filtered_p, function(x) {
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

filtered_p$CD14

do.call(cowplot::plot_grid, c(list(ncol = 3), filtered_p2))

plotPDF(plotList = filtered_p2,
        name = 'plot_filtered_umap_marker_genes.pdf',
        ArchRProj = filtered_proj, 
        addDOC = FALSE, # this adds dates to the plot
        width = 5, 
        height = 5)



