
# same-cell single-cell RNA- and ATAC-seq multiomic data integration
# doublet analysis using ArchR


setwd('Labs/wijst/')
rds <- readRDS('QualityControl/pbmc_granulocyte_sorted_10k_atac/pbmc_granulocyte_sorted_10k_atac-Doublet-Summary.rds')

rds
colnames(rds)
str(rds)
rds$doubletResults
rds$doubletResults$doubletUMAP[1,]
rds$doubletResults$doubletScoreUMAP
# I see many cells with doublet score == 0
# I read: "Most doublet detection methods provide a ‘doublet score’ that is higher on average in doublets than in singlets, and users are left to decide on a threshold above which droplets will be excluded as doublets."
# https://doi.org/10.12688%2Ff1000research.73600.2
# so I think only score == 0 are singlets
length(which(rds$doubletResults$doubletScoreUMAP==0))
# 11234
length(rds$doubletResults$doubletScoreUMAP)
# 12162

# filtered 12162-11234 = 928 doublets

# create a list for singlets
# I see that barcode names are stored in attr(obj, 'attribute name')
arch_barcodes <- attr(rds$doubletResults$doubletScoreUMAP,'names')
arch_barcodes[1:5]
# remove the prefix
arch_barcodes <- gsub("^.*?#", "", arch_barcodes)
arch_barcodes[1:5]


# now create list
# of singlets
#archr <- cbind(arch_barcodes[which(rds$doubletResults$doubletScoreUMAP==0)],
#      rep('False',length(which(rds$doubletResults$doubletScoreUMAP==0))))
# of doublets
archr <- cbind(arch_barcodes[which(rds$doubletResults$doubletScoreUMAP!=0)],
      rep('True',length(which(rds$doubletResults$doubletScoreUMAP!=0))))
archr[1:5,]
write.csv(archr, 'out/archr_doublets.csv', row.names=FALSE)

#archr <- read.csv('out/archr_singlets.csv', header=TRUE)
archr <- read.csv('out/archr_doublets.csv', header=TRUE)
archr[1:5,]
nrow(archr)
# 928

# import scrublet predictions
scrub <- read.csv('out/scrublet_predictions.csv', header=FALSE)
scrub[1:5,]
scrub[which(scrub[,2]=='True'),][1:5,]
nrow(scrub[which(scrub[,2]=='False'),])
# number of singlets: 10994

length(arch_barcodes)
# number of cells from atac: 12162
nrow(scrub)
# number of cells from rna: 11898
length(intersect(arch_barcodes, scrub[,1]))
# number of shared cells from both datasets: 11494

# doublets identified by scrublet 
scrub <- scrub[which(scrub[,2]=='True'),]

# compare archr and scrublet predictions
nrow(archr)
# 11234
nrow(scrub)
# 904

# Out of 12162 cells, archr found 11234 singlets, 928 doublets
# Out of 11898 cells, scrub found 10994 singlets, 904 doublets

# for singlets
length(which(archr[,1] %in% scrub[,1]))
# 10376 cells from archr are also in scrub
length(which(scrub[,1] %in% archr[,1]))
# 10376 cells from scrub are also in archr

# for doublets
length(intersect(archr[,1], scrub[,1]))
# 618 doublets are shared
# 10994 = 618 + 10376 (number of shared doublets + number of shared singlets = # of singlets found by scrub, weird coincidence?!)

# are these 10376 or 618 the same cells
all(archr[,1][which(archr[,1] %in% scrub[,1])] %in% scrub[,1][which(scrub[,1] %in% archr[,1])])
# TRUE

# plot scrublet doublets in ArchR obj

# make a table of scrub doublets
scrub <- read.csv('out/scrublet_predictions.csv', header=FALSE)
scrub <- scrub[which(scrub[,2]=='True'),]
scrub[1:5,]

# load archr obj
library(ArchR)
#install.packages('hexbin')
proj <- loadArchRProject(path = 'pbmc10x_archr', force = FALSE, showLogo = TRUE)
proj
# nCells: 10846; where did I lose cells??

# add prefix to scrub barcode names, so names match archr obj metadata names
scrub[,1] <- paste0("pbmc_granulocyte_sorted_10k_atac#", scrub[,1])
scrub[1:5,]

# add a col to archr obj metadata
colnames(proj@cellColData)
proj@cellColData$scrublet_doublet <- scrub[match(rownames(proj@cellColData), scrub[,1]),2]
#proj@cellColData$archr_doublet <- archr[match(rownames(proj@cellColData), archr[,1]),2]

proj@cellColData[1,]

# plot scrublet doublets
plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'scrublet_doublet', embedding = 'UMAP')
plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'archr_doublet', embedding = 'UMAP')
# Error: polygon edge not found
# okay, the scrublet_doublet col contains: NA and 'True', perhaps not plottable

# replace NA with 0, and 'True' with 1
#proj@cellColData$scrublet_doublet[is.na(proj@cellColData$scrublet_doublet)] <- 0
#proj@cellColData$scrublet_doublet[proj@cellColData$scrublet_doublet=='True'] <- 1
proj@cellColData$archr_doublet[proj@cellColData$scrublet_doublet != 'False'] <- 'True'
#proj@cellColData$scrublet_doublet <- as.numeric(proj@cellColData$scrublet_doublet)
#proj@cellColData$scrublet_doublet[1:100]
#plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'scrublet_doublet', embedding = 'UMAP')
# Error: polygon edge not found
#proj@cellColData$DoubletEnrichment[1:100]
#typeof(proj@cellColData$scrublet_doublet)
#typeof(proj@cellColData$DoubletEnrichment)

plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'DoubletScore', embedding = 'UMAP')
plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'DoubletEnrichment', embedding = 'UMAP')
plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'scrublet_doublet', embedding = 'UMAP')


# I realized there's a discrepancy between the number of cells labeled as doublets
# with rds$doubletResults$doubletScoreUMAP!=0, and the number of cells actually 
# thrown out by running filterDoublets()
# Total nCell is 12162, only 10846 are retained in the ArchR proj after running filterDoublets()
# as in message: removing 1316 of 12162 (10.8%)
# 10846+1316=12162
# therefore, it is not 928 cells removed, but 1316
# where did this data come from?
# I re-ran the ArchR protocol, and found out 928 cells are assigned non-zero DoubletScore

# after addDoubletScores() and creating ArchRProject()
length(which(proj@cellColData$DoubletScore!=0))
# 928
length(which(proj@cellColData$DoubletEnrichment==0))
# 45

# however, when removing doublets
proj <- filterDoublets(ArchRProj = proj)
# 1316 of 12162 (10.8%) are removed!
# this means that non-zero DoubletScore is not the only determining factor as to which cells are removed as doublets

length(which(rds$doubletResults$doubletScoreUMAP==0))
# 11234 = 12162 - 928, this is the right place to look for DoubletScore==0

# so I will use the doublet_summary rds to extract all cell barcodes, 
# and the proj obj to get the cell barcodes of the retained singlets

rds <- readRDS('QualityControl/pbmc_granulocyte_sorted_10k_atac/pbmc_granulocyte_sorted_10k_atac-Doublet-Summary.rds')
atac_barcodes <- attr(rds$doubletResults$doubletScoreUMAP,'names')

# create a table of cell barcodes with a col of all True's
archr <- cbind(atac_barcodes, rep('True', length(atac_barcodes)))
archr[1:5,]
nrow(archr)

# turn all retained siglets into False
archr[rownames(proj@cellColData) %in% archr[,1],2] <- 'False'

length(which(archr[,2]=='True'))
length(which(archr[,2]=='False'))
# why all True's turned into False's??

length(rownames(proj@cellColData) %in% archr[,1])
# 10846
length(archr[rownames(proj@cellColData) %in% archr[,1],2])
# 12162, why?

# try sth else
length(archr[match(rownames(proj@cellColData), archr[,1]),2])
# 10846

archr[match(rownames(proj@cellColData), archr[,1]),2] <- 'False'
length(which(archr[,2]=='True'))
# 1316
length(which(archr[,2]=='False'))
# 10846

write.csv(archr, 'out/archr_predictions.csv', row.names=FALSE)






