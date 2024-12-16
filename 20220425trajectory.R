
# 20220425
# Monocle 2 draft
# vignette: http://cole-trapnell-lab.github.io/monocle-release/docs/
# feature selection: look at interesting genes (not just noisy genes expressed at low values) + supplied markers (or no supply)
# dimension reduction: reversed graph embedding
# fitting the best tree: manifold learning


library(Seurat)
library(dplyr)
library(ggplot2)
# BiocManager::install('monocle')
library(monocle)

setwd('~/Labs/xinhua/scSeq/')
#set.seed(1)

data <- readRDS(file='rds/sample_all_cluster_1to6_integ.rds')

DefaultAssay(object=data) <- 'RNA'

# view seurat identities based current res
length(data@active.ident)
# 20343
data@active.ident[1:5]

# view all seurat identities based on all calculated res
dim(data@meta.data)
# 20343*10
data@meta.data[1:10,]

# assign desired seurat identities to cell_type
data[['ident']] <- data@meta.data[,9]

# trim meta.data
data@meta.data <- subset(data@meta.data, select=c(1,4,9))
data@meta.data <- data@meta.data %>% dplyr::rename(ident=integrated_snn_res.1)
data@meta.data[1:10,]

saveRDS(object=data, file='rds/data.rds')
data <- readRDS(file='rds/data.rds')


### prepare count matrix

# view celluar gene expression
dim(data@assays$RNA@counts)
# 36601*20343
data@assays$RNA@counts[1:10,1:2]
# genes as rownames, cells as colnames
# how many genes have non-zero expression?
length(which(rowSums(data@assays$RNA@counts)!=0))
# 27569 of 36601 have non-zero expression

# test method: select rows by line numbers
data@assays$RNA@counts[c(1,2,3),1:2]

# create expression matrix
exp <- data@assays$RNA@counts[which(rowSums(data@assays$RNA@counts)!=0),]
dim(exp)
# 27569*20343
# genes*cells

# choice of expressionFamily
# acc to: https://www.jianshu.com/p/5d6fd4561bc0

# FPKM/TPM are often normal dist: tobit()
# UMI and transcript counts can be modeled better with negative binomial dist with fixed variance: negbinomial.size()
# negbinomial() also works for UMIs, slightly more accurate, but much much slower; not recommended unless for very small datasets
cds <- newCellDataSet(as(as.matrix(exp),'sparseMatrix'),
					phenoData=new('AnnotatedDataFrame', data=data@meta.data),
					featureData=new('AnnotatedDataFrame', data.frame(gene_short_name=row.names(data), row.names=row.names(data))),
					lowerDetectionLimit=0.5,
					expressionFamily=negbinomial.size()
					)
# invalid class “CellDataSet” object: 1: feature numbers differ between assayData and featureData
# invalid class “CellDataSet” object: 2: featureNames differ between assayData and featureData

# acc to: https://www.jianshu.com/p/5d6fd4561bc0
# ncol(exp) must equal nrow(phenoData): cell number
# nrow(exp) must equal nrow(featureData): gene number

dim(new('AnnotatedDataFrame', data=data@meta.data))
# i.e. the dim of meta.data, 20343*3
# rows are cell names
dim(new('AnnotatedDataFrame', data.frame(gene_short_name=row.names(data), row.names=row.names(data))))
# i.e. nrow(data), 2000
length(row.names(data))# 2000
nrow(data)# 2000
# why is it 2000?
rownames(data)[1:5]
# gene names
# I suspect data now only has variable features
which(VariableFeatures(data)!=row.names(data))
identical(VariableFeatures(data), row.names(data))# TRUE

dim(data)
# 2000*20343
colnames(data)[1:2]
# cell names

row.names(exp)[1:5]
# I will change var data to meta.data
cds <- newCellDataSet(as(as.matrix(exp),'sparseMatrix'),
					phenoData=new('AnnotatedDataFrame', data=data@meta.data),
					featureData=new('AnnotatedDataFrame', data.frame(gene_short_name=row.names(exp), row.names=row.names(exp))),
					lowerDetectionLimit=0.5,
					expressionFamily=negbinomial.size()
					)

saveRDS(cds, file='rds/cds.rds')
# difference between row.names() and rownames()
# row.names() operates faster on data.frame; using rownames() will call row.names() as an intermediate step
# rownames() operates faster on other data objects; using row.names() will call rownames() as an intermediate step

dim(cds)
# features: 27560
# samples: 20343
cds[1:5, 1:2]



### import Seurat object directly as CellDataSet object
# (not tested)
# importCDS(data)


### estimate size factor and dispersion
cds <- estimateSizeFactors(cds)# necessary for calling dispersions
cds <- estimateDispersions(cds)

cds <- readRDS(file='rds/cds_after_dispersion.rds')

### filter genes
# after filtering multiplets with Seurat based on expression values
# filter genes based on number of cells expressing this gene

# set global expression detection threshold for this CellDataSet
# counts how many cells each feature in a CDS obj are detectably expressed above a min threshold
# counts also number of genes above this threshold are detectable in each cell
cds <- detectGenes(cds, min_expr=0.1)

# view number of genes
print(head(fData(cds)))

# for later use, if needed
# filter genes expressed in less than 10 cells
expressed_genes <- subset(fData(cds), num_cells_expressed >= 10) %>% row.names()
length(expressed_genes)
# 21577


### select features

pdf(file='ordering_genes_by_VariableFeatures.pdf', width=20, height=20)

# using variable features detected by Seurat
data <- FindVariableFeatures(data)
var_features <- VariableFeatures(data)
cds <- setOrderingFilter(cds, ordering_genes=var_features)

plot_ordering_genes(cds)

# using clusters DE genes
de_markers <- FindAllMarkers(data) %>% subset(p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds, ordering_genes=de_markers)

plot_ordering_genes(cds)

# using variable features selected by Monocle
disp_genes <- subset(dispersionTable(cds), subset=(mean_expression>=0.1 & dispersion_empirical>=1*dispersion_fit))$gene_id

disp_genes <- subset(dispersionTable(cds), subset=(mean_expression>=0.1))$gene_id

cds <- setOrderingFilter(cds, disp_genes)

plot_ordering_genes(cds)

# using dpFeature, detect important ordering genes from data, not based on biological knowledge
diff <- differentialGeneTest(cds[expressed_genes,], 
							fullModelFormulaStr='~ident', 
							cores=1)
head(diff)

deg <- subset(diff, qval<0.01)
deg <- deg[order(deg$qval, decreasing=F),]
head(deg)

write.table(deg, 
			file='train.monocle.DEG.xls', 
			col.names=T, 
			row.names=F, 
			sep='\t', 
			quote=F)

cds <- setOrderingFilter(cds, rownames(deg))

plot_ordering_genes(cds)

# number of genes is best at ~2000
# use only top genes if necessary
cds <- setOrderingFilter(cds, rownames(deg)[order(deg$qval)][1:400])

plot_ordering_genes(cds)



### dimension reduction
cds <- reduceDimension(cds, max_components=2, method='DDRTree')



### establish pseudotime trajectory & order cells in pseudotime
cds <- orderCells(cds)

# if needed, set the number of roots
# cds <- orderCells(cds, root_state=5)

# visualize
cds <- readRDS(file='rds/cds_after_orderCells.rds')

str(cds)
# see factors that can use for visualize

#print(plot_cell_trajectory(cds, color_by='Size_Factor', cell_size=0.3, show_tree=TRUE))
# size factor is a scaling factor used to divide the raw counts of a particular cell to obtain normalized expression values,thus allowing downstream comparisons between cells that are not affected by differences in library size or total RNA content.

# color by pseudotime
print(plot_cell_trajectory(cds, color_by='Pseudotime', cell_size=0.3, show_tree=TRUE))
print(plot_cell_trajectory(cds, color_by='Pseudotime', cell_size=0.3, show_tree=TRUE)+ facet_wrap('~State', ncol=3))

# color by state
print(plot_cell_trajectory(cds, color_by='State', cell_size=0.3, show_tree=TRUE))
print(plot_cell_trajectory(cds, color_by='State', cell_size=0.3, show_tree=TRUE) + facet_wrap('~State', ncol=3))

# color by samples
print(plot_cell_trajectory(cds, color_by='sample', cell_size=0.3, show_tree=TRUE))

print(plot_cell_trajectory(cds, color_by='sample', cell_size=0.3, show_tree=TRUE) + facet_wrap('~sample', nrow=5))

# color by clusters
print(plot_cell_trajectory(cds, color_by='ident', cell_size=0.3, show_tree=TRUE))

print(plot_cell_trajectory(cds, color_by='ident', cell_size=0.3, show_tree=TRUE) + facet_wrap('~ident', nrow=5))

print(plot_cell_trajectory(cds, color_by='cluster', cell_size=0.3, show_tree=TRUE))

print(plot_cell_trajectory(cds, color_by='cluster', cell_size=0.3, show_tree=TRUE) + facet_wrap('~cluster', nrow=2))


dev.off()



plot_pseudotime_heatmap
Time_diff <- differentialGeneTest(cds[var_features,], 
								cores = 1, 
								fullModelFormulaStr = "~sm.ns(Pseudotime)")
								
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改

write.csv(Time_diff, "Time_diff_all.csv", row.names = F)

Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)

ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)












