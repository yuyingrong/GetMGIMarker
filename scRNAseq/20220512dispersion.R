
# dispersion analysis using seurat and monocle

setwd('~/fx/')

#library(Seurat)
library(dplyr)
library(monocle)

#data <- readRDS(file='~/fx/rds/data.rds')

### estimate size factor and dispersion
#cds <- estimateSizeFactors(cds)# necessary for calling dispersions
#cds <- estimateDispersions(cds)

#saveRDS(cds, file='~/fx/rds/cds_after_dispersion.rds')

cds <- readRDS(file='~/fx/rds/cds_after_dispersion.rds')

print('loaded RDS')


### filter genes
# skipped


### select features

## using variable features detected by Seurat
#data <- FindVariableFeatures(data)
#var_features <- VariableFeatures(data)
#cds <- setOrderingFilter(cds, ordering_genes=var_features)

## using variable features selected by Monocle
#disp_genes <- subset(dispersionTable(cds), subset=(mean_expression>=0.1 & dispersion_empirical>=1*dispersion_fit))$gene_id
#cds <- setOrderingFilter(cds, disp_genes)

## using dpFeature to detect important ordering genes from data

cds <- detectGenes(cds, min_expr=0.1)

# view number of genes
print(head(fData(cds)))

# filter genes expressed in less than 10 cells
expressed_genes <- subset(fData(cds), num_cells_expressed >= 10) %>% row.names()
print(length(expressed_genes))

diff <- differentialGeneTest(cds[expressed_genes,], 
			fullModelFormulaStr='~ident', 
			cores=1)

saveRDS(cds, file='~/fx/rds/cds_after_dpFeature.rds')

print(head(diff))

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


###########

print('set ordering filter')

pdf(file='~/fx/ordering_genes_by_dpFeature.pdf', width=20, height=20)

print(plot_ordering_genes(cds))


### dimension reduction
cds <- reduceDimension(cds, max_components=2, method='DDRTree')

print('reduced dimension')

### establish pseudotime trajectory & order cells in pseudotime
cds <- orderCells(cds)

print('ordered cells')


## visualization

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

saveRDS(cds, file='~/fx/rds/cds_after_orderCells_by_dpFeature.rds')

