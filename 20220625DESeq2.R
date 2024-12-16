
# 20220625
# functions to produce DESeq results and common up/down reg genes

library(DESeq2)
library(ggplot2)
library(dplyr)

setwd('Labs/xinhua/bulk/')


get_raw <- function (filename, group1, group2) {
	# import raw data
	raw <- read.csv(paste0('data/', filename), header=T, sep=',')
	print(dim(raw))

	# simplify colnames
	colnames(raw)[5:(5+group1+group2-1)] <- sub('_count', '', colnames(raw)[5:(5+group1+group2-1)])
	
	# use unique() to eliminate duplicate rows
	print(nrow(raw))
	raw <- unique(raw)
	print(nrow(raw))
	
	print('get_raw: finished!')
	return(raw)

}

get_res <- function (raw, group1, group2) {

	# create var data
	data <- raw[-c(2,3,4)]
	print(colnames(data))

	# create var metadata
	meta <- data.frame(name=factor(colnames(raw)[5:(5+group1+group2-1)]), group=factor(c(rep('group1',group1), rep('group2',group2))))

	print(meta)

	# construct DESeqDataSet obj
	dds <- DESeqDataSetFromMatrix(countData=data, 
							colData=meta, 
							design=~ group, 
							tidy=T)

	print(dds)

	# run DESeq
	dds <- DESeq(dds)

	# generate results
	res <- results(dds, contrast=c('group', 'group1', 'group2'))

	# check results
	print(head(results(dds, tidy=T)))
	
	print('get_res: finished!')
	return(res)
	
}


get_processed <- function (raw, res) {

	# check if gene order is preserved
	identical(raw[,1], rownames(res))

	# append DESeq results to raw, so I have the gene names again
	processed <- cbind(raw[1:4], res[2:6])
	print(head(processed, 4))
	
	print('get_processed: finished!')
	return(processed)
	
}


write_res <- function (up, down, out_path, z) {
	print(z)
	# write to files: all up/down reg genes
	#write.csv(up, paste0(out_path, '/group1_up_2fc.csv'))
	#write.csv(down, paste0(out_path, '/group1_down_2fc.csv'))
	
	print('write_res: finished!')
}


get_up_down <- function (processed, z, out_path) {
	print(paste0('For group', i, ' ...'))
	
	# select log2FC > 1, which means raw_FC > 2
	# up: group1 up; down: group1 down
	up <- processed[which(res$log2FoldChange >= log(z, base=2)),]
	print(paste0('number of up regulated genes in group1: ', nrow(up)))
	# rank selected genes by padj
	up <- up[order(up$padj, decreasing=F),]
	
	print(nrow(up))
	print(nrow(unique(up)))
	
	# find intersect of prev common up reg genes and current up reg genes
	if (length(all_up)==0) {
		all_up <<- up$gene_name
		print('prev common up reg genes = 0----------')
		}
	all_up <<- intersect(all_up, up$gene_name)
	print(paste0('number of common up reg genes in prev and current samples: ', length(all_up)))

	down <- processed[which(res$log2FoldChange <= log(1/z, base-2)),]
	print(paste0('number of down regulated genes in group1: ', nrow(down)))
	down <- down[order(down$padj, decreasing=F),]
	
	#assign(paste0('sample', i, '_down_FC', z), down, inherits=T)
	print('get_down: finished!')
	
	# write res to csv
	#write_res(up, down, out_path, z)
	return(all_up)
		
	print('get_up: finished!')
	
}

#group <- data.frame(sample1=c(3,7), sample2=c(7,2), sample3=c(4,4), sample4=c(4,4), sample5=c(3,7))

group <- rbind(c(3,7), c(7,2), c(4,4), c(4,4), c(3,7))

#nrow(group)
for (i in 1:3) {
	filename <- list.files(path='data', pattern='_raw.csv$', full.names=F)[i]
	print(filename)
	
	out_path <- sub(pattern='_raw.csv', replacement='', x=filename)
	#dir.create(out_path)
	
	
	# format raw data
	raw <- get_raw(filename, group[i,1], group[i,2])
	
	# run DESeq2
	res <- get_res(raw, group[i,1], group[i,2])
		
	# add res to raw data, so genes have names
	processed <- get_processed(raw, res)
	
	# declare lists
	all_up <- list()
	all_down <- list()
	
	# select log2FC > log(1.5, base=2), which means raw_FC > 1.5
	get_up_down(processed, 2, out_path)
	#get_up_down(processed, 1.5, out_path)
	
}
















######################
# assign(paste0('sample', i, '_up_FC', z), up, inherits=T)
# my first method of obtaining common up/down reg genes from the loop directly was not making progress

	# compare with the prev group of up
	if (length(all_up)==0) {all_up <<- up$gene_name}
	all_up <<- intersect(all_up, up$gene_name)
	print(paste0('number of intersected up reg genes in samples: ', length(all_up)))
	

	# compare with the prev group of down
	if (length(all_down)==0) {all_down <<- down$gene_name}
	all_down <<- intersect(all_down, down$gene_name)
	print(paste0('number of intersected down reg genes in samples: ', length(all_down)))
	
	

	# write to files: intersections of up/down reg gene names from all samples
	cat(all_up, file=paste0('all_up_fc', z, '.csv'), sep='\n')
	cat(all_down, file=paste0('all_down_fc', z, '.csv'), sep='\n')
######################

# second method
# read in up/down reg genes

group1_up_2fc <- list()
common_up_2fc <- data.frame()

for (i in 1:5) {
	filename <- list.files(path='output', pattern='_up_2fc.csv$', full.names=T, recursive=T)[i]
	csv <- read.csv(file=filename, header=T, sep=',')
	# group1_up_2fc <- append(x=group1_up_2fc, values=csv$gene_name): this does not work
	group1_up_2fc[[i]] <- csv$gene_name
	print(paste0('num of genes for sample ', i, ': ', length(group1_up_2fc[[i]])))
	# start 
	if (i > 1) {
		common_up_2fc <- intersect(group1_up_2fc[[i]], group1_up_2fc[[i-1]])
		print(paste0('num of common genes for sample ', i, ' and prior: ', length(common_up_2fc)))
	}
}

write.csv(common_up_2fc, file='common_group1_up_2fc.csv')


# write this into a function

get_common_genes <- function (direction, fc) {
	sample_genes <- list()
	common_genes <- data.frame()
	for (i in 1:5) {
		filename <- list.files(path='output', pattern=paste0('_', direction, '_', fc, 'fc.csv$'), full.names=T, recursive=T)[i]
		csv <- read.csv(file=filename, header=T, sep=',')
		sample_genes[[i]] <- csv$gene_name
		print(paste0('num of genes for sample ', i, ': ', length(sample_genes[[i]])))
		# start 
		if (i == 2) {
			common_genes <- intersect(sample_genes[[i]], sample_genes[[i-1]])
			print(paste0('num of common genes for sample ', i, ' and prior: ', length(common_genes)))
		}
		if (i > 2) {
			common_genes <- intersect(common_genes, sample_genes[[i]])
			print(paste0('num of common genes for sample ', i, ' and prior: ', length(common_genes)))
		}
	}
	write.csv(common_genes, file=paste0('common_group1_', direction, '_', fc, 'fc.csv'))
}

get_common_genes('up', 1.5)
get_common_genes('down', 1.5)
get_common_genes('up', 2)
get_common_genes('down', 2)

###################################
# plot counts with plotCounts
par(mfrow=c(1,3))
plotCounts(dds, gene=data[1,1], intgroup='group')# Pax7
plotCounts(dds, gene=up[1,1], intgroup='group')
plotCounts(dds, gene=down[1,1], intgroup='group')

# volcano plot
par(mfrow=c(1,1))
# make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main='Volcano plot', xlim=c(-3,3)))

# color points
# blue if padj<0.01, red if log2FC>1 and padj<0.05
with(subset(res, padj<0.01), points(log2FoldChange, -log10(pvalue), pch=20, col='blue'))

with(subset(res, padj<0.01 & abs(log2FoldChange)>(log(1.5, base=2))), points(log2FoldChange, -log10(pvalue), pch=20, col='yellow'))

with(subset(res, padj<0.01 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col='red'))


# PCA plot
# transform the raw count data
# vst function will perform variance stabilizing transformation
vsdata <- vst(dds, blind=F)
plotPCA(vsdata, intgroup='group')


# pathway analysis

# KEGG pathways

# following: https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf
# mainly following: https://www.r-bloggers.com/2015/12/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/

BiocManager::install(c('pathview','gage','gageData'))
browseVignettes('pathview')

library(pathview)
# Error: ... there is no package called ‘org.Hs.eg.db’

BiocManager::install('org.Hs.eg.db')
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pathview)
library(gage)
# there is no package called ‘GO.db’

BiocManager::install('GO.db')
library(GO.db)
library(gage)
library(gageData)
# usage of gageData: https://www.bioconductor.org/packages/devel/data/experiment/manuals/gageData/man/gageData.pdf

res$symbol <- mapIds(org.Hs.eg.db,
					keys=row.names(res),
					column='SYMBOL',
					keytype='ENSEMBL',
					multiVals='first')
# 'select()' returned 1:many mapping between keys and columns

res$entrez <- mapIds(org.Hs.eg.db,
					keys=row.names(res),
					column='ENTREZID',
					keytype='ENSEMBL',
					multiVals='first')

res$name <- mapIds(org.Hs.eg.db,
					keys=row.names(res),
					column='GENENAME',
					keytype='ENSEMBL',
					multiVals='first')

head(res, 10)

# use gage (generally applicable gene-set enrichment for pathway analysis)
# to get a list of enriched pathways

# the gageData package has pre-compiled databases mapping genes to KEGG pathways and GO terms for common organisms

# kegg.sets.hs is a named list of 229 elements
# each element is a character vector of member gene Engrez IDs for a single KEGG pathway

data(kegg.sets.hs)
# Warning message: data set ‘kegg.sets.hs’ not found
# see: https://www.rdocumentation.org/packages/gage/versions/2.22.0/topics/kegg.gsets
# generate up-to-date KEGG pathway gene sets
# kegg.sets.hs <- kegg.gsets(species='hsa', id.type='kegg')
# only for 4 species, human, mouse, rat, yeast as well as KEGG Ortholog

data(kegg.sets.hs)
# still same warning; I'll just let it be now

data(sigmet.idx.hs)
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
# KEGG pathway include other types of pathway definitions, like "Global Map" and "Human Diseases", which may be undesirable in pathway analysis
# 'sigmet.idx.hs' is an index of only signaling and metabolic pathways

#head(kegg.sets.hs, 3)

# create a vector of fold changes, where the names of the values are the Entrez gene IDs
# foldchanges <- res$log2FoldChange
# names(foldchanges) <- res$entrez
# head(foldchanges)

# customized method:

# get rid of foldchanges 1/2 < rawFC < 2, or -1 < logFC < 1
colnames(res)
select_lfc <- data.frame(lfc=res$log2FoldChange, entrez=res$entrez)
typeof(select_lfc)
nrow(select_lfc)
summary(select_lfc)
select_lfc <- subset(select_lfc, lfc>=1|lfc<=(-1))
nrow(select_lfc)# 5567
summary(select_lfc)

foldchanges <- select_lfc$lfc
names(foldchanges) <- select_lfc$entrez


# same.dir=F: test for changes towards both directions simultaneously; in KEGG, BioCarta pathways, genes frequently are not co-regulated, hence same.dir=F
# same.dir=T: test for changes in a single direction only; all genes must be up or down regulated: for experimentally derived gene sets, GO term group, etc., co-regulation is commonly the case, hence same.dir=T
keggres <- gage(exprs=foldchanges, gsets=kegg.sets.hs, same.dir=T)

lapply(keggres, head)

# write tables
write.table(keggres$greater, file='kegg/lfc2/res.kegg.greater.tsv', sep='\t')
write.table(keggres$less, file='kegg/lfc2/res.kegg.less.tsv', sep='\t')
write.table(keggres$stats, file='kegg/lfc2/res.kegg.stats.tsv', sep='\t')


# process the results to pull out the top 5 upregulated pathways
# keggres_pathways <- data.frame(id=rownames(keggres$greater), keggres$greater) %>% tbl_df() %>% filter(row_number()<=5) %>% .$id %>% as.character()
# Warning: `tbl_df()` was deprecated in dplyr 1.0.0.
# Please use `tibble::as_tibble()` instead.

keggres_pathways <- data.frame(id=rownames(keggres$greater), keggres$greater) %>% tibble::as_tibble() %>% filter(row_number()<=5) %>% .$id %>% as.character()

keggres_pathways

# get just the IDs
keggresids <- substr(keggres_pathways, start=1, stop=8)
keggresids


# visualize
# pid???
pathview(gene.data <- foldchanges, pathway.id='pid', species='hsa', new.signature=F)

# define plotting function for applying later
plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species='hsa', new.signature=F)

getwd()
# setwd('Labs/xinhua/bulk/')

# plot multiple pathways
# download and save plots to disk, and return a throwaway list obj
tmp <- sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species='hsa'))

# view sig genes, default cutoff=0.1
keggres_sig <- sigGeneSet(keggres, outname='res.kegg')
str(keggres_sig)
# only 2 sig up-reg gene sets and 3 sig down-reg gene sets
# heatmap was not correctly produced

# heatmap(keggres$greater)
# cannot produce correct heatmaps, no sample info

write.table(keggres_sig$greater, file='res.kegg.sig.greater.tsv', sep='\t')
write.table(keggres_sig$less, file='res.kegg.sig.less.tsv', sep='\t')
write.table(keggres_sig$stats, file='res.kegg.sig.stats.tsv', sep='\t')



# Gene Ontology
# BP: biological process
# CC: cellular component
# MF: molecular function

data(go.sets.hs)
data(go.subs.hs)

gobpsets <- go.sets.hs[go.subs.hs$BP]
gobpres <- gage(foldchanges, gsets=gobpsets, same.dir=T)

lapply(gobpres, head)

# tutorial ends here
str(gobpres)
write.table(gobpres$greater, file='res.gobp.greater.tsv', sep='\t')
write.table(gobpres$less, file='res.gobp.less.tsv', sep='\t')
write.table(gobpres$stats, file='res.gobp.stats.tsv', sep='\t')

gobpres_sig <- sigGeneSet(gobpres, outname='res.gobp')
str(gobpres_sig)
write.table(gobpres_sig$greater, file='res.gobp.sig.greater.tsv', sep='\t')
write.table(gobpres_sig$less, file='res.gobp.sig.less.tsv', sep='\t')
write.table(gobpres_sig$stats, file='res.gobp.sig.stats.tsv', sep='\t')

goccsets <- go.sets.hs[go.subs.hs$CC]
goccres <- gage(foldchanges, gsets=goccsets, same.dir=T)
str(goccres)
write.table(goccres$greater, file='res.gocc.greater.tsv', sep='\t')
write.table(goccres$less, file='res.gocc.less.tsv', sep='\t')
write.table(goccres$stats, file='res.gocc.stats.tsv', sep='\t')

goccres_sig <- sigGeneSet(goccres, outname='res.gocc')
str(goccres_sig)
write.table(gobpres_sig$greater, file='res.gocc.sig.greater.tsv', sep='\t')
write.table(gobpres_sig$less, file='res.gocc.sig.less.tsv', sep='\t')
write.table(gobpres_sig$stats, file='res.gocc.sig.stats.tsv', sep='\t')

gomfsets <- go.sets.hs[go.subs.hs$MF]
gomfres <- gage(foldchanges, gsets=gomfsets, same.dir=T)
str(gomfres)
write.table(gomfres$greater, file='res.gomf.greater.tsv', sep='\t')
write.table(gomfres$less, file='res.gomf.less.tsv', sep='\t')
write.table(gomfres$stats, file='res.gomf.stats.tsv', sep='\t')

gomfres_sig <- sigGeneSet(gomfres, outname='res.gomf')
str(gomfres_sig)
write.table(gomfres_sig$greater, file='res.gomf.sig.greater.tsv', sep='\t')
write.table(gomfres_sig$less, file='res.gomf.sig.less.tsv', sep='\t')
write.table(gomfres_sig$stats, file='res.gomf.sig.stats.tsv', sep='\t')




