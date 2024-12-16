
# 20220628
# functions to run GO/KEGG analysis 

library(DESeq2)
library(ggplot2)
library(dplyr)

library(AnnotationDbi)
library(org.Hs.eg.db)
library(pathview)
library(GO.db)
library(gage)
library(gageData)

data(kegg.sets.hs)
data(sigmet.idx.hs)

data(go.sets.hs)
data(go.subs.hs)


kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]


setwd('Labs/xinhua/bulk/')


# 
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

	# construct DESeqDataSet obj
	dds <- DESeqDataSetFromMatrix(countData=data, 
							colData=meta, 
							design=~ group, 
							tidy=T)
	print('get_res: dds created!')

	# run DESeq
	dds <- DESeq(dds)
	print('get_res: DESeq run!')

	# generate results
	res <- results(dds, contrast=c('group', 'group1', 'group2'))
	
	# check results
	print(head(results(dds, tidy=T)))
	print('get_res: done!')
	return(res)	
}


map_id <- function (res) {
	# add more columns
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
	
	print('map_id: done!')
	return(res)
}

get_foldchanges <- function (res, z) {
	# create intermediate df
	select_lfc <- data.frame(lfc=res$log2FoldChange, entrez=res$entrez)
	select_lfc <- subset(select_lfc, lfc>=log(z, base=2) | lfc<=log(1/z, base=2))
	
	# create foldchanges
	foldchanges <- select_lfc$lfc
	names(foldchanges) <- select_lfc$entrez
	
	print('get_foldchanges: done!')
	return(foldchanges)
}


write_kegg <- function (foldchanges, z, out_path) {
	keggres <- gage(exprs=foldchanges, gsets=kegg.sets.hs, same.dir=T)
	print(lapply(keggres, head))
	
	# create output dir
	dir.create(paste0('output/', out_path, '/kegg'))
	
	# write tables
	write.table(keggres$greater, file=paste0('output/', out_path, '/kegg/lfc', z, '_res.kegg.greater.tsv'), sep='\t')
	write.table(keggres$less, file=paste0('output/', out_path, '/kegg/lfc', z, '_res.kegg.less.tsv'), sep='\t')
	write.table(keggres$stats, file=paste0('output/', out_path, '/kegg/lfc', z, '_res.kegg.stats.tsv'), sep='\t')
	
	print('write_kegg: getting sig genes...')
	keggres_sig <- sigGeneSet(keggres, outname='res.kegg')
	write.table(keggres_sig$greater, file=paste0('output/', out_path, '/kegg/lfc', z, '_res.kegg.sig.greater.tsv'), sep='\t')
	write.table(keggres_sig$less, file=paste0('output/', out_path, '/kegg/lfc', z, '_res.kegg.sig.less.tsv'), sep='\t')
	write.table(keggres_sig$stats, file=paste0('output/', out_path, '/kegg/lfc', z, '_res.kegg.sig.stats.tsv'), sep='\t')
	
	print('write_kegg: done!')
	print('returned keggres!')
	return(keggres)
}


plot_pathway <- function (keggres, foldchanges) {
	# pull out the top 5 upregulated pathways
	keggres_pathways <- data.frame(id=rownames(keggres$greater), keggres$greater) %>% tibble::as_tibble() %>% filter(row_number()<=5) %>% .$id %>% as.character()
	
	# get just the IDs
	keggresids <- substr(keggres_pathways, start=1, stop=8)
	
	# visualize
	pathview(gene.data <- foldchanges, pathway.id=pid, species='hsa', new.signature=F)
	
	# define plotting function for applying later
	plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species='hsa', new.signature=F)
	
	# plot, download, and save pathways
	tmp <- sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species='hsa'))
	
	print('plot_pathway: done!')	
}


write_go <- function (pathway, foldchanges, z, out_path) {
	# create output dir
	dir.create(paste0('output/', out_path, '/go'))
	
	# get BP/CC/MF pathways
	gosets <- go.sets.hs[go.subs.hs[[pathway]]]
	gores <- gage(foldchanges, gsets=gosets, same.dir=T)
	
	# write BP/CC/MF pathways
	write.table(gores$greater, file=paste0('output/', out_path, '/go/', pathway, '_lfc', z, '_res.go', pathway, '.greater.tsv'), sep='\t')
	write.table(gores$less, file=paste0('output/', out_path, '/go/', pathway, '_lfc', z, '_res.go', pathway, '.less.tsv'), sep='\t')
	write.table(gores$stats, file=paste0('output/', out_path, '/go/', pathway, '_lfc', z, '_res.go', pathway, '.stats.tsv'), sep='\t')
	
	# write sig BP/CC/MF pathways
	gores_sig <- sigGeneSet(gores, outname=paste0('res.go', pathway))
	write.table(gores_sig$greater, file=paste0('output/', out_path, '/go/', pathway, '_lfc', z, '_res.go', pathway, '.sig.greater.tsv'), sep='\t')
	write.table(gores_sig$less, file=paste0('output/', out_path, '/go/', pathway, '_lfc', z, '_res.go', pathway, '.sig.less.tsv'), sep='\t')
	write.table(gores_sig$stats, file=paste0('output/', out_path, '/go/', pathway, '_lfc', z, '_res.go', pathway, '.sig.stats.tsv'), sep='\t')
	
	print('write_go: done!')
}


write_results <- function(z, res, foldchanges, out_path) {
	# get only the fold change < 1/z and > z genes
	foldchanges  <- get_foldchanges(res, z)
	
	# kegg analysis
	keggres <- write_kegg(foldchanges, z, out_path)
	#plot_pathway(keggres, foldchanges)
	
	# go analysis
	for (pathway in c('BP', 'CC', 'MF')) {
		write_go(pathway, foldchanges, z, out_path)
	}
}


group <- rbind(c(3,7), c(7,2), c(4,4), c(4,4), c(3,7))

for (i in 3:5) {
	filename <- list.files(path='data', pattern='_raw.csv$', full.names=F)[i]
	print(filename)
	
	out_path <- sub(pattern='_raw.csv', replacement='', x=filename)
	#dir.create(paste0('output/', out_path))
	
	# format raw data
	raw <- get_raw(filename, group[i,1], group[i,2])
	
	# run DESeq2
	res <- get_res(raw, group[i,1], group[i,2])
	
	# map id
	res <- map_id(res)
	
	# enter interested foldchange cutoffs
	write_results(2, res, foldchanges, out_path)
	write_results(1.5, res, foldchanges, out_path)
}








