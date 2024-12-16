
# pseudo-bulk normalization
# using DESeq2
# method protocol: https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html


library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

setwd('~/Labs/xinhua/scSeq/')

merged <- readRDS(file='rds/merged.rds')

# get rid of unassigned cells
merged <- subset(merged, cells=rownames(merged@meta.data)[-which(merged@meta.data$cluster=='unassigned')])

# view clusters
VlnPlot(merged, features = c('nFeature_RNA','nCount_RNA'),ncol=2)
VlnPlot(merged, features = c('percent_mt','log10GenesPerUMI'),ncol=2)


# get cell:gene expression matrix
counts <- GetAssayData(merged,slot='data')
dim(counts)
# 36601 * 25320

# get rid of zero expression genes from cell:gene expression matrix
length(which(rowSums(counts)==0))
# 8750 genes have zero expression by any cell

counts <- counts[-which(rowSums(counts)==0),]
dim(counts)
# 27851 * 25320

# get metadata:cell matrix
metadata <- merged@meta.data[,1:4]


# normalize for all clusters and samples
for (clust in c(paste0('clust', 1:7))) {
	print(clust)
	x <- matrix(0L, nrow=nrow(counts), ncol=6)
	rownames(x) <- rownames(counts)
	colnames(x) <- c(paste0('fx',seq(1,6)))
	for (i in 1:6) {
		fx <- c(paste0('fx', i))
		x[,i] <- counts[,which(metadata$cluster==clust & metadata$sample==fx)] %>% rowSums()
		#x[,i] <- x[,i]/sum(colSums(counts[,which(metadata$cluster==clust & metadata$sample==fx)]))
	}
	assign(clust, x)
}

# save the result in a list
clusters <- list(clust1, clust2, clust3, clust4, clust5, clust6, clust7)

# view normalization results
for (i in 1:length(clusters)) {
	print(paste0('>>> cluster ',i, 'Pax7, Myf5, MyoD1'))
	print(clusters[[i]]['PAX7',])
	print(clusters[[i]]['MYF5',])
	print(clusters[[i]]['MYOD1',])
}

# remove vars from env
# rm(clust1, clust2, clust3, clust4, clust5, clust6, clust7)

str(clusters)

dim(clusters[[1]])

saveRDS(clusters, file='rds/clusters_normalized.rds')

# might update stuff above to work in a list directly one day



# install and load DESeq2
BiocManager::install('DESeq2')
library(DESeq2)


# try one pair: cluster 1 vs 4
# tutorials used
# bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html
# jianshu.com/p/b5541d695108


# create DESeq obj

# identical(colnames(counts), rownames(metadata))
# metadata[1,]
# actually I won't be working single cell counts and metadata

# I will work with the following input:
# data: obtain from clusters: rows are genes, cols are clust1_fx1, fx2, ... clust2_fx1, clust2_fx2, ...
# meta: use factor(): clusters and samples as colnames

# create matrix data
data <- cbind(clusters[[1]], clusters[[4]])

colnames(data)
# only sample names; cluster names need to be added
colnames(data) <- c(paste0('clust1_', colnames(data)[1:6]), 
					paste0('clust4_', colnames(data)[7:12]))

dim(data)
data[1,]

# create data.frame meta
meta <- data.frame(cluster=factor(rep(c('clust1', 'clust4'), each=6), 
							levels=c('clust1', 'clust4')),
					sample=factor(paste0('fx',1:6),
							levels=paste0('fx',1:6)))

meta

# create DESeq object
# counts must be integers
clust_1v4 <- DESeqDataSetFromMatrix(countData=data, 
									colData=meta, 
									design= ~ cluster+sample)

# normalize the DESeq obj
# turn on parallel when you have a large dataset
clust_1v4 <- DESeq(clust_1v4, parallel=TRUE)

# res <- results(clust_1v4)
# default contrasts sample fx6 vs fx1
# I want to contrast cluster 1 vs 4
res <- results(clust_1v4, contrast=c('cluster','clust1', 'clust4'))
res[1:10,]
# more about results() here: https://rdrr.io/bioc/DESeq2/man/results.html

# what if I contrast cluster 4 vs 1?
results(clust_1v4, contrast=c('cluster','clust4', 'clust1'))[1:10,]
# basically, signs before log2FoldChange and stat are flipped
# why?

res['PAX7',]
res['MYF5',]
res['MYOD1',]

# instead of using results(x, contrast=...), use name=... for continuous variables, indiv effects or indiv interaction terms
# the value provided to name must be an element of:
resultsNames(clust_1v4)


plotMA(clust_1v4, contrast=c('cluster','clust1', 'clust4'))
plotMA(res)# the same as above
summary(res)

# how many genes have fold_change >2, that is, log_fold_change >1?
length(which(res$log2FoldChange >1))
# 6440

# get all 6440 genes
out <- res[which(res$log2FoldChange >1),]
dim(out)
# 6440*6
class(out)
as.matrix(out)[1,]
out <- as.data.frame(out)
out[1:2,]

# order the 6440 genes by log2_fold_change, desc
out <- out[order(out$log2FoldChange, decreasing=TRUE),]
out[1:5,]

# write the 6440 into table in TXT
write.table(out, file='clust_1v4_normalized.txt', sep='\t', quote=F, col.names=NA)

# check the counts slot of normalized output
counts(clust_1v4,normalize=TRUE)[1,]
counts(clust_1v4,normalize=FALSE)[1,]

# add normalized expression values
resdata <- merge(as.data.frame(res), 
		as.data.frame(counts(clust_1v4,normalize=TRUE)), 
		by='row.names', 
		sort=FALSE)

resdata[1,]

# store up/down regulated DE genes in CSV
subset(resdata, (padj<0.05) & (log2FoldChange>1)) %>% write.csv('clust_1v4_up.csv')
subset(resdata, (padj<0.05) & (log2FoldChange<(-1))) %>% write.csv('clust_1v4_down.csv')


###

# now compare all pairs
fold_change <- function() {
	for (i in 1:length(clusters)) {
		for (j in seq(1, length(clusters))[-i]) {
			clust_i <- paste0('clust', i)
			clust_j <- paste0('clust', j)
			# create data
			data <- cbind(clusters[[i]], clusters[[j]])
			colnames(data) <- c(paste0(clust_i, '_', colnames(data)[1:6]), 
								paste0(clust_j, '_', colnames(data)[7:12]))
			
			# create meta
			meta <- data.frame(cluster=factor(rep(c(clust_i, clust_j), each=6), 
										levels=c(clust_i, clust_j)),
								sample=factor(paste0('fx',1:6),
										levels=paste0('fx',1:6)))
			
			# create DESeq object
			clust_pair <- DESeqDataSetFromMatrix(countData=data, 
												colData=meta, 
												design= ~ cluster+sample)
			
			# normalize DESeq object
			clust_pair <- DESeq(clust_pair, parallel=TRUE)
			
			# obtain results
			res <- results(clust_pair, contrast=c('cluster',clust_i, clust_j))
			
			# add normalized counts to results
			resdata <- merge(as.data.frame(res), 
							as.data.frame(counts(clust_pair,normalize=TRUE)), 
							by='row.names', 
							sort=FALSE)
							
			# store up/down regulated DE genes in CSV
			subset(resdata, (padj<0.05) & (log2FoldChange>1)) %>% write.csv(paste0('DESeq_output/', clust_i, '_vs_', clust_j, '_up.csv'))
			subset(resdata, (padj<0.05) & (log2FoldChange<(-1))) %>% write.csv(paste0('DESeq_output/', clust_i, '_vs_', clust_j, '_down.csv'))
		}
	}
}


fold_change()

# store up/down regulated DE genes in CSV






###

# now try compare clust1 vs clust2&3
i <- 1
j <- (1:3)[-i]

clust_i <- paste0('clust', i)
clust_j <- paste0('clust', j)

data <- cbind(clusters[[i]], clusters[[j[1]]], clusters[[j[2]]])

colnames(data) <- c(paste0(clust_i, '_', colnames(data)[1:6]), 
					paste0(clust_j[1], '_', colnames(data)[7:12]), 
					paste0(clust_j[2], '_', colnames(data)[13:18]))

meta <- data.frame(cluster=factor(c(rep(clust_i, each=6), rep(paste0(clust_j[1],clust_j[2]), each=12)), 
							levels=c(clust_i, paste0(clust_j[1],clust_j[2]))),
					sample=factor(paste0('fx',1:6),
							levels=paste0('fx',1:6)))

clust_pair <- DESeqDataSetFromMatrix(countData=data, 
									colData=meta, 
									design= ~ cluster+sample)

clust_pair <- DESeq(clust_pair, parallel=TRUE)

res <- results(clust_pair, contrast=c('cluster',clust_i, paste0(clust_j[1],clust_j[2])))
res[1,]

resdata <- merge(as.data.frame(res), 
				as.data.frame(counts(clust_pair,normalize=TRUE)), 
				by='row.names', 
				sort=FALSE)
							
subset(resdata, (padj<0.05) & (log2FoldChange>1)) %>% write.csv(paste0('DESeq_output2/', clust_i, '_vs_', clust_j[1], clust_j[2], '_up.csv'))

subset(resdata, (padj<0.05) & (log2FoldChange<(-1))) %>% write.csv(paste0('DESeq_output2/', clust_i, '_vs_', clust_j[1], clust_j[2], '_down.csv'))


# write a loop to compare grouped clusters

fold_change_grouped_clusters <- function() {
	for (i in 1:3) {
		j <- (1:3)[-i]
		clust_i <- paste0('clust', i)
		clust_j <- paste0('clust', j)
		
		# create data
		data <- cbind(clusters[[i]], clusters[[j[1]]], clusters[[j[2]]])
		colnames(data) <- c(paste0(clust_i, '_', colnames(data)[1:6]), 
							paste0(clust_j[1], '_', colnames(data)[7:12]), 
							paste0(clust_j[2], '_', colnames(data)[13:18]))

			
		# create meta
		meta <- data.frame(cluster=factor(c(rep(clust_i, each=6), rep(paste0(clust_j[1],clust_j[2]), each=12)), 
								levels=c(clust_i, paste0(clust_j[1],clust_j[2]))),
							sample=factor(paste0('fx',1:6),
								levels=paste0('fx',1:6)))
			
		# create DESeq object
		clust_pair <- DESeqDataSetFromMatrix(countData=data, 
											colData=meta, 
											design= ~ cluster+sample)
			
		# normalize DESeq object
		clust_pair <- DESeq(clust_pair, parallel=TRUE)
			
		# obtain results
		res <- results(clust_pair, contrast=c('cluster',clust_i, paste0(clust_j[1],clust_j[2])))
			
		# add normalized counts to results
		resdata <- merge(as.data.frame(res), 
						as.data.frame(counts(clust_pair,normalize=TRUE)), 
						by='row.names', 
						sort=FALSE)
							
		# store up/down regulated DE genes in CSV
		subset(resdata, (padj<0.05) & (log2FoldChange>1)) %>% write.csv(paste0('DESeq_output2/', clust_i, '_vs_', clust_j[1], clust_j[2], '_up.csv'))
		subset(resdata, (padj<0.05) & (log2FoldChange<(-1))) %>% write.csv(paste0('DESeq_output2/', clust_i, '_vs_', clust_j[1], clust_j[2], '_down.csv'))

	}
}

fold_change_grouped_clusters()




# rank fold_change by adjusted p-value
# for different cutoffs, just change the subset() criteria

fold_change <- function() {
	for (i in 1:length(clusters)) {
		for (j in seq(1, length(clusters))[-i]) {
			clust_i <- paste0('clust', i)
			clust_j <- paste0('clust', j)
			# create data
			data <- cbind(clusters[[i]], clusters[[j]])
			colnames(data) <- c(paste0(clust_i, '_', colnames(data)[1:6]), 
								paste0(clust_j, '_', colnames(data)[7:12]))
			# create meta
			meta <- data.frame(cluster=factor(rep(c(clust_i, clust_j), each=6), 
										levels=c(clust_i, clust_j)),
								sample=factor(paste0('fx',1:6),
										levels=paste0('fx',1:6)))
			# create DESeq object
			clust_pair <- DESeqDataSetFromMatrix(countData=data, 
												colData=meta, 
												design= ~ cluster+sample)
			# normalize DESeq object
			clust_pair <- DESeq(clust_pair, parallel=TRUE)
			# obtain results
			res <- results(clust_pair, contrast=c('cluster',clust_i, clust_j))
			# order results rows by padj i.e. adjusted p-value
			res <- res[order(res$padj),]	
			# add normalized counts to results
			resdata <- merge(as.data.frame(res), 
							as.data.frame(counts(clust_pair,normalize=TRUE)), 
							by='row.names', 
							sort=FALSE)				
			# store up/down regulated DE genes in CSV
			subset(resdata, (padj<0.05) & (log2FoldChange>0.5849625)) %>% write.csv(paste0('DESeq_output_fold_1.5/', clust_i, '_vs_', clust_j, '_up.csv'))
			#subset(resdata, log2FoldChange<(-1)) %>% write.csv(paste0('DESeq_output_fold_1.5/', clust_i, '_vs_', clust_j, '_down.csv'))
		}
	}
}

fold_change()

fold_change_grouped_clusters <- function() {
	for (i in 1:3) {
		j <- (1:3)[-i]
		clust_i <- paste0('clust', i)
		clust_j <- paste0('clust', j)
		# create data
		data <- cbind(clusters[[i]], clusters[[j[1]]], clusters[[j[2]]])
		colnames(data) <- c(paste0(clust_i, '_', colnames(data)[1:6]), 
							paste0(clust_j[1], '_', colnames(data)[7:12]), 
							paste0(clust_j[2], '_', colnames(data)[13:18]))
		# create meta
		meta <- data.frame(cluster=factor(c(rep(clust_i, each=6), rep(paste0(clust_j[1],clust_j[2]), each=12)), 
								levels=c(clust_i, paste0(clust_j[1],clust_j[2]))),
							sample=factor(paste0('fx',1:6),
								levels=paste0('fx',1:6)))
		# create DESeq object
		clust_pair <- DESeqDataSetFromMatrix(countData=data, 
											colData=meta, 
											design= ~ cluster+sample)
		# normalize DESeq object
		clust_pair <- DESeq(clust_pair, parallel=TRUE)
		# obtain results
		res <- results(clust_pair, contrast=c('cluster',clust_i, paste0(clust_j[1],clust_j[2])))
		# add normalized counts to results
		resdata <- merge(as.data.frame(res), 
						as.data.frame(counts(clust_pair,normalize=TRUE)), 
						by='row.names', 
						sort=FALSE)
		# store up/down regulated DE genes in CSV
		subset(resdata, (padj<0.05) & (log2FoldChange>0.5849625)) %>% write.csv(paste0('DESeq_output_fold_1.5/', clust_i, '_vs_', clust_j[1], clust_j[2], '_up.csv'))
		subset(resdata, (padj<0.05) & (log2FoldChange<(-1))) %>% write.csv(paste0('DESeq_output2_fold_1.5/', clust_i, '_vs_', clust_j[1], clust_j[2], '_down.csv'))
	}
}

fold_change_grouped_clusters()

# gotta find a way to enter eval into a func!



