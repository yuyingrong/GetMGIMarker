
# Yuying Rong
# 2022

# this script visualizes scRNA-seq output for quality control 
# before and after filtering single cells
# it also visualizes marker gene expression levels
# inspired by: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

setwd('~/Labs/xinhua/scSeq/')
data <- Read10X(data.dir = 'FX_6/filtered_feature_bc_matrix/')
data <- CreateSeuratObject(counts = data, project = 'fx6')


#### Quality Control #####

data$percent_mt <- PercentageFeatureSet(data, pattern = '^MT-')
data$log10GenesPerUMI <- log10(data@meta.data$nFeature_RNA)/log10(data@meta.data$nCount_RNA)
data@meta.data <- data@meta.data %>% dplyr::rename(sample=orig.ident)


## before trimming

# open PDF for QC plots
pdf(file='FX6_original.pdf',width=10,height=10)

# check if nCells matches expectation based on nCells loaded and lib prep method used
data@meta.data %>% ggplot(aes(x=sample, fill=sample))+
	geom_bar()+
	theme_classic()+
	theme(axis.text.x=element_text(angle=45,hjust=1))+
	theme(plot.title=element_text(hjust=0.5,face='bold'))+
	ggtitle('nCells')

# nFeature_RNA stands for the number of genes detected in each cell
# nCount_RNA stands for the number of molecules (mRNA) detected in a cell
# low nFeature_RNA: dying cell, empty droplet
# very high nFeature_RNA or nCount_RNA: multiplet?
# percent_mt: mitochondria_gene%: high: necrotic cell or higher metabolic activity cell
VlnPlot(data, features = c('nFeature_RNA','nCount_RNA'),ncol=2)
VlnPlot(data, features = c('percent_mt','log10GenesPerUMI'),ncol=2)

# another view of above matrix
data@meta.data %>% ggplot(aes(x=nCount_RNA,y=nFeature_RNA,color=percent_mt))+
	geom_point(alpha=0.5)+
	scale_color_gradient(low='gray',high='black')+
	stat_smooth(method=lm)+
	scale_x_log10()+
	scale_y_log10()+
	theme_classic()+
	geom_vline(xintercept=500)+
	geom_hline(yintercept=250)+
	facet_wrap(~sample)

# finding the best cutoff thresholds
ggplot(data@meta.data, aes(x=percent_mt))+
	geom_histogram(binwidth=0.5)
	
ggplot(data@meta.data, aes(x=nCount_RNA))+
	geom_histogram(binwidth=500)

ggplot(data@meta.data, aes(x=nFeature_RNA))+
	geom_histogram(binwidth=50)

dev.off()


## begin trimming

pdf(file='FX6_trimmed.pdf',width=10,height=10)

summary(data@meta.data$percent_mt)
summary(data@meta.data$nCount_RNA)
summary(data@meta.data$nFeature_RNA)

# trim the visible high percent_mt cluster
data <- subset(data,subset=percent_mt<30)
### I feel more could be trimmed!! view the summary()

summary(data@meta.data$percent_mt)

# plots
# check if nCells matches expectation based on nCells loaded and lib prep method used
data@meta.data %>% ggplot(aes(x=sample, fill=sample))+
	geom_bar()+
	theme_classic()+
	theme(axis.text.x=element_text(angle=45,hjust=1))+
	theme(plot.title=element_text(hjust=0.5,face='bold'))+
	ggtitle('nCells')

VlnPlot(data, features = c('nFeature_RNA','nCount_RNA'),ncol=2)
VlnPlot(data, features = c('percent_mt','log10GenesPerUMI'),ncol=2)

data@meta.data %>% ggplot(aes(x=nCount_RNA,y=nFeature_RNA,color=percent_mt))+
	geom_point(alpha=0.5)+
	scale_color_gradient(low='gray',high='black')+
	stat_smooth(method=lm)+
	scale_x_log10()+
	scale_y_log10()+
	theme_classic()+
	geom_vline(xintercept=500)+
	geom_hline(yintercept=250)+
	facet_wrap(~sample)

ggplot(data@meta.data, aes(x=percent_mt))+
	geom_histogram(binwidth=0.5)
	
ggplot(data@meta.data, aes(x=nCount_RNA))+
	geom_histogram(binwidth=500)

ggplot(data@meta.data, aes(x=nFeature_RNA))+
	geom_histogram(binwidth=50)


### custom cell sorting ###

sum(GetAssayData (object = data, slot = "data") ['PAX7',]>0)
# 1321
sum(GetAssayData (object = data, slot = "data") ['MYF5',]>0)
# 2054
sum(GetAssayData (object = data, slot = "data") ['MYOG',]>0)
# 5
sum(GetAssayData (object = data, slot = "data") ['MYOD1',]>0)
# 2349

# scatterplots
plot(sort(GetAssayData (object = data, slot = "data") ['PAX7',]), xlab='cells', ylab='Pax7 Levels')
plot(sort(GetAssayData (object = data, slot = "data") ['MYF5',]), xlab='cells', ylab='Myf5 Levels')
plot(sort(GetAssayData (object = data, slot = "data") ['MYOG',]), xlab='cells', ylab='MyoG Levels')
plot(sort(GetAssayData (object = data, slot = "data") ['MYOD1',]), xlab='cells', ylab='MyoD Levels')


# histograms
hist(sort(GetAssayData (object = data, slot = "data") ['PAX7',]), breaks=100, xlab='Levels', ylab='Frequency', main='Histogram of Pax7 Levels (FX1)')

hist(sort(GetAssayData (object = data, slot = "data") ['PAX7',]), breaks=100, xlab='Levels', ylab='Frequency', main='Histogram of Pax7 Levels (FX1)', ylim=c(0,500))

hist(sort(GetAssayData (object = data, slot = "data") ['MYF5',]), breaks=100, xlab='Levels', ylab='Frequency', main='Histogram of Myf5 Levels (FX1)')

hist(sort(GetAssayData (object = data, slot = "data") ['MYF5',]), breaks=100, xlab='Levels', ylab='Frequency', main='Histogram of Myf5 Levels (FX1)', ylim=c(0,1600))

hist(sort(GetAssayData (object = data, slot = "data") ['MYOG',]), breaks=100, xlab='Levels', ylab='Frequency', main='Histogram of MyoG Levels (FX1)')

hist(sort(GetAssayData (object = data, slot = "data") ['MYOG',]), breaks=100, xlab='Levels', ylab='Frequency', main='Histogram of MyoG Levels (FX1)', ylim=c(0,20))

hist(sort(GetAssayData (object = data, slot = "data") ['MYOD1',]), breaks=100, xlab='Levels', ylab='Frequency',main='Histogram of MyoD Levels (FX1)')

hist(sort(GetAssayData (object = data, slot = "data") ['MYOD1',]), breaks=100, xlab='Levels', ylab='Frequency',main='Histogram of MyoD Levels (FX1)', ylim=c(0,500))


# dotplot
DotPlot(data, features=c('PAX7','MYF5','MYOD1','MYOG'))


# find where interested markers are stored
match('PAX7', rownames(GetAssayData (object = data, slot = "data")))
# 427
match('MYF5', rownames(GetAssayData (object = data, slot = "data")))
# 22200
match('MYOD1', rownames(GetAssayData (object = data, slot = "data")))
# 19328
match('MYOG', rownames(GetAssayData (object = data, slot = "data")))
# 2774

# create a smaller expression matrix with only the concerned markers
exp <- rbind(GetAssayData (object = data, slot = "data")[427,], 
GetAssayData (object = data, slot = "data")[22200,], 
GetAssayData (object = data, slot = "data")[19328,], 
GetAssayData (object = data, slot = "data")[2774,])
rownames(exp) <- c('PAX7','MYF5','MYOD1','MYOG')

# sort matrix by a column in ascending order
# M1[order(M1[,1],decreasing=FALSE),]
# this can be applied to sort all rows by the order of one particular row
exp_ordered <- exp[,order(exp[1,],decreasing=TRUE)]

# plot expression levels
plot(exp[1,], xlab='Cells', ylab='Levels', main='Pax7')
plot(exp[2,], xlab='Cells', ylab='Levels', main='Myf5')
plot(exp[3,], xlab='Cells', ylab='Levels', main='MyoD')
plot(exp[4,], xlab='Cells', ylab='Levels', main='MyoG')

dev.off()


### Summarizing data points in a table ###
length(exp[1,])
length(exp['PAX7',])
# number of trimmed cells: 11964 (fx1: 6262)
length(which(exp[1,]==0))
# number of Pax7 expression == 0: 11581
length(which(exp[1,]>0))
# number of Pax7 expression > 0: 383
# 11581+383 == 11964
# good

# view how many cells return 0 and how many >0
for ( m in c('PAX7','MYF5','MYOD1','MYOG') ) { 
	print(paste0(m,'==0: ',length(which(exp[m,]==0))))
	print(paste0(m,'>0: ',length(which(exp[m,]>0))))
	}

# create a logical expression matrix (TRUE==1 and FALSE==0)
# create a matrix containing all 16 situations
mtx <- cbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1),c(1,1,0,0),c(1,0,1,0),c(1,0,0,1),c(0,1,1,0),c(0,1,0,1),c(0,0,1,1),c(1,1,1,0),c(0,1,1,1),c(1,0,1,1),c(1,1,0,1),c(1,1,1,1),c(0,0,0,0))

# create an count vector storing the number of cells counted for each situation
output <- rep(0L, 16)

# count cells for all 16 scenarios (O(n)=16n)
for ( cell in 1:ncol(exp) ) {
	for ( i in 1:ncol(mtx) ) {
		if (identical((exp[,cell]>rbind(0,0,0,0))==mtx[,i],rbind(TRUE, TRUE, TRUE, TRUE))) {
			output[i] <- output[i] + 1
		}
	}
}

# print out the results
for (c in output) {
	print(c)
	}






