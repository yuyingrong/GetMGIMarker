
# This is a cleaned version of the split-to-clusters-and-merge script
# Splitting merged dataset based on presence/absence of markers and re-merge

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

setwd('~/Labs/xinhua/scSeq/')

for (n in 1:6) {
	obj <- paste0('fx',n)
	assign(obj, Read10X(data.dir = sprintf('FX_%d/filtered_feature_bc_matrix/',n))) 
	assign(obj, CreateSeuratObject(counts=get(obj), project=obj))
	# should I normalize?
}

obj <- NULL

merged <- merge(
	x=fx1, y=c(fx2,fx3,fx4,fx5,fx6), 
	add.cell.ids=c('fx1','fx2','fx3','fx4','fx5','fx6'), 
	project='fx'
)

fx1 <- NULL
fx2 <- NULL
fx3 <- NULL
fx4 <- NULL
fx5 <- NULL
fx6 <- NULL

# rename orig.ident as samples
merged@meta.data <- merged@meta.data %>% dplyr::rename(sample=orig.ident)

# add a col of cellnames; unnecessary
# merged@meta.data$cells <- rownames(merged@meta.data)

# add a col of clusters
merged@meta.data$cluster <- NA

merged@meta.data[1:5,]
length(rownames(merged@meta.data))
# 38164
merged@meta.data[38160:38164,]

# add cluster info

# for cells in $sample fx1:fx6

merged@meta.data$cluster[
	which(GetAssayData(object = merged, slot = "data")['PAX7',]>0 &
		GetAssayData(object = merged, slot = "data")['MYF5',]==0 &
		GetAssayData(object = merged, slot = "data")['MYOD1',]==0 &
		GetAssayData(object = merged, slot = "data")['MYOG',]==0)
] <- 'clust1'

length(which(merged@meta.data$cluster=='clust1'))
# 813, good

# add the rest of clusters

merged@meta.data$cluster[
	which(GetAssayData(object = merged, slot = "data")['PAX7',]==0 &
		GetAssayData(object = merged, slot = "data")['MYF5',]>0 &
		GetAssayData(object = merged, slot = "data")['MYOD1',]==0 &
		GetAssayData(object = merged, slot = "data")['MYOG',]==0)
] <- 'clust2'

merged@meta.data$cluster[
	which(GetAssayData(object = merged, slot = "data")['PAX7',]==0 &
		GetAssayData(object = merged, slot = "data")['MYF5',]==0 &
		GetAssayData(object = merged, slot = "data")['MYOD1',]>0 &
		GetAssayData(object = merged, slot = "data")['MYOG',]==0)
] <- 'clust7'

merged@meta.data$cluster[
	which(GetAssayData(object = merged, slot = "data")['PAX7',]>0 &
		GetAssayData(object = merged, slot = "data")['MYF5',]>0 &
		GetAssayData(object = merged, slot = "data")['MYOD1',]==0 &
		GetAssayData(object = merged, slot = "data")['MYOG',]==0)
] <- 'clust3'

merged@meta.data$cluster[
	which(GetAssayData(object = merged, slot = "data")['PAX7',]>0 &
		GetAssayData(object = merged, slot = "data")['MYF5',]==0 &
		GetAssayData(object = merged, slot = "data")['MYOD1',]>0 &
		GetAssayData(object = merged, slot = "data")['MYOG',]==0)
] <- 'clust4'

merged@meta.data$cluster[
	which(GetAssayData(object = merged, slot = "data")['PAX7',]==0 &
		GetAssayData(object = merged, slot = "data")['MYF5',]>0 &
		GetAssayData(object = merged, slot = "data")['MYOD1',]>0 &
		GetAssayData(object = merged, slot = "data")['MYOG',]==0)
] <- 'clust5'

merged@meta.data$cluster[
	which(GetAssayData(object = merged, slot = "data")['PAX7',]>0 &
		GetAssayData(object = merged, slot = "data")['MYF5',]>0 &
		GetAssayData(object = merged, slot = "data")['MYOD1',]>0 &
		GetAssayData(object = merged, slot = "data")['MYOG',]==0)
] <- 'clust6'

merged@meta.data$cluster[
	which(is.na(merged@meta.data$cluster))
] <- 'unassigned'

merged@meta.data[1:5,]
length(which(merged@meta.data$cluster=='clust2'))
# 3584
length(which(merged@meta.data$cluster=='clust3'))
# 481
length(which(merged@meta.data$cluster=='clust4'))
# 3983
length(which(merged@meta.data$cluster=='clust5'))
# 4917
length(which(merged@meta.data$cluster=='clust6'))
# 6581
length(which(merged@meta.data$cluster=='clust7'))
# 4961

ls()

DotPlot(merged, features=c('PAX7','MYF5','MYOD1','MYOG'), group.by='cluster')
# returned error! why?

# randomly looked into data
merged@meta.data$cluster[37900:37920]
# there are NA's! perhaps they are the cause

merged@meta.data[37983,]
# yes, the cluster slot is NA

# use is.na() to get the indices of NA's; ...==NA won't work
merged@meta.data$cluster[
	which(is.na(merged@meta.data$cluster))
] <- 'unassigned'

merged@meta.data[38159,]
# okay, now there it has 'unassigned' for $cluster, intead of NA

DotPlot(merged, features=c('PAX7','MYF5','MYOD1','MYOG'), group.by='cluster')
# worked


### Quality Control ###

merged$percent_mt <- PercentageFeatureSet(merged, pattern = '^MT-')
merged$log10GenesPerUMI <- log10(merged@meta.data$nFeature_RNA)/log10(merged@meta.data$nCount_RNA)

# looking at dataset by cluster

VlnPlot(merged, features = c('nFeature_RNA','nCount_RNA'), group.by='cluster', ncol=2)
VlnPlot(merged, features = c('percent_mt','log10GenesPerUMI'), group.by='cluster', ncol=2)

merged@meta.data %>% ggplot(aes(x=nCount_RNA,y=nFeature_RNA,color=percent_mt))+
	geom_point(alpha=0.5)+
	scale_color_gradient(low='gray',high='black')+
	stat_smooth(method=lm)+
	scale_x_log10()+
	scale_y_log10()+
	theme_classic()+
	geom_vline(xintercept=500)+
	geom_hline(yintercept=250)+
	facet_wrap(~cluster)

# trim cells	
# merged <- subset(merged,subset=percent_mt<50)
# does not seem necessary, since most high mt% cells are in 'unassigned' cluster

# looking at dataset by sample
	
DotPlot(merged, features=c('PAX7','MYF5','MYOD1','MYOG'), group.by='sample')
VlnPlot(merged, features = c('nFeature_RNA','nCount_RNA'), group.by='sample', ncol=2)
VlnPlot(merged, features = c('percent_mt','log10GenesPerUMI'), group.by='sample', ncol=2)

merged@meta.data %>% ggplot(aes(x=nCount_RNA,y=nFeature_RNA,color=percent_mt))+
	geom_point(alpha=0.5)+
	scale_color_gradient(low='gray',high='black')+
	stat_smooth(method=lm)+
	scale_x_log10()+
	scale_y_log10()+
	theme_classic()+
	geom_vline(xintercept=500)+
	geom_hline(yintercept=250)+
	facet_wrap(~sample)

# save obj
saveRDS(object=merged, file='rds/merged.rds')


# I feel this could be done an easier way
subset(merged,cells=rownames(merged@meta.data)[which(merged@meta.data$sample!=c('fx2','fx6'))]) %>% DotPlot(features=c('CDKN1A','CDKN1B','CDKN1C','TP53','RBL2','CCND1','CCND3','MKI67','ANLN','BIRC5','CCNA2','CCNB1','CCNE2','CGOL1','CENPE','TOP2A','UBE2C','CDK1','CDK2','CDC6','CDC20'), group.by='cluster') + coord_flip() + ggtitle(label='all but fx2,6')

