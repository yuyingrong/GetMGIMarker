
# 20231008
# Yuying Rong
# clean code for all tasks


# import input file

file_path <- '/Users/shell/MSc/LMU/Schneeberger'

import_cellsnp_output <- function (x) {
  # import cellsnp output mmmatrix
  mtx <- Matrix::readMM(paste0(file_path, '/data/', x, '/cellSNP.tag.DP.mtx'))
  m <- cbind.data.frame(row=mtx@i+1, col=mtx@j+1, x=mtx@x)
  m <- m[,-3]
  colnames(m) <- c('position_index','barcode_index')
  rm(mtx)
  
  # write m to file
  #write.table(m, 
  #            file=paste0(file_path, '/data/', x, '/m.csv'), 
  #            quote=F, row.names=F, col.names=F)
  return(m)
}



# process input file

make_count_table <- function (x,m) {
  # make count table t for number of markers covered (position_indices) for each barcode_index
  t <- table(m$barcode_index)
  # rownames of t are barcode_indices
  t <- data.frame(t)
  colnames(t) <- c('barcode_index','occurrence')
  
  # added: additional tables to write
  write.table(t[t$occurrence >= 1 & t$occurrence < 10,1], 
              file=paste0(file_path, '/data/', x, '/t1to10_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  write.table(t[t$occurrence >= 10 & t$occurrence < 100,1], 
              file=paste0(file_path, '/data/', x, '/t10to100_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  write.table(t[t$occurrence >= 300 & t$occurrence < 700,1], 
              file=paste0(file_path, '/data/', x, '/t300to700_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  write.table(t[t$occurrence >= 700 & t$occurrence < 1000,1], 
              file=paste0(file_path, '/data/', x, '/t700to1000_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  write.table(t[t$occurrence >= 10000,1], 
              file=paste0(file_path, '/data/', x, '/t10000_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  write.table(t[t$occurrence >= 1000 & t$occurrence < 5000,1], 
              file=paste0(file_path, '/data/', x, '/t1000to5000_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  write.table(t[t$occurrence >= 5000 & t$occurrence < 10000,1], 
              file=paste0(file_path, '/data/', x, '/t5000to10000_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  #return(t)
  
  # write groups of barcode_indices to files
  write.table(t[t$occurrence >= 100,1], 
              file=paste0(file_path, '/data/', x, '/t100_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  write.table(t[t$occurrence >= 100 & t$occurrence < 300,1], 
              file=paste0(file_path, '/data/', x, '/t100to300_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  write.table(t[t$occurrence >= 300 & t$occurrence < 1000,1], 
              file=paste0(file_path, '/data/', x, '/t300to1000_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  write.table(t[t$occurrence >= 1000,1], 
              file=paste0(file_path, '/data/', x, '/t1000_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  return(t)
}



# Task 1

plot_num_markers_per_barcode <- function (x,t) {
  # change barcode_index into factors
  t$barcode_index <- as.factor(t$barcode_index)
  # plot histogram of number of markers mapped by (occurrence per) each barcode_index
  pdf(paste0(file_path, '/out/', x, '_hist_marker_mappings_per_barcode.pdf'))
  
  hist(log10(t[t$occurrence>=1,2]),
       xlab='log10(markers per barcode)',
       main='Distribution of marker counts')
  hist(log10(t[t$occurrence>1,2]),
       xlab='log10(markers per barcode)',
       main='Distribution of marker counts\n(for barcodes with >1 markers)')
  
  dev.off()
}



# Task 6

# plot marker density: number of markers per region

# input: mp, t_group barcode indices (txt generated from function make_count_table)
# inherit: m
plot_marker_position <- function (x,m,t) {
  pdf(paste0(file_path, '/out/', x, '_marker_density0.pdf'))
  mp <- read.table(paste0(file_path, '/data/', x, '/cellSNP.markers.tsv'), header=F, sep='\n')
  mp <- mp$V1
  
  # max(mp)=85139182 and 85139182/100000=851.4
  # each bin represents 100000 nucleotides
  hist(mp, breaks=852, labels=F,
       xlab='Markers on 100,000-bp bins',
       main='Marker density for all barcodes')
  
  # trim t and m (if m is not by default 100+)
  #t <- t[t$occurrence>=100,]
  #m <- subset(m, barcode_index %in% t[t$occurrence>=100,1])
  
  # for 1~10 group
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t1to10_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(mp[subset(m, barcode_index %in% barcode_indices)$position_index], breaks=852, labels=F,
       xlab='Markers on 100,000-bp bins',
       main='Marker density\n(for barcodes with 1~9 markers)')
  
  # for 10~100 group
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t10to100_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(mp[subset(m, barcode_index %in% barcode_indices)$position_index], breaks=852, labels=F,
       xlab='Markers on 100,000-bp bins',
       main='Marker density\n(for barcodes with 10~99 markers)')
  
  # plot for 100+ group
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t100_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(mp[subset(m, barcode_index %in% barcode_indices)$position_index], breaks=852, labels=F,
       xlab='Markers on 100,000-bp bins',
       main='Marker density\n(for barcodes with >=100 markers)')
  
  # plot for 300+ group
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t300_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(mp[subset(m, barcode_index %in% barcode_indices)$position_index], breaks=852, labels=F,
       xlab='Markers on 100,000-bp bins',
       main='Marker density\n(for barcodes with >=300 markers)')
  
  # plot for 100~300 group
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t100to300_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(mp[subset(m, barcode_index %in% barcode_indices)$position_index], breaks=852, labels=F,
       xlab='Markers on 100,000-bp bins',
       main='Marker density\n(for barcodes with 100~299 markers)')
  
  # 300~1000
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t300to1000_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(mp[subset(m, barcode_index %in% barcode_indices)$position_index], breaks=852, labels=F,
       xlab='Markers on 100,000-bp bins',
       main='Marker density\n(for barcodes with 300~999 markers)')
  
  # 1000+
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t1000_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(mp[subset(m, barcode_index %in% barcode_indices)$position_index], breaks=852, labels=F,
       xlab='Markers on 100,000-bp bins',
       main='Marker density\n(for barcodes with >=1000 markers)')
  
  # 10000+
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t10000_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(mp[subset(m, barcode_index %in% barcode_indices)$position_index], breaks=852, labels=F,
       xlab='Markers on 100,000-bp bins',
       main='Marker density\n(for barcodes with >=10000 markers)')

  # 1000-5000
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t1000to5000_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with 1000~4999 markers)')
  
  # 5000-10000
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t5000to10000_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with 5000~9999 markers)')
  
  dev.off()
}



# Task 2 part 1

# crossover resolution, or distance between adjacent markers

prep_resolution <- function (x,m) {
  # import marker_positions array mp from bash formatted table
  # row_number of mp correspond to position_indices of m
  mp <- read.table(paste0(file_path, '/data/', x, '/cellSNP.markers.tsv'), header=F, sep='\n')
  mp <- mp$V1
  
  # swap position_indices in m with marker_positions from mp
  # duplicate m (m is not sorted yet)
  m_sorted <- m
  m_sorted$position <- mp[as.numeric(m_sorted$position_index)]
  m_sorted$position <- as.integer(m_sorted$position)
  m_sorted <- m_sorted[,-1]
  
  # sort position_index and barcode_index table m by barcode_index then by position_index
  m_sorted <- m_sorted[order(m_sorted$barcode_index, m_sorted$position),]
  
  # write m_sorted to file for import to python
  write.table(m_sorted, file=paste0(file_path, '/data/', x, '/m_sorted.csv'), sep=',', row.names=F)
  rm(m_sorted)
}



# Task 2 part 2

# after running python: "20230731_calc_resolution.py"
# plot from output

plot_resolution <- function (x) {
  # import res from python output
  res <- read.table(paste0(file_path, '/out/', x, '_res.csv'), header=T, sep=',')
  
  pdf(paste0(file_path, '/out/', x, '_marker_distance.pdf'))

  # plot for all
  hist(log10(res$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of marker distances for all barcodes')
  
  #hist(res$marker_distance, breaks=50,
  #     xlab='Nucleotide distances between adjacent markers',
  #     main='Distribution of marker distances for all barcodes')
  
  # plot for groups
  
  # for multiplets
  multiplet_barcode_indices <- read.table(paste0(file_path, '/data/', x, '/multiplet_barcode_indices.txt'), header=F, sep='\n')
  multiplet_barcode_indices <- multiplet_barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% multiplet_barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for multiplets by AMULET)')
  # for all but multiplets
  #`%notin%` <- function(x,y) !(x %in% y)
  #hist(log10(subset(res, barcode_index %notin% multiplet_barcode_indices)$marker_distance), breaks=50,
  #     xlab='log10(nucleotide distances between adjacent markers)',
  #     main='Distribution of CO marker resolution\n(for all but multiplets labeled by AMULET)')
  
  # for >=2/>2 mean coverage depth
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/avg_mkr_cvg_dep_over_equal_2_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with >=2 mean marker coverage depth)')
  
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/avg_mkr_cvg_dep_over_2_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with >2 mean marker coverage depth)')
  
  # for 1~10
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t1to10_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with 1~9 markers)')
  
  # for 10~100
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t10to100_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with 10~99 markers)')
  
  # for 100+ group
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t100_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  # making res smaller(!)
  res <- subset(res, barcode_index %in% barcode_indices)
  hist(log10(res$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with >=100 markers)')

  # for 100~300 group
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t100to300_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with 100~299 markers)')
  
  # for 300+ group
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t300_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  # making res smaller(!)
  res <- subset(res, barcode_index %in% barcode_indices)
  hist(log10(res$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with >=300 markers)')
  
  # for 300+ group without multiplets
  `%notin%` <- function(x,y) !(x %in% y)
  hist(log10(subset(res, barcode_index %notin% multiplet_barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with >=300 markers without AMULET labeled multiplets)')
  
  # for 300~1000
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t300to1000_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with 300~999 markers)')
  
  # for 300~700
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t300to700_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with 300~699 markers)')
  
  # for 700~1000
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t700to1000_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with 700~999 markers)')
  
  # for 1000+
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t1000_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with >=1000 markers)')
  
  # for 10000+
  barcode_indices <- read.table(paste0(file_path, '/data/', x, '/t10000_barcode_indices.txt'), header=F, sep='\n')
  barcode_indices <- barcode_indices$V1
  hist(log10(subset(res, barcode_index %in% barcode_indices)$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of adjacent marker distances\n(for barcodes with >=10000 markers)')
  
  # save all plots
  dev.off()
  rm(res)
}



# Task 5
# barcode multiplet detection (python)
# cellular multiplet detection
# check marker distributions of barcodes of cellular multiplets

# inherit var m (position_id,barcode_id,count; without subsetting to >=100 mappings per barcode)
# inherit var mp (array of marker positions, index corresponds to position_index in var m)
# additional input: avgs
plot_marker_coverage_depth <- function (x, m) {
  
  # import marker positions, avgs
  mp <- read.table(paste0(file_path, '/data/', x, '/cellSNP.markers.tsv'), header=F, sep='\n')
  mp <- mp$V1
  m$position <- mp[m$position_index]
  print(x)
  print(identical(length(unique(m$position_index)),length(mp)))# should print TRUE

  avgs <- read.table(paste0(file_path, '/out/', x, '_avgs.csv'), header=T, sep=',')
  colnames(avgs) <- c('barcode_index','avg_depth')# it was mislabeled in the shell script
  
  # write mean_marker_cvg_depth barcode_indices to file
  write.table(avgs[avgs$avg_depth>=2,1], 
              file=paste0(file_path, '/data/', x, '/avg_mkr_cvg_dep_over_equal_2_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  
  write.table(avgs[avgs$avg_depth>2,1], 
              file=paste0(file_path, '/data/', x, '/avg_mkr_cvg_dep_over_2_barcode_indices.txt'), 
              quote=F, row.names=F, col.names=F)
  
  # plot mean marker coverage depth per barcode and marker position along the genome for groups
  pdf(paste0(file_path, '/out/', x, '_mean_marker_coverage_depth.pdf'))
  
  # for all
  hist(avgs$avg_depth, 
       labels = F,
       xlab = 'Mean marker coverage depth per barcode', 
       main = 'Histogram of mean marker coverage depth per barcode')
  
  # for cells with >2 coverage depth
  hist(avgs[avgs$avg_depth>2,2], 
       labels = T,
       xlab = 'Mean marker coverage depth per barcode', 
       main = 'Histogram of mean (>2) marker coverage depth per barcode')
  
  bi <- subset(avgs, avg_depth>2)[,1]
  hist(subset(m, barcode_index %in% bi)[,3],
       breaks=852,
       xlab = 'Marker positions along the genome',
       main = 'Marker positions for barcodes with >2 mean marker coverage depth')
  
  # for cells with >=2 coverage depth
  hist(avgs[avgs$avg_depth>=2,2], 
       labels = T,
       xlab = 'Mean marker coverage depth per barcode', 
       main = 'Histogram of mean (>=2) marker coverage depth per barcode')
  
  bi <- subset(avgs, avg_depth>=2)[,1]
  hist(subset(m, barcode_index %in% bi)[,3],
       breaks=852,
       xlab = 'Marker positions along the genome',
       main = 'Marker positions for barcodes with >=2 mean marker coverage depth')
  
  # for AMULET multiplets
  bi <- read.table(paste0(file_path, '/data/', x, '/multiplet_barcode_indices.txt'), header=F, sep='\n')
  bi <- bi$V1
  
  hist(subset(avgs, barcode_index %in% bi)[,2],
       labels = T,
       xlab = 'Mean marker coverage depth per barcode', 
       main = 'Histogram of mean marker coverage depth per barcode\nfor barcodes labeled as multiplets by AMULET')
  
  hist(subset(m, barcode_index %in% bi)[,3],
       breaks=852,
       xlab = 'Marker positions along the genome',
       main = 'Marker positions for barcodes labeled as multiplets by AMULET')

  # for 100+ group
  bi <- read.table(paste0(file_path, '/data/', x, '/t100_barcode_indices.txt'), header=F, sep='\n')
  bi <- bi$V1
  
  hist(subset(avgs, barcode_index %in% bi)[,2],
       labels = T,
       xlab = 'Mean marker coverage depth per barcode', 
       main = 'Histogram of mean marker coverage depth per barcode\nfor barcodes with >=100 markers')
  
  hist(subset(m, barcode_index %in% bi)[,3],
       breaks=852,
       xlab = 'Marker positions along the genome',
       main = 'Marker positions for barcodes with >=100 markers')
  
  # for 100~300 group
  bi <- read.table(paste0(file_path, '/data/', x, '/t100to300_barcode_indices.txt'), header=F, sep='\n')
  bi <- bi$V1
  
  hist(subset(avgs, barcode_index %in% bi)[,2],
       labels = T,
       xlab = 'Mean marker coverage depth per barcode', 
       main = 'Histogram of mean marker coverage depth per barcode\nfor barcodes with 100~300 markers')
  
  hist(subset(m, barcode_index %in% bi)[,3],
       breaks=852,
       xlab = 'Marker positions along the genome',
       main = 'Marker positions for barcodes with 100~300 markers')
  
  # for 300+ group
  bi <- read.table(paste0(file_path, '/data/', x, '/t300_barcode_indices.txt'), header=F, sep='\n')
  bi <- bi$V1
  
  hist(subset(avgs, barcode_index %in% bi)[,2],
       labels = T,
       xlab = 'Mean marker coverage depth per barcode', 
       main = 'Histogram of mean marker coverage depth per barcode\nfor barcodes with >=300 markers')
  
  hist(subset(m, barcode_index %in% bi)[,3],
       breaks=852,
       xlab = 'Marker positions along the genome',
       main = 'Marker positions for barcodes with >=300 markers')
  
  # for 300~1000
  bi <- read.table(paste0(file_path, '/data/', x, '/t300to1000_barcode_indices.txt'), header=F, sep='\n')
  bi <- bi$V1
  
  hist(subset(avgs, barcode_index %in% bi)[,2],
       labels = T,
       xlab = 'Mean marker coverage depth per barcode', 
       main = 'Histogram of mean marker coverage depth per barcode\nfor barcodes with 300~1000 markers')
  
  hist(subset(m, barcode_index %in% bi)[,3],
       breaks=852,
       xlab = 'Marker positions along the genome',
       main = 'Marker positions for barcodes with 300~1000 markers')
  
  # for 1000+
  bi <- read.table(paste0(file_path, '/data/', x, '/t1000_barcode_indices.txt'), header=F, sep='\n')
  bi <- bi$V1
  
  hist(subset(avgs, barcode_index %in% bi)[,2],
       labels = T,
       xlab = 'Mean marker coverage depth per barcode', 
       main = 'Histogram of mean marker coverage depth per barcode\nfor barcodes with >=1000 markers')
  
  hist(subset(m, barcode_index %in% bi)[,3],
       breaks=852,
       xlab = 'Marker positions along the genome',
       main = 'Marker positions for barcodes with >=1000 markers')
  
  dev.off()
}



# Task 3&4
# based on cellranger output,
# calc breadth of coverage and depth of coverage,
# make histograms

# create var tb: barcode,occurrence
tb <- readRDS('/Users/shell/MSc/LMU/Schneeberger/tb.rds')
head(tb)
tb <- tb[,-3]# rm third col breadth, which has not finished calc due to memory
saveRDS(tb, '/Users/shell/MSc/LMU/Schneeberger/tb.rds')
nrow(tb[tb$occurrence >= 1000,])# 2842
nrow(tb[tb$occurrence >= 300,])# 12501
nrow(tb[tb$occurrence >= 100,])# 46090
# cluster >=100
write.table(tb[tb$occurrence >= 100,1], 
            file='/Users/shell/MSc/LMU/Schneeberger/tb100_barcodes.txt', 
            quote=F, row.names=F, col.names=F)
# cluster 100-300
write.table(tb[tb$occurrence >= 100 & tb$occurrence < 300,1], 
            file='/Users/shell/MSc/LMU/Schneeberger/tb100to300_barcodes.txt', 
            quote=F, row.names=F, col.names=F)
# cluster 300-1000
write.table(tb[tb$occurrence >= 300 & tb$occurrence < 1000,1], 
            file='/Users/shell/MSc/LMU/Schneeberger/tb300to1000_barcodes.txt', 
            quote=F, row.names=F, col.names=F)
# cluster >=1000
write.table(tb[tb$occurrence >= 1000,1], 
            file='/Users/shell/MSc/LMU/Schneeberger/tb1000_barcodes.txt', 
            quote=F, row.names=F, col.names=F)



# Task 7
# prep for AMULET run; viz AMULET results
# AMULET output: MultipletBarcodes.txt is an array of barcodes
# translate to barcode_indices and viz marker positions + marker resolutions

# CO resolution for multiplets
prep_res_mult <- function(x,m) {
  barcodes <- read.table(paste0(file_path, '/out/amulet_', x, '/MultipletBarcodes_01.txt'), header=F, sep='\n')
  barcodes <- barcodes$V1
  print(length(barcodes))# number of multiplets should correspond to AMULET report
  
  # create mapping table for barcode_indices and barcodes
  mb <- read.table(paste0(file_path, '/data/', x, '/cellSNP.samples.tsv'), header=F, sep='\n')
  colnames(mb) <- 'barcode'
  mb$barcode_index <- rownames(mb)
  
  # subset mb with wanted barcodes
  barcode_indices <- subset(mb, barcode %in% barcodes)[,2]
  # subset m with wanted barcode_indices, the 2nd col of mb_subset
  #m <- import_cellsnp_output(x)# if m is not given#####
  m_multiplets <- subset(m, barcode_index %in% barcode_indices)
  # now use m_multiplets as m as arg for prep_resolution() function!-----
  
  # or continue below:
  length(unique(m_multiplets$barcode_index))
  # swap position_index with position
  mp <- read.table(paste0(file_path, '/data/', x, '/cellSNP.markers.tsv'), header=F, sep='\n')
  mp <- mp$V1
  m_multiplets$position <- mp[as.numeric(m_multiplets$position_index)]
  m_multiplets$position <- as.integer(m_multiplets$position)
  m_multiplets <- m_multiplets[,-1]
  # sort by barcode_index, then by position
  m_multiplets <- m_multiplets[order(m_multiplets$barcode_index, m_multiplets$position),]#added
  # write to file
  write.table(m_multiplets, file=paste0(file_path, '/data/', x, '/m_multiplets.csv'), sep=',', row.names=F)
  
  return(m_multiplets)
}

# run python code calc_distance()
# after obtaining file "X_res_multiplets.csv"

plot_res_mult <- function(x) {
  # import python output
  res <- read.table(paste0(file_path, '/out/', x, '_res_multiplets.csv'), header=T, sep=',')
  
  # plot
  pdf(paste0(file_path, '/out/', x, '_cellular_multiplets.pdf'))
  hist(res$marker_distance, breaks=50,
       xlab='Nucleotide distances between adjacent markers',
       main='Distribution of CO marker resolution for barcodes\nlabeled as multiplets by AMULET')
  
  hist(log10(res$marker_distance), breaks=50,
       xlab='log10(nucleotide distances between adjacent markers)',
       main='Distribution of CO marker resolution for barcodes\nlabeled as multiplets by AMULET')
  
  # import R output: m_multiplets
  # m_multiplets <- prep_res_mult(x)
  m_multiplets <- read.table(paste0(file_path, '/data/', x, '/m_multiplets.csv'), header=T, sep=',')
  
  # also plot multiplet
  hist(m_multiplets$position, breaks=852, labels=F,
       xlab='Distribution of marker positions',
       main='Marker density for barcodes\nlabeled as multiplets by AMULET')
  
  tm <- table(m_multiplets$barcode_index)
  tm <- data.frame(tm)
  colnames(tm) <- c('barcode_index','occurrence')
  hist(tm$occurrence, labels=F, xlab='Number of mapped markers', main='Histogram of number of markers mapped per barcode')
  
  hist(log10(tm$occurrence),
       xlab='log10(number of mapped markers)',
       main='Histogram of number of markers mapped per barcode')
  
  dev.off()
}


### subset t and m to retain only occurrence>=100 ###
t <- t[t$occurrence>=100,]
m <- subset(m, barcode_index %in% t[t$occurrence>=100,1])
# subsequent analyses will not need occurrence<100


### run function here! ###
file_path <- '/Users/shell/MSc/LMU/Schneeberger'
for (x in c('A','B','C')) {
  m <- import_cellsnp_output(x)
  #plot_marker_coverage_depth(x,m)
  t <- make_count_table(x,m)
  #plot_num_markers_per_barcode(x,t)
  #plot_marker_position(x,m,t)
  # trim t and m (if m is not by default 100+)
  #t <- t[t$occurrence>=100,]
  #m <- subset(m, barcode_index %in% t[t$occurrence>=100,1])
  prep_resolution(x,m)
}

for (x in c('A','B','C')) {
  m <- import_cellsnp_output(x)
  t <- make_count_table(x,m)
  #t <- t[t$occurrence>=100,]
  #m <- subset(m, barcode_index %in% t[t$occurrence>=100,1])
  #prep_resolution(x,m)
  #plot_marker_coverage_depth(x,m)
  #plot_marker_position(x, m, t)
  #plot_num_markers_per_barcode(x,t)
  plot_resolution(x)
}

for (x in c('A','B','C')) {
  #plot_res_mult(x)
  plot_resolution(x)
}

