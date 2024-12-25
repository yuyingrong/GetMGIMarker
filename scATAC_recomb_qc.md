# Quality control scATAC-seq and cellsnp-lite output

Yuying Rong
20231008

Due to the termination of production of the Chromium Single Cell CNV kit (10x Genomics, 2020), we have identified no suitable alternative single-cell whole genome sequencing (scWGS) protocols for plant anthers.

Therefore, we probed at the possibility of applying single-cell Assay of Transposase Accessible Chromatin sequencing (scATAC-seq) to emulate scWGS. 

Since ATAC-seq preferentially cleaves and sequences nucleosome-free chromatin regions, to increase representation, we adopted the lithium-assisted nucleosome depletion (LAND; Vitak and Torkenczy et al., 2017) method to unwrap the nucleosomes from the chromatin. 

By mapping the scATAC-seq reads back to the reference genome and calling variant at crossover marker sites, we can identify crossover events in the genome of the single cell. 

Full report: https://docs.google.com/document/d/1Cm1XDc0TmMHJJIGBLbFP2rR463NN_ZfrIEGrjmzgFBc/edit?usp=sharing.

## Import cellsnp-lite output file as input

```r
file_path <- '/Users/shell/MSc/LMU/Schneeberger'

# import cellsnp output mmmatrix

import_cellsnp_output <- function (x) {
  # read matrix market file into df
  mtx <- Matrix::readMM(paste0(file_path, '/data/', x, '/cellSNP.tag.DP.mtx'))
  m <- cbind.data.frame(row=mtx@i+1, col=mtx@j+1, x=mtx@x)
  m <- m[,-3]
  colnames(m) <- c('position_index','barcode_index')
  rm(mtx)
  
  # write m to intermediate file
  #write.table(m, 
  #            file=paste0(file_path, '/data/', x, '/m.csv'), 
  #            quote=F, row.names=F, col.names=F)
  return(m)
}
```


## Visualize number of markers per single cell

```r
# make count table t for number of markers covered (position_indices) for each barcode_index

make_count_table <- function (x,m) {
  # count frequency of each barcode_index
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


# plot number of SNPs (markers) per barcode

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
```

## Visualize number of markers per 100 kb windows

```r
# plot marker density: number of markers per region
# plot for each t_group (based on total number of markers per barcode)

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
```

## Calculate nucleotide distance between adjacent markers

Part 1/3. 
R code:

```r
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
```

Part 2/3. 
Python code:

```python
#!/dss/dsshome1/08/ra65mav/.conda/envs/py3.11/bin/python

'''
input file: modified cellsnp output table "m" from R
input arg: file path
input arg: sample

Output: csv file with barcode and distances between adjacent markers for each barcode
csv format: barcode, distance
'''

import sys

import pandas as pd
from collections import defaultdict


# sys.argv[0] is the name of the python script
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

#input_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/data/A/m.csv'
#output_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/out/res.txt'


# convert cellsnp-lite output to intermediate dict
def table_to_interval_dict(input_file_path, output_file_path):
    # read input file to table
    t = pd.read_table(input_file_path, sep=',')
    # first order table by barcode, then by position
    t = t.sort_values(by=['barcode_index', 'position'])
    # reset indices so enumerate will start from 0, not original indices
    t = t.reset_index(drop=True)
    # open output file
    with open(output_file_path, 'w+') as fw:
        # write header
        fw.write('barcode_index,marker_distance\n')
        for i,barcode_index in enumerate(t['barcode_index']):
            try:
                if t['barcode_index'][i+1]==barcode_index:
                    # calc distance for each pair of positions within the same barcode
                    distance = t['position'][i+1]-t['position'][i]
                    fw.write(f'{barcode_index},{distance}\n')
            except KeyError:
                # when t['barcode_index'][i+1] iterates out of range
                print('done')

    return 0


table_to_interval_dict(input_file_path, output_file_path)
```

Part 3/3. 
R code:

```r
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
```


## Detect cellular multiplets based on marker coverage depth.

Cellular multiplet: during single-cell sequencing library preparation, >1 cells may enter one emulsion bead, resulting in one barcode representing the genome of multiple cells. 

Each marker should only be covered once per single cell, because a single pollen cell is haploid. A marker coverage depth >1 represents a multiplet. 

```r
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
```


## Run AMULET for cellular multiplet detection.

AMULET is designed for diploid genomes, not suitable for our haploid data. Nevertheless we tried it. 

Part 1/3. 
R code:

```r
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
```

Part 2/3. 
Python code:

```python
#!/dss/dsshome1/08/ra65mav/.conda/envs/py3.11/bin/python

'''
input file: modified cellsnp output table "m" from R
input arg: file path
input arg: sample

Output: csv file with barcode and distances between adjacent markers for each barcode
csv format: barcode, distance
'''

import sys
import csv


# sys.argv[0] is the name of the python script
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

#input_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/data/A/m.csv'
#output_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/out/res.txt'


def calc_distance(input_file_path, output_file_path):
    # open output file
    with open(output_file_path, 'w+') as fw:
        # write header
        fw.write('barcode_index,marker_distance\n')
        
        # initiate barcode_index
        prev_barcode_index = ''
        prev_position = ''
        counter = 0
        # input file contains only valid barcodes, and is already ordered by barcode then by position
        with open(input_file_path) as f:
            # skip header
            next(f)
            # read one row at a time; instead of reading in the whole file
            rows = csv.reader(f, delimiter=',', quotechar='"')
            
            for row in rows:
                barcode_index = row[0]
                position = int(row[1])
                if barcode_index == prev_barcode_index:
                    distance = position - prev_position
                    fw.write(f'{barcode_index},{distance}\n')
                    prev_position = position
                elif barcode_index != prev_barcode_index:
                    # if the first row is being iterated
                    # or if all rows for prev_barcode_index has been iterated
                    prev_barcode_index = barcode_index
                    prev_position = position
                    counter += 1
                    continue
            print(f'Iterated {counter} barcodes.')
        print('done')

    return 0


calc_distance(input_file_path, output_file_path)
```

Part 3/3: visualize AMULET output. 
R code:

```r
# after running python calc_distance()
# and obtaining "X_res_multiplets.csv"

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
```


## Detect barcode multiplets. 

Barcode multiplet: during single-cell sequencing library preparation, if an emulsion bead received >1 barcodes, then multiple barcodes should show the same marker pattern. 

I decided not to pursue this detection because single-cell data is already sparse and multiple true single cells may show the same marker pattern by chance alone.  

```python
#!/dss/dsshome1/08/ra65mav/.conda/envs/py3.11/bin/python
'''
This script generates avg depth with and without removal of barcode multiplets
# input: unmodified cellsnp output
# part 1: create defaultdict with d[barcode_index]((position, count)); output the dict
# part 2: reduce dict with duplicate vals (barcode multiplets)
# part 3: write avg depth per barcode_index, write to file
# to contain only the keys with unique values; output the cleaned dict
'''

import sys
from collections import defaultdict


# sys.argv[0] is the name of the python script
cellsnp_file_path = sys.argv[1]
output1_file_path = sys.argv[2]
output2_file_path = sys.argv[3]

# Example file paths
# cellsnp_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/cellsnp_lite_results/B_nigra_paternal_5733_B/cellSNP.tag.DP.mtx'
# output1_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/out/B_avgs.csv'
# output2_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/out/B_avgs_dup_rm.csv'

print(f"cellsnp_file_path = '{cellsnp_file_path}';")
print(f"output1_file_path (all cells) = '{output1_file_path}';")
print(f"output2_file_path (after barcode multiplet removal) = '{output2_file_path}'.")



def create_dict(in_file):
    # create dict of sets for each barcode_id
    d = defaultdict(set)
    
    with open(in_file, 'r') as f:
        # skip first three rows
        next(f)
        next(f)
        next(f)
        # format of in_file: position_id \t barcode_id \t count \n
        
        for line in f:
            line = line.strip('\n').split('\t')
            # add (position_id,count) to each barcode_id
            # using a set means no (position_id,count) duplicates will be added to the same barcode
            d[line[1]].add((line[0], int(line[2])))
    
    f.close()
    print(f'Number of cells before dup and zero removal: {len(d)}')
    
    return d


# remove duplicates based on 1) marker positions and 2) per marker position mapping counts
# if barcode_indices have identical 1) and 2) then they are considered barcode multiplets and removed

def remove_duplicates_and_zeros(d):
    # empty duplicates and retain unique
    count_dup = 0#
    for i,value in enumerate(d.values()):
        # search for (position_id,count) duplicates between barcodes
        # if counted more than one pair between barcodes
        if list(d.values()).count(value) > 1:
            # update counter
            count_dup += 1
            # remove the content of the current barcode and preserve the content of the other barcode(s)
            # this will remove the content of all barcodes with duplicate contents, until one remain
            d[list(d.keys())[i]] = ()
    
    # how many cells had multiple barcodes?
    print(f'Number of duplicate keys: {count_dup}')#
    
    # create tuple to contain barcodes with no content
    # either due to no content or removal in prev step
    zero_keys = ()
    for i,value in enumerate(d.values()):
        # if a key's value is empty, add key to the tuple
        # the only way to grab key from evaluating value is by enumerating the keys/values
        if not value:
            zero_keys += (list(d.keys())[i],)
            
    print(f'Number of zero keys after emptying duplicates: {len(zero_keys)}')#
    
    # rm zero_keys from output dict
    for key in zero_keys:
        d.pop(key)
    
    print(f'Number of cells after dup and zero removal: {len(d)}')#
    
    return d


# write output to file
d = create_dict(cellsnp_file_path)
avgs = [(int(key), sum([tup[1] for tup in d[key]])/len(d[key])) for key in d]
with open(output1_file_path, 'w+') as fw:
    fw.write('position_index,avg_depth\n')
    for pair in avgs:
        fw.write(f'{pair[0]},{round(pair[1],2)}\n')
print('done 1/2')

# (optional)
d = remove_duplicates_and_zeros(d)
avgs = [(int(key), sum([tup[1] for tup in d[key]])/len(d[key])) for key in d]
with open('/Users/shell/MSc/LMU/Schneeberger/out/{}_avgs1.csv', 'w+') as fw:
    fw.write('position_index,avg_depth\n')
    for pair in avgs:
        fw.write(f'{pair[0]},{round(pair[1],2)}\n')
print('done 2/2')
```


## Calculate breadth and depth of coverage by fragments.

A fragment is represented as the region between the start and end of each pair of reads.

Part 1/2. 
R code:

```r
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
```

Part 2/2. 
Python code:

```python

#!/dss/dsshome1/08/ra65mav/.conda/envs/py3.11/bin/python

'''
Input 1: cellranger output fragments.tsv, after
- subsetting for wanted barcodes, which also removes headers starting with '#'
- removing unwanted cols
- sorting rows by barcodes

Output: txt tile containing max_depth, mean_depth, median_depth, and breadth
'''

import sys

import csv
import numpy as np

# sys.argv[0] is the name of the python script
#barcodes_file_path = sys.argv[1]
fragments_file_path = sys.argv[1]
stats_file_path = sys.argv[2]

# Example file paths
# barcodes_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/tb1000_barcodes.txt'
# fragments_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/data/A/frags1000.tsv'
# stats_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/out/stats_tb1000.csv'

#print(f"barcodes_file_path = '{barcodes_file_path}';")
print(f"fragments_file_path = '{fragments_file_path}';")
print(f"stats_file_path = '{stats_file_path}'.")



def calc_stats_per_barcode(fragments_file_path, stats_file_path):
    # open output file for writing
    fw = open(stats_file_path, 'w+')
    fw.write('barcode,max,mean,median,breadth\n')
    counter = 0
    # initiate barcode
    barcode = 'barcode'
    # input file contains only valid barcodes and is already ordered by barcode
    with open(fragments_file_path) as f:
        # read one row at a time; instead of reading in the whole file
        rows = csv.reader(f, delimiter="\t", quotechar="'")

        for row in rows:
            # format of a row in frags.tsv: 7\t43\tGGGTGTCAGCGCATTT-1\n
            
            # the start/end indices are not nucleotide positions, but (nucleotide position - 1)
            start = int(row[0])-4-1
            end = int(row[1])+5-1
            # convert start/end indices into readpair indicies
            read1 = np.arange(start, min(start+50, end), dtype='int32')
            read2 = np.arange(max(start, end-49), end, dtype='int32')
            
            current_barcode = row[2]
            
            if barcode == current_barcode:
                # update the arr
                arr[np.union1d(read1, read2)]+=1
            elif barcode != current_barcode:
                # this means that either all frags for this barcode have been iterated, since rows are ordered by barcode
                # or the first row is being iterated
                try: 
                    # write stats for this barcode to output file
                    fw.write(f'{barcode},{np.max(arr)},{round(np.average(arr[arr!=0]),2)},{np.median(arr[arr!=0])},{arr[arr!=0].shape[0]}\n')
                    counter += 1
                except NameError:
                    # this should only happen when the first row is iterated, so arr has not been declared
                    print(barcode)### should print 'barcode' ###
                    pass
                # update barcode
                barcode = current_barcode
                # (re-)initiate arr of all genomic positions
                arr = np.zeros(shape=(85145900,), dtype='int32')
                # update the arr
                arr[np.union1d(read1, read2)]+=1
        # write the last line to file
        fw.write(f'{barcode},{np.max(arr)},{round(np.average(arr[arr!=0]),2)},{np.median(arr[arr!=0])},{arr[arr!=0].shape[0]}\n')
        counter += 1
        fw.close()
        print(f'Wrote {counter} lines to file.')
    
    return 0

calc_stats_per_barcode(fragments_file_path, stats_file_path)
```

