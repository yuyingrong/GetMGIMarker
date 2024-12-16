This repository contains some scripts I wrote for analyses and efficiency. 

# List of scripts
Tools (details see bottom)
- GetGeneID.py
- GetMGIMarker.py

Single-cell RNA-seq
- 20220315qc.R (quality control visualizations)
- 20220325merge.R (merge datasets sequenced from multiple samples)
- 20220407normalize.R
- 20220413clustering.R
- 20220425trajectory.R
- 20220518trajectory_function.R (clean code organized into functions)
- 20220512dispersion_on_server.R
- 20220614loom_viz.R
- 20220701loom_heatmap.R

Bulk RNA-seq
- 20220625DESeq2.R
- 20220628GO.KEGG.R
- 20220720loom_viz_func.R

Single-cell same-cell RNA/ATAC-seq integration
- 20230206seurat.R
- 20230206signac.R
- 20230208archR.R
- 20230224doublet.R
- 20230316portal.ipynb (see: https://github.com/YangLabHKUST/Portal/issues/6)
- 20230223scrublet.ipynb

Nucleosome-removed single-cell ATAC-seq data analysis
- 20231008clean_code.R (functions and guide to the analysis)
- 20230731calc_resolution.py
- 20230810calc_depth_breadth_per_barcode.py
- 20230820count_barcode_duplicates.py
- 20230923calc_distance.py
- 20230928calc_avg_mkr_cvg_dep.py (average marker coverage depth, number of times an SNP marker is covered by reads)

Brassica rapa trio-binning genome assembly, recombination landscape
- 20240722extract_seq_append_to_fasta.py
- 20240722rm_seq_from_fasta.py

(below not uploaded)
- cellranger.sh
- jellyfish.sh
- busco.sh
- cellsnp-lite.sh
- hifiasm.sh
- merqury.sh
- minimap2.sh
- mmseqs2.sh
- mosdepth.sh
- mummer.sh
- omni-c.sh
- purge_dups.sh
- quast.sh
- ragtag.sh
- syri.sh
- yak.sh

Misc
- 20220315ms.R: mass spectrometry data clustering
- 2018-2019: plant duplicate gene retention scripts: https://drive.google.com/drive/folders/1vftU-BLgMZeUw42-xUxTfKDPOa_Gy6qy?usp=sharing

# Tools

**GetGeneID.py**
Last update: 2022 May 09.
Required package: selenium.
Required installation: chromedriver. 
Written with Python 3.10.0.

Purpose:
Bulk convert human marker IDs to Rhesus monkey gene ID and Ensembl ID from NCBI Gene (https://www.ncbi.nlm.nih.gov/gene).

Usage:
- Enter the path where you installed the chromedriver.
- Replace the human marker IDs in the input list `l`.
- Replace the 'rhesus' if you are searching for another organism (organism must be present in the NCBI Gene database).
- Run the script.
- Results will be printed to the console. 

**GetMGIMarker.py**
Last update: 2020 Dec 13 (Script may be outdated due to selenium update, see GetGeneID for updated usage).
Required package: selenium.
Written with Python 3.8.3.

Purpose:
Bulk convert protein names to marker IDs for Seurat input.
Scrapes info from: http://www.informatics.jax.org/batch/summary.

Usage:
- Use text editor to replace the list of mRNA/protein names with your search terms and run the script.
- Results will be returned to the console. 




