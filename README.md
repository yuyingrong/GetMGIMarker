This repository contains scripts I wrote for efficiency. Feel free to adopt them for your own use. 

# GetGeneID.py
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

# GetMGIMarker.py
Last update: 2020 Dec 13 (Script may be outdated due to selenium update, see GetGeneID for updated usage).
Required package: selenium.
Written with Python 3.8.3.

Purpose:
Bulk convert protein names to marker IDs for Seurat input.
Scrapes info from: http://www.informatics.jax.org/batch/summary.

Usage:
- Use text editor to replace the list of mRNA/protein names with your search terms and run the script.
- Results will be returned to the console. 
