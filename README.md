# ST_Human_SI

Code and functions to process spatial transcriptomic (10X Visium) data (RAW) of the human small intestine. Run the main script to produce the processed data that was used to support the results in the manuscript "A spatial expression atlas of *he adult human small intestine".

# Processed data includes:
- RAW and Normalized Visium Data.
- Patient metadata.
- Spatial metadata.
- Crypt-Villus + Stromal-Epithlial zonation tables.
- scRNAseq-based cell type signature tables.

## Prerequisites

* Matlab R2023a 

## Installation

* Download RAW and intermediary input data files: [DOI 10.5281/zenodo.8396559](https://zenodo.org/record/8396559/files/Data.zip?download=1) 
* Unpack the zip to the source directory of the main script (...\main\).
* Run the main script: import_and_process_visium_data.m      

## Notes

* Crypt-villus gene expression zonation profiles can be browsed via our web-app: https://itzkovitzwebapps.weizmann.ac.il/webapps/home/session.html?app=Human_villus_zonation1_1
* processed zonation tables are provided as supplementary tables in the manuscript. 
* Code to generate intermediary results can be provided upon request.
