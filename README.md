# ST_Human_SI

Code and functions to process spatial transcriptomic (10X Visium) data (RAW) of the human small intestine. Run the main script to produce the processed data that was used to support the results in the manuscript "A spatial expression atlas of *he adult human small intestine". To generate figues run script in the "scr" directory.

### Processed data includes:
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
* run time 3-5 min.
* Data will automatically by saved to the "Data/Processsed" directory.
* To procude figues in the manuscirpt run: figues.m       

## Notes

* Crypt-villus gene expression zonation profiles can be browsed via our web-app: https://itzkovitzwebapps.weizmann.ac.il/webapps/home/session.html?app=Human_villus_zonation1_1
* processed zonation tables are provided as supplementary tables in the manuscript. 
