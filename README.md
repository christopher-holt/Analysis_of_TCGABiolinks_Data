# BiolinksAnalysis
These files will download, manipulate TCGA data obtained through the R package TCGAbiolinks

#### Scripts
This folder contains all R scripts that will be used in this project

functions.R - This will define all functions used in this project

Download.R - This will download, format, and save all clinical and mutational data to csv files. Files downloaded as "\t" separated flat files and saved to a different
folder. 

Analysis1.R - Read in the csv data from Download.R and remove any sites that do not contain enough data to do an analysis
