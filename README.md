# BiolinksAnalysis
These files will download, manipulate TCGA data obtained through the R package TCGAbiolinks

## Directory Structure
#### The Datasets and Output directories have been ignored for space reasons. Output may be added later
~/BiolinksAnalysis/Scripts  
~/BiolinksAnalysis/Datasets  
~/BiolinksAnalysis/Output  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/Output/Smoke  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/Smoke/total_mut_graphs  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/Smoke/total_mut_pval    

## Scripts
####This folder contains all R scripts that will be used in this project

functions.R - This will define all  general functions used in this project

functions1.R - Script that creates functions used in Analysis1.R

Download.R - This will download, format, and save all clinical and mutational data to csv files. Files downloaded as "\t" separated flat files and saved to Datasets
directory. 

Analysis1.R - Read in the csv data from Download.R and generate pValues/Quartile data for each Cancer, location, and mutational pipeline only comparing smokers and nonsmokers

select_col.py/select_col.ipynb - This a python script and corresponding jupyter notebook script that will take the tab-delimated flat files created by Download.R and will extract important columns to make smaller, more efficient files. 

## Datasets
####This folder contains all flat files and extra data downloaded due to TCGAbiolinks from Download.R

## Output
####This folder will contain all info generated as a result of Analysis\*.R files such as pValues and graphs

Within total_mut_graphs and total_mut_pval are four directories created from Analysis1.R that separate the output data by the different mutational pipelines
