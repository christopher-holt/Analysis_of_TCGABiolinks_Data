# BiolinksAnalysis
These scripts will download, manipulate, and analyse TCGA data obtained through the R package TCGAbiolinks

## Table of Contents
* [Languages](#languages)
* [Packages](#packages)  
* [Scripts](#Scripts)
* [Datasets](#Datasets)
* [Output](#Output)



<div id='languages'/>  

## Languages
```console
$ python3 --version
Python 3.6.7

$ R --version
R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under the terms of the
GNU General Public License versions 2 or 3.
For more information about these matters see
http://www.gnu.org/licenses/.
```
<div id='packages'/>  

## Packages  
#### These packages need to be installed (R)
```r
> install.packages(c("tidyverse", "BiocManager"))
> BiocManager::install("TCGAbiolinks")
```
#### These packages need to be installed (python3) 
```bash
$ pip3 install --user pandas
$ pip3 install --user scipy
$ pip3 install --user os
$ pip3 install --user glob
```
<div id='Scripts'/>  

## Scripts
#### This folder contains all R scripts that will be used in this project

##### General scripts

functions.R - This will define all  general functions used in this project  

Download.R - This will download, format, and save all clinical and mutational data to csv files. Files downloaded as "\t" separated flat files and saved to Datasets
directory  

select_col.py - This a python script that will take the tab-delimated flat files created by Download.R and will extract important columns to make smaller, more efficient files (~/BiolinksAnalysis/Datasets/\*\_select.csv)  

MannU_Test.py - This python file will calculate the pValues for all flat files written out by Analysis*_2.R. For processing reasons, you must specify the individual comparison directory (Smoke/Race/Gender) and then a valid cancer cohort within that directory (more info in file)  

combine_pValues.R - This file will take the output of MannU_Test.py and combine them into one flat file per cancer to make later analysis easier  

significant_genes.R - This will file will read in the outputs of combine_pValues.R and create two flat files. One has all the genes for each pipeline with a pvalue less than 0.05 and the other is a summary of the number of significant genes for each pipeline and cancer

##### Smoking Data
functions1.R - Script that creates functions used in Analysis1.R  

Analysis1.R - Read in the csv data from Download.R and generate pValues/Quartile data for each Cancer, location, and mutational pipeline only comparing smokers and nonsmokers  

Analysis1_1.R - This will perform the same action as Analysis1.R except that it will look at the whole cancer, not specific sites

Analysis1_2.R - This script will calculate the frequencies of Genes in each population and write out a tab delim file

Analysis1_3.R - This script will perform the same analysis as Analysis1_1.R except that it will focus only on stage I cancers

Analysis1_4.R - This script will calculate the number of mutations for each age and look for a relation

##### Race Data
functions2.R - Script that creates functions used in Analysis2.R  

Analysis2.R - Read in the csv data from Download.R and generate pValues/Quartile data for each Cancer, location, and mutational pipeline  comparing AfrAmr and EurAmr  

Analysis2_1.R - This will perform the same action as Analysis2.R except that it will look at the whole cancer, not specific sites\

Analysis2_2.R - This script will calculate the frequencies of Genes in each population and write out a tab delim file

Analysis2_4.R - This script will calculate the number of mutations for each age and look for a relation

##### Gender Data
functions3.R - Script that creates functions used in Analysis3.R  

Analysis3.R - Read in the csv data from Download.R and generate pValues/Quartile data for each Cancer, location, and mutational pipeline  comparing Male and Female  

Analysis3_1.R - This will perform the same action as Analysis3.R except that it will look at the whole cancer, not specific sites

Analysis3_2.R - This script will calculate the frequencies of Genes in each population and write out a tab delim file

Analysis3_4.R - This script will calculate the number of mutations for each age and look for a relation

<div id='Datasets'/>  

## Datasets
#### This folder contains all flat files and extra data downloaded due to TCGAbiolinks from Download.R

<div id='Output'/>  

## Output
#### This folder will contain all info generated as a result of R/py files such as pValues and graphs and frequency tables 

#### The data was calculated for each group, for each cancer, for each valid site in each cancer, and for each somatic varaint pipeline (muse, mutect, somaticsniper, varscan2)  
Graphs/nucChange_Graphs contains boxplots showing the different nucleotide changes (eg A > G) and their frequency distrubution per person in the population  
pValues/nucChange_pVal shows the pValue, calculated by wilcox.text, between each respective population for each nucleotide change

Graphs/TiTv_Graphs shows the frequency distribution of Transitions and Transversions per person between two populations  
pValues/TiTv_pVal shows the pValue between each population for Transitions and Transversions

Graphs/total_mut_Graphs shows the distribution of the total number of somatic point mutations per person in each population  
pValues/total_mut_pVal shows the pValue between each group for the distribution of somatic point mutations  

Age contains graphs showing the relation between number of mutations and age.   
Age/Files contains flat files with the rsquared values of the linear regression line  
Age/Graphs contains the .jpg files of the graphs

Genes_Pvalues contains gene frequencies flat files and the files containg pvalues (\*_FINAL_PVALUES) and the combined files of FINAL_PVALUES (\*_pValues_Combined)  
Genes_Pvalues/Summary contains flat files with all the genes that have pValues less than 0.05 and a combined summary table showing the number of significant genes for each pipeline  