# BiolinksAnalysis
These scripts will download, manipulate, and analyse TCGA data obtained through the R package TCGAbiolinks

## Table of Contents
* [Directory Structure](#directory_structure)
* [Languages](#languages)
* [Packages](#packages)  
* [Scripts](#Scripts)
* [Datasets](#Datasets)
* [Output](#Output)

<div id='directory_structure'/>  

## Directory Structure
#### Datasets and Ouput/Smoke/Graphs have been ignored for space reasons
```console
BiolinksAnalysis
├── Datasets
├── Output
│   └── Smoke
│       ├── Graphs
│       │   ├── nucChange_Graph
│       │   │   ├── muse
│       │   │   ├── mutect
│       │   │   ├── somaticsniper
│       │   │   └── varscan2
│       │   ├── TiTv_Graph
│       │   │   ├── muse
│       │   │   ├── mutect
│       │   │   ├── somaticsniper
│       │   │   └── varscan2
│       │   └── total_mut_Graph
│       │       ├── muse
│       │       ├── mutect
│       │       ├── somaticsniper
│       │       └── varscan2
│       └── pValues
│           ├── nucChange_pVal
│           │   ├── muse
│           │   ├── mutect
│           │   ├── somaticsniper
│           │   └── varscan2
│           ├── TiTv_pVal
│           │   ├── muse
│           │   ├── mutect
│           │   ├── somaticsniper
│           │   └── varscan2
│           └── total_mut_pval
│               ├── muse
│               ├── mutect
│               ├── somaticsniper
│               └── varscan2
└── Scripts

```
 

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
> install.packages(c("tidyverse", "BiocManager", "magrittr"))
> BiocManager::install("TCGAbiolinks")
```
#### These packages need to be installed (python3) 
```console
$ pip3 install --user pandas
```
<div id='Scripts'/>  

## Scripts
#### This folder contains all R scripts that will be used in this project

##### General scripts

functions.R - This will define all  general functions used in this project  

Download.R - This will download, format, and save all clinical and mutational data to csv files. Files downloaded as "\t" separated flat files and saved to Datasets
directory  

select_col.py - This a python script that will take the tab-delimated flat files created by Download.R and will extract important columns to make smaller, more efficient files  


##### Smoking Data
functions1.R - Script that creates functions used in Analysis1.R  

Analysis1.R - Read in the csv data from Download.R and generate pValues/Quartile data for each Cancer, location, and mutational pipeline only comparing smokers and nonsmokers


<div id='Datasets'/>  

## Datasets
#### This folder contains all flat files and extra data downloaded due to TCGAbiolinks from Download.R

<div id='Output'/>  

## Output
#### This folder will contain all info generated as a result of Analysis\*.R files such as pValues and graphs