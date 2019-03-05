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
#### Datasets/GDCdata/* has been ignored for space reasons
```console
BiolinksAnalysis
├── BiolinksAnalysis.Rproj
├── Datasets
│   ├── BLCA.csv
│   ├── BLCA_select.csv
│   ├── ESCA.csv
│   ├── ESCA_select.csv
│   ├── HNSC.csv
│   ├── HNSC_select.csv
│   ├── KICH.csv
│   ├── KICH_select.csv
│   ├── KIRC.csv
│   ├── KIRC_select.csv
│   ├── KIRP.csv
│   ├── KIRP_select.csv
│   ├── LAML.csv
│   ├── LAML_select.csv
│   ├── LIHC.csv
│   ├── LIHC_select.csv
│   ├── LUAD.csv
│   ├── LUAD_select.csv
│   ├── LUSC.csv
│   ├── LUSC_select.csv
│   ├── PAAD.csv
│   ├── PAAD_select.csv
│   ├── STAD.csv
│   └── STAD_select.csv
├── Output
│   ├── Gender
│   │   ├── Graphs
│   │   │   ├── nucChange_Graph
│   │   │   │   ├── muse
│   │   │   │   ├── mutect
│   │   │   │   ├── somaticsniper
│   │   │   │   └── varscan2
│   │   │   ├── TiTv_Graph
│   │   │   │   ├── muse
│   │   │   │   ├── mutect
│   │   │   │   ├── somaticsniper
│   │   │   │   └── varscan2
│   │   │   └── total_mut_Graph
│   │   │       ├── muse
│   │   │       ├── mutect
│   │   │       ├── somaticsniper
│   │   │       └── varscan2
│   │   └── pValues
│   │       ├── nucChange_pVal
│   │       │   ├── muse
│   │       │   ├── mutect
│   │       │   ├── somaticsniper
│   │       │   └── varscan2
│   │       ├── TiTv_pVal
│   │       │   ├── muse
│   │       │   ├── mutect
│   │       │   ├── somaticsniper
│   │       │   └── varscan2
│   │       └── total_mut_pval
│   │           ├── muse
│   │           ├── mutect
│   │           ├── somaticsniper
│   │           └── varscan2
│   ├── Race
│   │   ├── Graphs
│   │   │   ├── nucChange_Graph
│   │   │   │   ├── muse
│   │   │   │   ├── mutect
│   │   │   │   ├── somaticsniper
│   │   │   │   └── varscan2
│   │   │   ├── TiTv_Graph
│   │   │   │   ├── muse
│   │   │   │   ├── mutect
│   │   │   │   ├── somaticsniper
│   │   │   │   └── varscan2
│   │   │   └── total_mut_Graph
│   │   │       ├── muse
│   │   │       ├── mutect
│   │   │       ├── somaticsniper
│   │   │       └── varscan2
│   │   └── pValues
│   │       ├── nucChange_pVal
│   │       │   ├── muse
│   │       │   ├── mutect
│   │       │   ├── somaticsniper
│   │       │   └── varscan2
│   │       ├── TiTv_pVal
│   │       │   ├── muse
│   │       │   ├── mutect
│   │       │   ├── somaticsniper
│   │       │   └── varscan2
│   │       └── total_mut_pval
│   │           ├── muse
│   │           ├── mutect
│   │           ├── somaticsniper
│   │           └── varscan2
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
├── README.md
└── Scripts
    ├── Analysis1.R
    ├── Analysis2.R
    ├── Analysis3.R
    ├── Download.R
    ├── functions1.R
    ├── functions2.R
    ├── functions3.R
    ├── functions.R
    └── select_col.py


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

select_col.py - This a python script that will take the tab-delimated flat files created by Download.R and will extract important columns to make smaller, more efficient files (~/BiolinksAnalysis/Datasets/\*\_select.csv)  


##### Smoking Data
functions1.R - Script that creates functions used in Analysis1.R  

Analysis1.R - Read in the csv data from Download.R and generate pValues/Quartile data for each Cancer, location, and mutational pipeline only comparing smokers and nonsmokers

##### Race Data
functions2.R - Script that creates functions used in Analysis2.R  

Analysis2.R - Read in the csv data from Download.R and generate pValues/Quartile data for each Cancer, location, and mutational pipeline  comparing AfrAmr and EurAmr

##### Gender Data
functions3.R - Script that creates functions used in Analysis3.R  

Analysis3.R - Read in the csv data from Download.R and generate pValues/Quartile data for each Cancer, location, and mutational pipeline  comparing Male and Female

<div id='Datasets'/>  

## Datasets
#### This folder contains all flat files and extra data downloaded due to TCGAbiolinks from Download.R

<div id='Output'/>  

## Output
#### This folder will contain all info generated as a result of Analysis\*.R files such as pValues and graphs

#### The data was calculated for each group, for each cancer, for each valid site in each cancer, and for each somatic varaint pipeline (muse, mutect, somaticsniper, varscan2)  
nucChange_Graphs contains boxplots showing the different nucleotide changes (eg A > G) and their frequency distrubution per person in the population  
nucChange_pVal shows the pValue, calculated by wilcox.text, between each respective population for each nucleotide change

TiTv_Graphs shows the frequency distribution of Transitions and Transversions per person between two populations  
TiTv_pVal shows the pValue between each population for Transitions and Transversions

total_mut_Graphs shows the distribution of the total number of somatic point mutations per person in each population  
total_mut_pVal shows the pValue between each group for the distribution of somatic point mutations  