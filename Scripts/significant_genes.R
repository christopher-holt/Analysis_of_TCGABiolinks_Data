## Author: Chris Holt
## Purpose: Analyses pValues 
## Date Created: 28/Mar/2019
## Date of Last Update: 28/Mar/2019

##-----------------------------------------------------------------------------------------------------------
## This file will read in flat files outputed by combine_pValues.R
##-----------------------------------------------------------------------------------------------------------

##-----------------------------------------------------------------------------------------------------------
## Clear Global Env
##-----------------------------------------------------------------------------------------------------------

rm(list = ls())

##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------
library(tidyverse)


##-----------------------------------------------------------------------------------------------------------
## Defining main function
##-----------------------------------------------------------------------------------------------------------

count_sig_genes <- function(compare, cancers){
  Main_Data <- tribble(~Pipeline, ~num_sig_genes, ~Cancer)
  for (i in 1:length(cancers)){
    setwd(paste0("~/Research/BiolinksAnalysis/Output/", compare, "/Genes_Pvalues/"))
    df <- read_delim(paste0(cancers[i], "_pValues_Combined.tsv"), delim = "\t")
    
    pipelines <- unique(df$pipeline)
    sig_genes <- df %>% filter(pvalue < 0.05)
    write_delim(sig_genes, file.path(paste0("~/Research/BiolinksAnalysis/Output/", compare, "/Genes_Pvalues/Summary/", 
                                            cancers[i], "_sig_genes.tsv")), delim = "\t")
    Data <- tribble(~Pipeline, ~num_sig_genes, ~Cancer)
    for (j in 1:length(pipelines)){
      df1 <- df %>% filter(pipeline == pipelines[j], pvalue < 0.05)
      temp <- tribble(~Pipeline, ~num_sig_genes, ~Cancer,
                      paste0(pipelines[j]), nrow(df1), cancers[i])
      
      Data <- rbind(Data, temp)
    }
    Main_Data <- rbind(Main_Data, Data)
  }
  write_delim(Main_Data, file.path(paste0("~/Research/BiolinksAnalysis/Output/", compare, "/Genes_Pvalues/Summary/", 
                                          compare, "_sig_genes_summary.tsv")), delim = "\t")
  setwd("~/Research/BiolinksAnalysis/")
}


##-----------------------------------------------------------------------------------------------------------
## Defining parametres and running main
##-----------------------------------------------------------------------------------------------------------
Data_Names <- c("BLCA", "HNSC", "KICH", "LUAD", "LUSC", "PAAD", "KIRP")
Data_Names_1 <- c("BLCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "STAD")


  
count_sig_genes("Smoke", Data_Names)  
count_sig_genes("Race", Data_Names_1) 
count_sig_genes("Gender", Data_Names_1)  


  
















