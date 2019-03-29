## Author: Chris Holt
## Purpose: Looks at mutational burden and age between Smokers and Non-smokers
## Date Created: 28/Mar/2019
## Date of Last Update: 28/Mar/2019

##-----------------------------------------------------------------------------------------------------------
## This file will read in flat files containing gene frequency pvalues for each comparison
## , combines the pipelines and writes out one flat file for each comparison
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
## Read in the data
##-----------------------------------------------------------------------------------------------------------
combine_pvalues <- function(compare, cancers){
  Full_Data <- data.frame()
  pipelines <- c("muse", "mutect", "somaticsniper", "varscan2")
  for (i in 1:length(cancers)){
    for (j in 1:length(pipelines)){
      setwd(paste0("~/Research/BiolinksAnalysis/Output/", compare, "/Genes_Pvalues/"))
      file = Sys.glob(paste0(cancers[i], "_", pipelines[j], '_FINAL_PVALUES.tsv'))
      df = read_delim(file, delim = "\t")
      df$pipeline <- pipelines[j]
      Full_Data <- rbind(Full_Data, df)
    }
    write_delim(Full_Data, file.path("~/Research/BiolinksAnalysis/Output/", paste0(compare, "/Genes_Pvalues/", 
                                                                                   cancers[i], "_pValues_Combined.tsv")), delim = "\t")
  }
  setwd("~/Research/BiolinksAnalysis/")
}
 
Data_Names <- c("BLCA", "HNSC", "KICH", "LUAD", "LUSC", "PAAD", "KIRP")
Data_Names_1 <- c("BLCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "STAD")

combine_pvalues("Smoke", Data_Names) 
combine_pvalues("Race", Data_Names_1)
combine_pvalues("Gender", Data_Names_1)
















