##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

library(tidyverse)
library(TCGAbiolinks)
source(paste0(getwd(),"/../functions.R"))

##-----------------------------------------------------------------------------------------------------------
## These are the cancers that have been selected as potentially having a relationship with smoking
## and occurence
##-----------------------------------------------------------------------------------------------------------

Cancers <- c("LAML", "BLCA", "COAD", "ESCA","KICH", "KIRC", "KIRP", "LUAD", "LUSC", "PAAD", "STAD", "LIHC" )
Cancers <- "HNSC"
##-----------------------------------------------------------------------------------------------------------
## This will utilise the created merged function from functions.R to download and merge
## this clinical and mutational data from TCGA biolinks for each cancer
##-----------------------------------------------------------------------------------------------------------

for (i in 1:length(Cancers)){
  total <- merged(Cancers[i])
  assign(paste(Cancers[i], "1", sep = "_"), total)
}


