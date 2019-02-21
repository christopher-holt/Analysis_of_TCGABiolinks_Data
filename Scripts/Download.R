##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

source("Scripts/functions.R")

##-----------------------------------------------------------------------------------------------------------
## set working directory to datasets directory
##-----------------------------------------------------------------------------------------------------------

setwd('Datasets')

##-----------------------------------------------------------------------------------------------------------
## These are the cancers that have been selected as potentially having a relationship with smoking
## and occurence
##-----------------------------------------------------------------------------------------------------------

Cancers <- c("LAML", "BLCA", "ESCA","KICH", "KIRC","KIRP", "LUAD","LUSC", "PAAD","STAD", "LIHC" )

##-----------------------------------------------------------------------------------------------------------
## This will utilise the created merged function from functions.R to download and merge
## this clinical and mutational data from TCGA biolinks for each cancer
##-----------------------------------------------------------------------------------------------------------

for (i in 1:length(Cancers)){
  clin <- download_clinical(Cancers[i])
  mut <- download_mutational(Cancers[i])
  mut <- classify_Changes(mut)
  total <- merge_ind(clin,mut)
  write.csv(total, paste0(Cancers[i], ".csv"))
  assign(paste(Cancers[i]), total)
}