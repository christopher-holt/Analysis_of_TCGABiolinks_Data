##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

library(tidyverse)
library(TCGAbiolinks)
source("Scripts/functions.R")
##-----------------------------------------------------------------------------------------------------------
## set working directory to datasets directory
##-----------------------------------------------------------------------------------------------------------

setwd('Datasets')

##-----------------------------------------------------------------------------------------------------------
## These are the cancers that have been selected as potentially having a relationship with smoking
## and occurence
##-----------------------------------------------------------------------------------------------------------

Cancers <- c("LAML", "BLCA")
Cancers1 <- c("COAD", "ESCA")
Cancers2 <- c("KICH", "KIRC")
Cancers3 <- c("KIRP", "LUAD")
Cancers4 <- c("LUSC", "PAAD")
Cancers5 <- c("STAD", "LIHC" )

##-----------------------------------------------------------------------------------------------------------
## This will utilise the created merged function from functions.R to download and merge
## this clinical and mutational data from TCGA biolinks for each cancer
##-----------------------------------------------------------------------------------------------------------
for (i in 1:length(Cancers)){
  total <- merged(Cancers[i])
  write.csv(total, paste0(Cancers[i], ".csv"))
  assign(paste(Cancers[i]), total)
}
for (i in 1:length(Cancers1)){
  total <- merged(Cancers1[i])
  write.csv(total, paste0(Cancers1[i], ".csv"))
  assign(paste(Cancers1[i]), total)
}
for (i in 1:length(Cancers2)){
  total <- merged(Cancers2[i])
  write.csv(total, paste0(Cancers2[i], ".csv"))
  assign(paste(Cancers2[i]), total)
}

for (i in 1:length(Cancers3)){
  total <- merged(Cancers3[i])
  write.csv(total, paste0(Cancers3[i], ".csv"))
  assign(paste(Cancers3[i]), total)
}
for (i in 1:length(Cancers4)){
  total <- merged(Cancers4[i])
  write.csv(total, paste0(Cancers4[i], ".csv"))
  assign(paste(Cancers4[i]), total)
}
for (i in 1:length(Cancers5)){
  total <- merged(Cancers5[i])
  write.csv(total, paste0(Cancers5[i], ".csv"))
  assign(paste(Cancers5[i]), total)
}
rm(total)
## -------------------------------------------------------------------------------------------
for (i in 1:length(Cancers)){
  clin <- download_clinical(Cancers[i])
  mut <- download_mutational(Cancers[i])
  mut <- classify_Changes(mut)
  total <- merge_ind(clin,mut)
  write.csv(total, paste0(Cancers[i], ".csv"))
  assign(paste(Cancers[i]), total)
}
for (i in 1:length(Cancers1)){
  clin <- download_clinical(Cancers1[i])
  mut <- download_mutational(Cancers1[i])
  mut <- classify_Changes(mut)
  total <- merge_ind(clin,mut)  write.csv(total, paste0(Cancers1[i], ".csv"))
  assign(paste(Cancers1[i]), total)
}
for (i in 1:length(Cancers2)){
  total <- merged(Cancers2[i])
  write.csv(total, paste0(Cancers2[i], ".csv"))
  assign(paste(Cancers2[i]), total)
}

for (i in 1:length(Cancers3)){
  total <- merged(Cancers3[i])
  write.csv(total, paste0(Cancers3[i], ".csv"))
  assign(paste(Cancers3[i]), total)
}
for (i in 1:length(Cancers4)){
  total <- merged(Cancers4[i])
  write.csv(total, paste0(Cancers4[i], ".csv"))
  assign(paste(Cancers4[i]), total)
}
for (i in 1:length(Cancers5)){
  total <- merged(Cancers5[i])
  write.csv(total, paste0(Cancers5[i], ".csv"))
  assign(paste(Cancers5[i]), total)
}
