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

Cancers <- c("LAML", "BLCA")
Cancers1 <- c("ESCA","KICH")
Cancers2 <- c("KIRC","KIRP")
Cancers3 <- c("LUAD","LUSC")
Cancers4 <- c("PAAD","STAD")
Cancers5 <- c("LIHC" )

##-----------------------------------------------------------------------------------------------------------
## This will utilise the created merged function from functions.R to download and merge
## this clinical and mutational data from TCGA biolinks for each cancer
##-----------------------------------------------------------------------------------------------------------

for (i in 1:length(Cancers)){
  clin <- download_clinical(Cancers[i])
  mut <- download_mutational(Cancers[i])
  mut <- classify_Changes(mut)
  total <- merge_ind(clin,mut)
  write_delim(total, paste0(Cancers[i], ".csv"), delim = "\t")
  assign(paste(Cancers[i]), total)
}
for (i in 1:length(Cancers1)){
  clin <- download_clinical(Cancers1[i])
  mut <- download_mutational(Cancers1[i])
  mut <- classify_Changes(mut)
  total <- merge_ind(clin,mut)
  write_delim(total, paste0(Cancers1[i], ".csv"), delim = "\t")
  assign(paste(Cancers1[i]), total)
}
for (i in 1:length(Cancers2)){
  clin <- download_clinical(Cancers2[i])
  mut <- download_mutational(Cancers2[i])
  mut <- classify_Changes(mut)
  total <- merge_ind(clin,mut)
  write_delim(total, paste0(Cancers2[i], ".csv"), delim = "\t")
  assign(paste(Cancers2[i]), total)
}
for (i in 1:length(Cancers3)){
  clin <- download_clinical(Cancers3[i])
  mut <- download_mutational(Cancers3[i])
  mut <- classify_Changes(mut)
  total <- merge_ind(clin,mut)
  write_delim(total, paste0(Cancers3[i], ".csv"), delim = "\t")
  assign(paste(Cancers3[i]), total)
}
for (i in 1:length(Cancers4)){
  clin <- download_clinical(Cancers4[i])
  mut <- download_mutational(Cancers4[i])
  mut <- classify_Changes(mut)
  total <- merge_ind(clin,mut)
  write_delim(total, paste0(Cancers4[i], ".csv"), delim = "\t")
  assign(paste(Cancers4[i]), total)
}
for (i in 1:length(Cancers5)){
  clin <- download_clinical(Cancers5[i])
  mut <- download_mutational(Cancers5[i])
  mut <- classify_Changes(mut)
  total <- merge_ind(clin,mut)
  write_delim(total, paste0(Cancers5[i], ".csv"), delim = "\t")
  assign(paste(Cancers5[i]), total)
}

rm(clin)
rm(mut)
rm(total)



