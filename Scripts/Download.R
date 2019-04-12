## Author: Chris Holt
## Purpose: To download and standardise TCGAbiolinks data
## Date Created: Feb/2019
## Date of Last Update: 26/Mar/2019

##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

source("Scripts/functions.R")

##-----------------------------------------------------------------------------------------------------------
## set working directory to datasets directory
##-----------------------------------------------------------------------------------------------------------

setwd('Datasets')
columnsofInterest <- c("tumor_barcode", "primary_diagnosis", "tumor_stage", "gender", "ethnicity", "race", "age_at_diagnosis", 
                       "pipeline", "Mut_Type", "nucChange", "hugo_symbol", "disease", "years_smoked", "cigarettes_per_day",
                       "chromosome", "consequence", "one_consequence", "gene", "exon_number", "tissue_or_organ_of_origin",
                       "variant_classification", "mutation_status")
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
  total <- inner_join(clin,mut, by = "tumor_barcode")
  total_sub <- total %>% select(one_of(columnsofInterest))
  write_delim(total_sub, paste0(Cancers[i], "_select.csv"), delim = "\t")
  #write_delim(total, paste0(Cancers[i], ".csv"), delim = "\t")
  assign(paste(Cancers[i]), total)
  assign(paste(Cancers[i], "_select"), total)
}
for (i in 1:length(Cancers1)){
  clin <- download_clinical(Cancers1[i])
  mut <- download_mutational(Cancers1[i])
  mut <- classify_Changes(mut)
  total <- inner_join(clin,mut, by = "tumor_barcode")
  total_sub <- total %>% select(one_of(columnsofInterest))
  write_delim(total_sub, paste0(Cancers1[i], "_select.csv"), delim = "\t")
  #write_delim(total, paste0(Cancers1[i], ".csv"), delim = "\t")
  assign(paste(Cancers1[i]), total)
  assign(paste(Cancers1[i], "_select"), total)
  
}
for (i in 1:length(Cancers2)){
  clin <- download_clinical(Cancers2[i])
  mut <- download_mutational(Cancers2[i])
  mut <- classify_Changes(mut)
  total <- inner_join(clin,mut, by = "tumor_barcode")
  total_sub <- total %>% select(one_of(columnsofInterest))
  write_delim(total_sub, paste0(Cancers2[i], "_select.csv"), delim = "\t")
  #write_delim(total, paste0(Cancers2[i], ".csv"), delim = "\t")
  assign(paste(Cancers2[i]), total)
  assign(paste(Cancers2[i], "_select"), total)
  
}
for (i in 1:length(Cancers3)){
  clin <- download_clinical(Cancers3[i])
  mut <- download_mutational(Cancers3[i])
  mut <- classify_Changes(mut)
  total <- inner_join(clin,mut, by = "tumor_barcode")
  total_sub <- total %>% select(one_of(columnsofInterest))
  write_delim(total_sub, paste0(Cancers3[i], "_select.csv"), delim = "\t")
  #write_delim(total, paste0(Cancers3[i], ".csv"), delim = "\t")
  assign(paste(Cancers3[i]), total)
  assign(paste(Cancers3[i], "_select"), total)
  
}
for (i in 1:length(Cancers4)){
  clin <- download_clinical(Cancers4[i])
  mut <- download_mutational(Cancers4[i])
  mut <- classify_Changes(mut)
  total <- inner_join(clin,mut, by = "tumor_barcode")
  total_sub <- total %>% select(one_of(columnsofInterest))
  write_delim(total_sub, paste0(Cancers4[i], "_select.csv"), delim = "\t")
  #write_delim(total, paste0(Cancers4[i], ".csv"), delim = "\t")
  assign(paste(Cancers4[i]), total)
  assign(paste(Cancers4[i], "_select"), total)
  
}
for (i in 1:length(Cancers5)){
  clin <- download_clinical(Cancers5[i])
  mut <- download_mutational(Cancers5[i])
  mut <- classify_Changes(mut)
  total <- inner_join(clin,mut, by = "tumor_barcode")
  total_sub <- total %>% select(one_of(columnsofInterest))
  write_delim(total_sub, paste0(Cancers5[i], "_select.csv"), delim = "\t")
 # write_delim(total, paste0(Cancers5[i], ".csv"), delim = "\t")
  assign(paste(Cancers5[i]), total)
  assign(paste(Cancers5[i], "_select"), total)
  
}

rm(clin)
rm(mut)
rm(total)
rm(total_sub)
