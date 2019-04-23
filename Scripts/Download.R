## Author: Chris Holt
## Purpose: To download and standardise TCGAbiolinks data
## Date Created: Feb/2019
## Date of Last Update: 23/Apr/2019

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

total_Cancers <- getGDCprojects()$project_id %>%
  enframe() %>% select(value) %>% rename("cancer" = "value") %>% 
  filter(grepl("^TCGA.*$", cancer))

total_Cancers$cancer <- total_Cancers$cancer %>% str_remove("TCGA-") %>% sort()
total_Cancers <- total_Cancers %>% pull()

## A select number of cancers

Cancers <- c("LAML", "BLCA", "ESCA","KICH", "KIRC","KIRP", "LUAD","LUSC", "PAAD","STAD", "LIHC")
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
}


rm(clin)
rm(mut)
rm(total)
rm(total_sub)
