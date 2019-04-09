## Author: Chris Holt
## Purpose: Compares gene frequencies between smokers and non smokers
## Date Created: 13/Mar/2019
## Date of Last Update: 9/Apr/2019

##-------------------------------------------------------------------------------------
## This script will find the frequencies of genes w/ somatic, point mutations 
## between smokers and non smokers and calculate pValues
##-------------------------------------------------------------------------------------

##-----------------------------------------------------------------------------------------------------------
## Clear Global Env
##-----------------------------------------------------------------------------------------------------------

rm(list = ls())

##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

source("Scripts/functions.R")

##----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------
## Read in the data
##-----------------------------------------------------------------------------------------------------------
main <- function(Data_Names, category){
  pipelines <- c("muse", "mutect", "somaticsniper", "varscan2")
  for (i in 1:length(Data_Names)){
    setwd("~/Research/BiolinksAnalysis/Datasets")
    df <- read_delim(paste0(Data_Names[i], "_select.csv"), delim = "\t")
    setwd("~/Research/BiolinksAnalysis/")
    df <- df %>% filter(mutation_status == "Somatic",
                        variant_classification %in% c("Missense_Mutation", 
                                                      "Nonsense_Mutations", "Silent",
                                                      "Frame_Shift_Del",
                                                      "Frame_Shift_Ins", "In_Frame_Del", 
                                                      "In_Frame_Ins", "Indel"))
    df <- df %>% mutate(smoke_status = ifelse(is.na(cigarettes_per_day), "non-smoker", "smoker"))
    colnames(df) <- str_to_lower(colnames(df))
    
    if(category == "Smoke"){
      df <- df %>% rename("status" = "smoke_status")
      df <- df %>% mutate(abbr = ifelse(status == "non-smoker", "NS", "S"))
    } else if(category == "Race"){
      df <- df %>% filter(race == "white" | race == "black or african american") %>% rename("status" = "race")
      df <- df %>% mutate(abbr = ifelse(status == "white", "EurAmr", "AfrAmr"))
    } else if(category == "Gender"){
      df <- df %>% rename("status" = "gender")
      df <- df %>% mutate(abbr = ifelse(status == "male", "M", "F"))
    }
    
    cats <- unique(df$status)
    abbr <- unique(df$abbr)
  for(j in 1:length(pipelines)){
    ## Subset the data into smokers and non-smokers. These will be in a loop that separates cancer, pipeline, and location
    ## ## Set_1 are smokers/male/EurAmr, Set_2 are non-smokers/female/AfrAmr
    Set_1 <- df %>% filter(status == cats[1], pipeline == pipelines[j])
    Set_2 <- df %>% filter(status == cats[2] , pipeline == pipelines[j])
    
    Set_1_people <- length(unique(Set_1$tumor_barcode))
    Set_2_people <- length(unique(Set_2$tumor_barcode))
    
    ## for each person, count the number of times a gene is mutated per person
    genes_1 <- as.data.frame(gene_count(Set_1)) 
    genes_2 <- as.data.frame(gene_count(Set_2))
    
    same_genes <- as.data.frame(intersect(genes_1$hugo_symbol, genes_2$hugo_symbol)) %>% 
                rename("hugo_symbol" = `intersect(genes_1$hugo_symbol, genes_2$hugo_symbol)`) %>% mutate_if(is.factor, as.character)
    
    diff_1 <- diff_genes(genes_1, genes_2) ## Genes in genes_1 that are not in genes_2
    diff_2 <-diff_genes(genes_2, genes_1) ## Genes in genes_2 that are not in genes_1
    
    diff_1$n <- as.integer(0)
    diff_2$n <- as.integer(0)
    
    ## number of people who have a mutation in that gene
    num_people_1 <- genes_1 %>% select(hugo_symbol, tumor_barcode) %>% group_by(hugo_symbol) %>% count() 
    num_people_2 <- genes_2 %>% select(hugo_symbol, tumor_barcode) %>% group_by(hugo_symbol) %>% count()
    
    genes_1 <- genes_1 %>% select(hugo_symbol, n)
    genes_2 <- genes_2 %>% select(hugo_symbol, n)
    
    genes_1 <- rbind(genes_1, diff_2)
    genes_2 <- rbind(genes_2, diff_1)
    
    
    genes_1$status <- abbr[1]
    genes_2$status <- abbr[2]
    
    genes_1 <- genes_1 %>% left_join(num_people_1, by = "hugo_symbol")
    genes_2 <- genes_2 %>% left_join(num_people_2, by = "hugo_symbol")
    genes_1$n.y[is.na(genes_1$n.y)] <- 0
    genes_2$n.y[is.na(genes_2$n.y)] <- 0
    
    ## n.x is the number of mutations a specific person has for gene
    ## n.y is the number of people in the whole population who have >=1 mutation in that gene
    
    genes_1 <- genes_1 %>% group_by(hugo_symbol) %>% mutate(perc = (n.x/n.y)*100)
    genes_2 <- genes_2 %>% group_by(hugo_symbol) %>% mutate(perc = (n.x/n.y)*100)
    
    genes <- rbind(genes_1, genes_2)
    genes$perc[is.nan(genes$perc)] <- 0
    
    setwd(paste0("~/Research/BiolinksAnalysis/Output/", category, "/Genes_Pvalues"))
    write_delim(genes, paste0(pipelines[j], "_", Data_Names[i], "_Gene_pvalues", ".tsv"), delim = "\t")
    setwd("~/Research/BiolinksAnalysis/")
  }
}
}

Data_Names <- c("BLCA", "HNSC", "KICH", "LUAD", "LUSC", "PAAD", "KIRP")

## Gender/Race analysis only
extra_data_names <- c("KIRC", "PAAD", "LIHC")

main(Data_Names, "Smoke")
main(Data_Names, "Gender")
main(Data_Names, "Race")

main(extra_data_names, "Gender")
main(extra_data_names, "Race")


##-----------------------------------------------------------------------------------------------------------
## End of Script
##-----------------------------------------------------------------------------------------------------------
















