##-------------------------------------------------------------------------------------
## This script will find the frequencies of genes w/ somatic, point mutations 
## between smokers and non smokers
##-------------------------------------------------------------------------------------

##-----------------------------------------------------------------------------------------------------------
## Clear Global Env
##-----------------------------------------------------------------------------------------------------------

rm(list = ls())

##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

source("Scripts/functions.R")
source("Scripts/functions1.R")

##----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------
## Read in the data
##-----------------------------------------------------------------------------------------------------------
Data_Names <- c("BLCA", "HNSC")
Data_Names <- c("KICH", "LUAD")
Data_Names <- c("LUSC")
Data_Names <- c("PAAD")
Data_Names <- c("KIRP")

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
  
  for(j in 1:length(pipelines)){
    ## Subset the data into smokers and non-smokers. These will be in a loop that separates cancer, pipeline, and location
    ## Set_1 are non-smokers, Set_2 are smokers
    Set_1 <- df %>% filter(is.na(cigarettes_per_day), pipeline == pipelines[j])
    Set_2 <- df %>% filter(!(is.na(cigarettes_per_day)), pipeline == pipelines[j])
    
    
    genes_1 <- gene_count(Set_1)
    genes_2 <- gene_count(Set_2)
    
    
    
    same_genes <- as.data.frame(intersect(genes_1$hugo_symbol, genes_2$hugo_symbol)) %>% 
                rename("genes" = `intersect(genes_1$hugo_symbol, genes_2$hugo_symbol)`) %>% mutate_if(is.factor, as.character)
    
    diff_1 <- diff_genes(genes_1, genes_2) ## Genes in genes_1 that are not in genes_2
    diff_2 <-diff_genes(genes_2, genes_1) ## Genes in genes_2 that are not in genes_1
    
    
    

  }
}



