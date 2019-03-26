## Author: Chris Holt
## Purpose: Compares gene frequencies between males and females
## Date Created: 13/Mar/2019
## Date of Last Update: 26/Mar/2019
##-------------------------------------------------------------------------------------
## This script will find the frequencies of genes w/ somatic, point mutations 
## between males and females and calculate pValues between males and females
##-------------------------------------------------------------------------------------

##-----------------------------------------------------------------------------------------------------------
## Clear Global Env
##-----------------------------------------------------------------------------------------------------------

rm(list = ls())

##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------
source("~/Research/BiolinksAnalysis/Scripts/functions.R")
source("~/Research/BiolinksAnalysis/Scripts/functions3.R")
##-----------------------------------------------------------------------------------------------------------
## Read in the data
##-----------------------------------------------------------------------------------------------------------
Data_Names <- c("BLCA", "HNSC")
Data_Names <- c("KICH", "LUAD")
Data_Names <- c("KIRC", "LUSC")
Data_Names <- c("PAAD", "STAD")
Data_Names <- c("LIHC", "KIRP")

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
    Set_1 <- df %>% filter(gender == "male", pipeline == pipelines[j])
    Set_2 <- df %>% filter(gender == "female", pipeline == pipelines[j])
    
    Set_1_people <- length(unique(Set_1$tumor_barcode))
    Set_2_people <- length(unique(Set_2$tumor_barcode))
    
    genes_1 <- as.data.frame(gene_count(Set_1))
    genes_2 <- as.data.frame(gene_count(Set_2))
    
    same_genes <- as.data.frame(intersect(genes_1$hugo_symbol, genes_2$hugo_symbol)) %>% 
      rename("hugo_symbol" = `intersect(genes_1$hugo_symbol, genes_2$hugo_symbol)`) %>% mutate_if(is.factor, as.character)
    
    diff_1 <- diff_genes(genes_1, genes_2) ## Genes in genes_1 that are not in genes_2
    diff_2 <-diff_genes(genes_2, genes_1) ## Genes in genes_2 that are not in genes_1
    
    diff_1$n <- as.integer(0)
    diff_2$n <- as.integer(0)
    
    num_people_1 <- genes_1 %>% select(hugo_symbol, tumor_barcode) %>% group_by(hugo_symbol) %>% count()
    num_people_2 <- genes_2 %>% select(hugo_symbol, tumor_barcode) %>% group_by(hugo_symbol) %>% count()
    
    genes_1 <- genes_1 %>% select(hugo_symbol, n)
    genes_2 <- genes_2 %>% select(hugo_symbol, n)
    
    genes_1 <- rbind(genes_1, diff_2)
    genes_2 <- rbind(genes_2, diff_1)
    
    num_same_genes <- nrow(same_genes)
    num_1_genes <- genes_1 %>% filter(n > 0) %>% nrow()
    num_2_genes <- genes_2 %>% filter(n > 0) %>% nrow()
    
    
    genes_1$status <- "M"
    genes_2$status <- "F"
    
    genes_1 <- genes_1 %>% left_join(num_people_1, by = "hugo_symbol")
    genes_2 <- genes_2 %>% left_join(num_people_2, by = "hugo_symbol")
    genes_1$n.y[is.na(genes_1$n.y)] <- 0
    genes_2$n.y[is.na(genes_2$n.y)] <- 0
    
    
    genes_1 <- genes_1 %>% group_by(hugo_symbol) %>% mutate(perc = (n.x/n.y)*100)
    genes_2 <- genes_2 %>% group_by(hugo_symbol) %>% mutate(perc = (n.x/n.y)*100)
    
    genes <- rbind(genes_1, genes_2)
    genes$perc[is.nan(genes$perc)] <- 0
    
    setwd("~/Research/BiolinksAnalysis/Output/Gender/Genes_Pvalues")
    write_delim(genes, paste0(pipelines[j], "_", Data_Names[i], "_Gene_pvalues", ".tsv"), delim = "\t")
    setwd("~/Research/BiolinksAnalysis/")
    
    # 
    # uniq_genes <- as.character(unique(genes$hugo_symbol))
    # 
    # genes_pval <- data.frame()
    # data <- tibble("hugo_symbol", "pval")
    # for(k in 1:length(uniq_genes)){
    #   gene_data <- genes %>% filter(hugo_symbol == uniq_genes[k])
    #   pVal <- wilcox.test(perc ~ status, data = gene_data, paired = F)$p.value ##Keep on getting pValues of 1
    #   
    #   
    #   data$hugo_symbol <- uniq_genes[k]
    #   data$pVal <- pVal
    #   
    #   genes_pval <- rbind(genes_pval, data)
    #   genes_pval$pVal[is.nan(genes_pval$pVal)] <- 1.1
    # }
    # genes_pavl <- genes_pval %>% select(hugo_symbol, pVal)
  }
}



