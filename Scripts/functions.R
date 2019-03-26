## Author: Chris Holt
## Purpose: Generate all functions used in analysis for this project
## Date Created: Feb/2019
## Date of Last Update: 26/Mar/2019

##-----------------------------------------------------------------------------------------------------------
## This file will generate all general functions need throughout this project 
##-----------------------------------------------------------------------------------------------------------


##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

library(tidyverse)
library(TCGAbiolinks)

##-----------------------------------------------------------------------------------------------------------
## Downloads the clinical data
##-----------------------------------------------------------------------------------------------------------

download_clinical <- function(Data){
  clinical <- GDCquery_clinic(paste0("TCGA-",Data))
  clinical$tumor_barcode <- str_replace_all(clinical$submitter_id, pattern = '-',replacement = '.')
  clinical$tumor_barcode <- str_sub(clinical$tumor_barcode, start = 1, end = 12)
  clinical$tumor_barcode <- sapply(clinical$tumor_barcode, tolower)
  clinical$submitter_id <- NULL
  colnames(clinical) <- tolower(colnames(clinical))
  
  return(clinical)
}

##-----------------------------------------------------------------------------------------------------------
## Dowloads and merges the mutational data from each of the four pipelines
## as well as formatting two columns to remove characters
##-----------------------------------------------------------------------------------------------------------

download_mutational<- function(file){
  pipelines <- c("somaticsniper", "muse", "mutect", "varscan2")
  combined <- data.frame()
  
  for (i in 1:length(pipelines)){
    Data <- GDCquery_Maf(paste0(file), pipelines = pipelines[i])
    Data$tumor_barcode <- str_replace_all(Data$Tumor_Sample_Barcode, pattern = '-',replacement = '.')
    Data$tumor_barcode <- str_sub(Data$tumor_barcode, start = 1, end = 12)
    Data$tumor_barcode <- sapply(Data$tumor_barcode, tolower)
    Data$pipeline = paste0(pipelines[i])
    colnames(Data) <- tolower(colnames(Data))


    hgvsc1 <- str_replace_all(Data$hgvsc, pattern = "^c.", replacement = "")
    c <- as.data.frame(hgvsc1)
    c <- type.convert(c, as.is = T)
    Data_1 <- cbind(Data, c) 
    
    combined <- rbind(combined, Data_1)
    
  }
  return(combined)
}



##-----------------------------------------------------------------------------------------------------------
## This will take two columns containing nucleotide changes and create their own column
## as well as classify the nucleotide changes ad transitions and transversion
##-----------------------------------------------------------------------------------------------------------

classify_Changes <- function(Data){
  ti_c = c("^.*G>A$", "^.*A>G$", "^.*C>T$", "^.*T>C$")
  tv_c = c("^.*C>A$", "^.*C>G$", "^.*A>C$", "^.*G>C$", "^.*G>T$", "^.*T>G$", "^.*T>A$","^.*A>T$")
  
  Data <- Data %>% mutate(Mut_Type = ifelse(grepl(paste(ti_c, collapse = '|'), Data$hgvsc1) , "Ti", 
                                                              ifelse(grepl(paste(tv_c, collapse = '|'), Data$hgvsc1), "Tv",
                                                                     ifelse(grepl("^.*ins.*$", Data$hgvsc1), "insertion",
                                                                            ifelse(grepl("^.*del.*$", Data$hgvsc1), "deletion", 
                                                                                   ifelse(grepl("^.$", Data$hgvsc1), "None",
                                                                                          ifelse(is.na(Data$hgvsc1), NA, "other")))))))
  
  
  Data <- Data %>% mutate(nucChange = ifelse(grepl("^.*G>A$",Data$hgvsc1) , "G > A", 
                                                               ifelse(grepl("^.*A>G$",Data$hgvsc1), "A > G", 
                                                                      ifelse(grepl("^.*C>T$",Data$hgvsc1), "C > T",
                                                                             ifelse(grepl("^.*T>C$",Data$hgvsc1), "T > C",
                                                                                    ifelse(grepl("^.*C>A$",Data$hgvsc1), "C > A" ,
                                                                                           ifelse(grepl("^.*C>G$",Data$hgvsc1), "C > G", 
                                                                                                  ifelse(grepl("^.*A>C$",Data$hgvsc1), "A > C",
                                                                                                         ifelse(grepl("^.*G>C$",Data$hgvsc1), "G > C",
                                                                                                                ifelse(grepl("^.*G>T$",Data$hgvsc1), "G > T",
                                                                                                                       ifelse(grepl("^.*T>G$",Data$hgvsc1), "T > G",
                                                                                                                              ifelse(grepl("^.*T>A$",Data$hgvsc1), "T > A",
                                                                                                                                     ifelse(grepl("^.*A>T$",Data$hgvsc1), "A > T",
                                                                                                                                            ifelse(grepl("deletion",Data$hgvsc1), "deletion",
                                                                                                                                                   ifelse(grepl("insertion", Data$hgvsc1), "insertion", "other" )))))))))))))))
  
  
  
  
}

##-----------------------------------------------------------------------------------------------------------
## Does each cancer site have enough data? This will only consider somatic point mutations. 
## This function only considers race and gender
##-----------------------------------------------------------------------------------------------------------

sites <- function(Data, category){
  pot_sites <- as.character(unique(Data$tissue_or_organ_of_origin))
  new_site <- character()
  category <- as.name(paste0(category))
  if ("white" %in% Data[[category]]){
    try <- Data%>% filter(race %in% c("white", "black or african american"))
    pot_sites <- as.character(unique(try$tissue_or_organ_of_origin))
    cats <- as.character(unique(try[[category]]))
    for (i in 1:length(pot_sites)){
      set1 <- try %>% filter(mutation_status == "Somatic",
                             variant_classification %in% c("Missense_Mutation", 
                                                           "Nonsense_Mutations", "Silent",
                                                           "Frame_Shift_Del",
                                                           "Frame_Shift_Ins", "In_Frame_Del", 
                                                           "In_Frame_Ins", "Indel"), 
                             tissue_or_organ_of_origin == pot_sites[i])
      
      set1_a <- set1 %>% filter(race == cats[1])
      set1_b <- set1 %>% filter(race == cats[2])
      
      if((nrow(set1_b) > 0 & nrow(set1_a) > 0) & ((length(unique(set1_b$tumor_barcode)) >= 3 & length(unique(set1_a$tumor_barcode)) >= 3))){
        new_site <- c(new_site, pot_sites[i])
      } 
    }
  } else {
    cats <- as.character(unique(Data[[category]]))
    for (i in 1:length(pot_sites)){
      set1 <- try %>% filter(mutation_status == "Somatic",
                             variant_classification %in% c("Missense_Mutation", 
                                                           "Nonsense_Mutations", "Silent",
                                                           "Frame_Shift_Del",
                                                           "Frame_Shift_Ins", "In_Frame_Del", 
                                                           "In_Frame_Ins", "Indel"), 
                             tissue_or_organ_of_origin == pot_sites[i])
      
      set1_a <- set1 %>% filter(gender == cats[1])
      set1_b <- set1 %>% filter(gender == cats[2])
      
      if((nrow(set1_b) > 0 & nrow(set1_a) > 0) & ((length(unique(set1_b$tumor_barcode)) >= 3 & length(unique(set1_a$tumor_barcode)) >= 3))){
        new_site <- c(new_site, pot_sites[i])
      } 
    }
  }
  return(new_site)
}
##-----------------------------------------------------------------------------------------------------------
## This function will take a dataset and create a dataframe of the summary data to find quartile info
## and add in pValue
##-----------------------------------------------------------------------------------------------------------
summary_data <- function(Data){
  summary_1 <- as.data.frame(do.call(cbind, lapply(Data, summary)))
  summary_1 <- add_rownames(summary_1, "Values")
  Overview <- summary_1 %>% select(Values, total)
  return(Overview)
  
}

##-----------------------------------------------------------------------------------------------------------
## This function will take a dataset and determine how many people have >= 1 mutation in that gene
##-----------------------------------------------------------------------------------------------------------

gene_count <- function(df){
  new_df <- df %>% select(tumor_barcode,hugo_symbol) %>% 
    group_by(tumor_barcode, hugo_symbol) %>% count()
  
  return(new_df)
}

##-----------------------------------------------------------------------------------------------------------
## This function will take two datasets and identify the different genes between the two
##-----------------------------------------------------------------------------------------------------------

diff_genes <- function(df, df2){
  df3 <- as.data.frame(setdiff(df$hugo_symbol, df2$hugo_symbol)) %>%
    rename("hugo_symbol" = `setdiff(df$hugo_symbol, df2$hugo_symbol)`) %>% mutate_if(is.factor, as.character)
  
  return(df3)
}


##-----------------------------------------------------------------------------------------------------------
## End of Script
##-----------------------------------------------------------------------------------------------------------
















