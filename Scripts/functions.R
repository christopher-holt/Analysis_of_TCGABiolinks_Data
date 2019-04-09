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
  pot_sites <- unique(Data$tissue_or_organ_of_origin)
  new_site <- character()
  for(i in 1:length(pot_sites)){
    df <- Data  %>% filter(mutation_status == "Somatic",
                           variant_classification %in% c("Missense_Mutation", 
                                                         "Nonsense_Mutations", "Silent",
                                                         "Frame_Shift_Del",
                                                         "Frame_Shift_Ins", "In_Frame_Del", 
                                                         "In_Frame_Ins", "Indel"), 
                           tissue_or_organ_of_origin == pot_sites[i])
    if(category == "Smoke"){
      df_1 <- df %>% filter(is.na(cigarettes_per_day))
      df_2 <- df %>% filter(!(is.na(cigarettes_per_day)))
      
      if((nrow(df_1) > 0 & nrow(df_2) > 0) & ((length(unique(df_1$tumor_barcode)) >= 3 & length(unique(df_2$tumor_barcode)) >= 3))){
        new_site <- c(new_site, pot_sites[i])
      }
    } else if(category == "Race"){
      df_1 <- df %>% filter(race == "black or african american")
      df_2 <- df %>% filter(race == "white")
      
      if((nrow(df_1) > 0 & nrow(df_2) > 0) & ((length(unique(df_1$tumor_barcode)) >= 3 & length(unique(df_2$tumor_barcode)) >= 3))){
        new_site <- c(new_site, pot_sites[i])
      }
    } else{
      df_1 <- df %>% filter(gender == "male")
      df_2 <- df %>% filter(gender == "female")
      
      if((nrow(df_1) > 0 & nrow(df_2) > 0) & ((length(unique(df_1$tumor_barcode)) >= 3 & length(unique(df_2$tumor_barcode)) >= 3))){
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
## This function will take a dataset and, for each unique person, count the number of times a specific genes
## appears
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
## This function will take a dataset, and return the number of specific nucleotide changes per person
##-----------------------------------------------------------------------------------------------------------

nucChange_Sum <- function(df, status, Initial){
  df_1 <- df %>% filter(!(nucChange %in% c("deletion", "insertion", "other")) )  %>% mutate_if(is.factor,as.character)
  df_1_table <- as.data.frame(table(df_1$tumor_barcode,df_1$nucChange)) %>% mutate_if(is.factor,as.character)
  df_1_table$status <- paste0(length(unique(df_1_table$Var1)), status)
  df_1_table$abbr <- paste0(Initial)
  return(df_1_table)
}

##-----------------------------------------------------------------------------------------------------------
## This function will take a dataset, count the number of TiTv per person
##-----------------------------------------------------------------------------------------------------------
TiTv_count <- function(df, status, Initial){
  df_1 <- df %>% filter(Mut_Type %in% c("Ti", "Tv")) %>% mutate_if(is.factor,as.character)
  df_1_table <- as.data.frame(table(df_1$tumor_barcode,df_1$Mut_Type)) %>% mutate_if(is.factor,as.character)
  
  df_sum <- df_1_table %>% group_by(Var1) %>% mutate(pc = Freq/sum(Freq)*100)
  
  df_sum$status <- paste0(length(unique(df_1_table$Var1)), "-", status)
  df_sum$abbr <- paste0(Initial)
  
  return(df_sum)
}
##-----------------------------------------------------------------------------------------------------------
## This function will take a dataset and return a data table that consists of the unique ID and the
## raw number of somatic point mutations for each person
## type refers to attribute such as male vs female, AfrAmr vs EurAmr or Smoker vs Non-smoker.
## Files must be subsetted first. So file must only contain values of one interest group eg. 
## either only male or only female
##-----------------------------------------------------------------------------------------------------------

number_mut <- function(File, type){
  File_1 <- File %>% filter( mutation_status == "Somatic",
                             variant_classification %in% c("Missense_Mutation", 
                                                           "Nonsense_Mutations", "Silent",
                                                           "Frame_Shift_Del",
                                                           "Frame_Shift_Ins", "In_Frame_Del", 
                                                           "In_Frame_Ins", "Indel"))
  File_1_table <- as.data.frame(table(File_1$tumor_barcode, File_1$variant_classification)) %>%mutate_if(is.factor,as.character)
  File_Sum <- File_1_table %>% group_by(Var1) %>% summarise(total = sum(Freq))
  File_Sum$status <- paste0(length(unique(File_1_table$Var1)), "-",type)
  
  
  return(File_Sum)
}

##-----------------------------------------------------------------------------------------------------------
## End of Script
##-----------------------------------------------------------------------------------------------------------
















