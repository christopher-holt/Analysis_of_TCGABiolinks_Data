##-----------------------------------------------------------------------------------------------------------
## This file will generate all functions used in specifically Analysis1.R
##-----------------------------------------------------------------------------------------------------------


##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

library(tidyverse)
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
## This function will take a dataset and determine which locations have enough data for smokers and 
## non-smokers
##-----------------------------------------------------------------------------------------------------------

smoke_sites <- function(Data){
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
    
    df_1 <- df %>% filter(is.na(cigarettes_per_day))
    df_2 <- df %>% filter(!(is.na(cigarettes_per_day)))
    
    if((nrow(df_1) > 0 & nrow(df_2) > 0) & ((length(unique(df_1$tumor_barcode)) >= 3 & length(unique(df_2$tumor_barcode)) >= 3))){
      new_site <- c(new_site, pot_sites[i])
    }
  }
  return(new_site)
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
## End of Script
##-----------------------------------------------------------------------------------------------------------


