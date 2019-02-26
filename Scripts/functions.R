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
## This will utilise the functions above and download and merge the entire mutational and clinical datasets for 
## each cancer
##-----------------------------------------------------------------------------------------------------------

merged <- function(Cancers){
  pipelines <- c("somaticsniper", "muse", "mutect", "varscan2")
  for (i in 1:length(Cancers)){
    combined <- data.frame()
    clinical <- download_clinical(Cancers[i])
    for(j in 1:length(pipelines)){
      mutational <- download_mutational(Cancers[i], pipelines[j])
      mutational <- classify_Changes(mutational)
      combined <- rbind(combined, mutational)
    }
    total <- merge(clinical, combined, by = "tumor_barcode")
  return(total)
  }
}

merge_ind <- function(clin, mut){
  total <- merge(clin, mut, by = "tumor_barcode")
  return(total)
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
    message("category does not equal gender")
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
## This function will take a dataset and return a data table that consists of the unique ID and the
## raw number of somatic point mutations for each person
## type refers to either "smoker" or "non-smoker"
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
  File_Sum$status <- type
  
  return(File_Sum)
}
##-----------------------------------------------------------------------------------------------------------
## This function will take a dataset and create a dataframe of the summary data to find quartile info
## and add in pValue
##-----------------------------------------------------------------------------------------------------------
summary_data <- function(Data, type){
  
}

##-----------------------------------------------------------------------------------------------------------
## End of Script
##-----------------------------------------------------------------------------------------------------------
















