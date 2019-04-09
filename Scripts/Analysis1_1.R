## Author: Chris Holt
## Purpose: Compares nucleotide change frequencies, total amount of somatic point mutations, 
## and Transition vs Transversion Frequencies between non-smokers and smokers
## Date Created: Feb/2019
## Date of Last Update: 8/Apr/2019

##-----------------------------------------------------------------------------------------------------------
## This file will compare smoker vs nonsmoker data from 10 cancer datasets
## In this script, non smokers are defined as all participants who have NA for cigarettes per day and
## Smokers are all patients who have a value for cigrattes per day. This will perform the same analysis
## as Analysis1.R except that it will exclude valid_sites and will only look at the overall TCGA Cancer
##-----------------------------------------------------------------------------------------------------------

##-----------------------------------------------------------------------------------------------------------
## Clear Global Env
##-----------------------------------------------------------------------------------------------------------

rm(list = ls())

##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

source("~/Research/BiolinksAnalysis/Scripts/functions.R")
##-----------------------------------------------------------------------------------------------------------
## Read in the data
##-----------------------------------------------------------------------------------------------------------

## LAML was removed for having no smoking data, ESCA has no valid sites for smoking data
## Kirc, STAD, LIHC has no smoking data
##-----------------------------------------------------------------------------------------------------------
## The cancers were split into sets of 2 for processing reason. To download the data, select the 
## Data_Names vector and run the script below.
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
          ## Set_1 are smokers/male/EurAmr, Set_2 are non-smokers/female/AfrAmr
          Set_1 <- df %>% filter(status == cats[1], pipeline == pipelines[j])
          Set_2 <- df %>% filter(status == cats[2] , pipeline == pipelines[j])
          
          
          ##-----------------------------------------------------------------------------------------------------------
          ## For each Cancer, determine the valid sites, then for each location and each pipeline,
          ## calculate the pValue for the distribution of all somatic point mutations
          ##-----------------------------------------------------------------------------------------------------------
          
          ## For each dataset, find the total number of somatic point mutations per person 
          Set_2_Sum <- number_mut(Set_2, cats[2])
          Set_1_Sum <- number_mut(Set_1, cats[1])
          Sum <- rbind(Set_2_Sum, Set_1_Sum)
          total_mut <- ggplot(Sum) + geom_boxplot(aes(x = status, y = total))
          
          ## Calculate the pValue and quartile data for both groups
          tryCatch({
            pVal_Sum <- wilcox.test(total ~ status, data = Sum, paired = F)$p.value
            Overview_Set_2 <- summary_data(Set_2_Sum)
            Overview_Set_1 <- summary_data(Set_1_Sum)
            
            ## Combine summary data
            Overview_Set_2$status <- cats[2]
            Overview_Set_1$status <- cats[1]
            Overview <- rbind(Overview_Set_2, Overview_Set_1)
            Overview$pValue <- pVal_Sum
            
            ## Save summary data and plots
            write_delim(Overview, file.path(paste0("~/Research/BiolinksAnalysis/Output/",category, "/pValues/total_mut_pval/", pipelines[j], "/", "Overall", 
                                                                                                                       "_", Data_Names[i],
                                                                                                                       ".tsv")), delim  = "\t")
          }, error = function(e){})
          ggsave(total_mut, file = file.path(paste0("~/Research/BiolinksAnalysis/Output/",category, "/pValues/total_mut_pval/", pipelines[j], "/", "Overall", 
                                                                                                                        "_",
                                                                                                                        Data_Names[i], 
                                                                                                                        ".jpg")), width = 6,
                 height = 6, units = "in")
          
          ##-----------------------------------------------------------------------------------------------------------
          ## For each Cancer, determine the valid sites, then for each location and each pipeline,
          ## calculate the pValue between the percentage of nucleotide changes
          ##-----------------------------------------------------------------------------------------------------------
          Set_1_nuc <- nucChange_Sum(Set_1, cats[1], abbr[1])
          Set_2_nuc <- nucChange_Sum(Set_2, cats[2], abbr[2])
         
          combined_nuc <- rbind(Set_1_nuc, Set_2_nuc)
          combined_nuc_1 <- combined_nuc %>%  group_by(status, Var1, Var2, abbr) %>% 
            summarise(total = sum(Freq)) %>% group_by(Var1) %>% mutate(perc = (total/sum(total))*100)
          
          
          combined_nuc_1 <- combined_nuc_1 %>% group_by(Var1) %>% mutate(perc = (total/sum(total))*100)
          nucChanges <- as.character(unique(combined_nuc_1$Var2))
          
          nuc_graph <- ggplot(combined_nuc_1) + geom_boxplot(aes(x = Var2, y = perc, fill = status)) + coord_flip()
          
          ggsave(nuc_graph, file = file.path(paste0("~/Research/BiolinksAnalysis/Output/",category,"/Graphs/nucChange_Graph/", pipelines[j], "/","Overall", 
                                                                                                                        "_", 
                                                                                                                        Data_Names[i],
                                                                                                                        ".jpg")), width = 6,
                 height = 6, units = "in")
          for(n in 1:length(nucChanges)){
            change <- combined_nuc_1 %>% filter(Var2 == nucChanges[n])
            tryCatch({
              pValue <- wilcox.test(total~status, data = change, paired = F)$p.value
              
              
              change_1_summary <- change %>% filter(abbr == abbr[1])
              change_1_summary_data <- summary_data(change_1_summary)
              change_1_summary_data$status <- abbr[1]
              
              change_2_summary <- change %>% filter(abbr == abbr[2])
              change_2_summary_data <- summary_data(change_2_summary)
              change_2_summary_data$status <- abbr[2]
              
              change_summary <- rbind(change_1_summary_data, change_2_summary_data)
              change_summary$pValue <- pValue
              write_delim(change_summary, file.path(paste0("~/Research/BiolinksAnalysis/Output/",category,"/pValues/nucChange_pVal/",pipelines[j], "/", "Overall", 
                                                                                                                               "_",
                                                                                                                               Data_Names[i],
                                                                                                                               "_", nucChanges[n],
                                                                                                                               ".tsv")), delim  = "\t")
              
            }, error = function(e){})
            
          }
          
          ##-----------------------------------------------------------------------------------------------------------
          ## For each Cancer, determine the valid sites, then for each location and each pipeline,
          ## calculate the pValue  for Ti vs Tv
          ##-----------------------------------------------------------------------------------------------------------
          
          Set_1_TiTv <- TiTv_count(Set_1, cats[1], abbr[1])
          Set_2_TiTv <- TiTv_count(Set_2, cats[2], abbr[2])
          
          combined_Ti_Tv <- rbind(Set_1_TiTv, Set_2_TiTv)
          
          TiTv <- ggplot(combined_Ti_Tv) + geom_boxplot(aes(x = status, y = pc, fill = Var2))
          ggsave(TiTv, file = 
                   file.path(paste0("~/Research/BiolinksAnalysis/Output/", category, "/Graphs/TiTv_Graph/",
                                    pipelines[j], "/", "Overall", "_",Data_Names[i], ".jpg")), width = 6, height = 6, units = "in")
          
          
          combined_Ti_Tv <- combined_Ti_Tv %>% rename("total" = "pc")
          type <- unique(combined_Ti_Tv$Var2)
          
          for(l in 1:length(type)){
            type_data <- combined_Ti_Tv %>% filter(Var2 == type[l])
            tryCatch({
              type_pval <- wilcox.test(total ~ status, data = type_data, paired = F)$p.value
              
              type_data_1 <- type_data %>% filter(abbr == abbr[1])
              type_data_1_summary <- summary_data(type_data_1)
              type_data_1_summary$status <- abbr[1]
              
              type_data_2 <- type_data %>% filter(abbr == abbr[2])
              type_data_2_summary <- summary_data(type_data_2)
              type_data_2_summary$status <- abbr[2]
              
              type_data_combined <- rbind(type_data_1_summary,type_data_2_summary)
              
              type_data_combined$pvalue <- type_pval
              
              write_delim(change_summary, file.path(paste0("~/Research/BiolinksAnalysis/Output/",category, "/pValues/TiTv_pVal/", pipelines[j], "/", "Overall", 
                                                                                                                          "_",
                                                                                                                          Data_Names[i],
                                                                                                                          "_", type[l],
                                                                                                                          ".tsv")), delim  = "\t")
            }, error = function(e){})
          }
    }
  }
}

Data_Names <- c("BLCA", "HNSC", "KICH", "LUAD", "LUSC", "PAAD", "KIRP")

main(Data_Names, "Smoke")
main(Data_Names, "Gender")
main(Data_Names, "Race")












