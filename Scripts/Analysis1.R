##-----------------------------------------------------------------------------------------------------------
## This file will compare smoker vs nonsmoker data from 12 cancer datasets
## In this script, non smokers are defined as all participants who have NA for cigarettes per day and
## Smokers are all patients who have a value for cigrattes per day
##-----------------------------------------------------------------------------------------------------------

##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

source("../Scripts/functions.R")
setwd("~/Research/BiolinksAnalysis/Datasets")
##-----------------------------------------------------------------------------------------------------------
## Read in the data
##-----------------------------------------------------------------------------------------------------------
Data_Names <- c("BLCA", "HNSC")#, "ESCA", "KICH", "KIRC", "LUSC", "PAAD", "STAD", "LIHC", "KIRP", "LUAD")

## LAML was removed for having no smoking data

pipelines <- unique(HNSC$pipeline)
## This will read in each Cancer and determine the important sites
for (i in 1:length(Data_Names)){
  df <- read_delim(paste0(Data_Names[i], "_select.csv"), delim = "\t")
  df <- df %>% filter(mutation_status == "Somatic",
                             variant_classification %in% c("Missense_Mutation", 
                                                           "Nonsense_Mutations", "Silent",
                                                           "Frame_Shift_Del",
                                                           "Frame_Shift_Ins", "In_Frame_Del", 
                                                           "In_Frame_Ins", "Indel")) 
  valid_sites <- smoke_sites(df)
  for(j in 1:length(pipelines)){
    for(k in 1:length(valid_sites)){
    ## Subset the data into smokers and non-smokers. These will be in a loop that separates cancer, pipeline, and location
    ## Set_1 are non-smokers, Set_2 are smokers
      Set_1 <- df %>% filter(is.na(cigarettes_per_day), pipeline == pipelines[j], tissue_or_organ_of_origin == valid_sites[k])
      Set_2 <- df %>% filter(!(is.na(cigarettes_per_day)), pipeline == pipelines[j], tissue_or_organ_of_origin == valid_sites[k])
    
    ## For each dataset, find the total number of somatic point mutations per person 
      Set_2_Sum <- number_mut(Set_2, "smoker")
      Set_1_Sum <- number_mut(Set_1, "non-smoker")
      Sum <- rbind(Set_2_Sum, Set_1_Sum)
      total_mut <- ggplot(Sum) + geom_boxplot(aes(x = status, y = total))

    ## Calculate the pValue and quartile data for both groups
      pVal_Sum <- wilcox.test(total ~ status, data = Sum, paired = F)$p.value
      Overview_Set_2 <- summary_data(Set_2_Sum)
      Overview_Set_1 <- summary_data(Set_1_Sum)

    ## Combine summary data
      Overview_Set_2$status <- "smoker"
      Overview_Set_1$status <- "non-smoker"
      Overview <- rbind(Overview_Set_2, Overview_Set_1)
      Overview$pValue <- pVal_Sum
    
    ## Save summary data and plots
      write_delim(Overview, file.path("../Output/Smoke/total_mut_pval/", paste0(pipelines[j], "/", 
                                                                                Data_Names[i],"_", valid_sites[k], 
                                                                                ".csv")), delim  = "\t")
      ggsave(total_mut, file = file.path("../Output/Smoke/total_mut_graphs/", paste0(pipelines[j], "/", 
                                                                                     Data_Names[i], "_", valid_sites[k], 
                                                                                     ".jpg")), width = 6,
             height = 6, units = "in")
    }
  }
}



















