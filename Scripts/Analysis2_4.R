##-----------------------------------------------------------------------------------------------------------
## This file will compare african_american and EurAmr data from 10 cancer datasets
## In this script, non smokers are defined as all participants who have NA for cigarettes per day and
## Smokers are all patients who have a value for cigrattes per day. This script will look at 
## the mutational burden compared to age
##-----------------------------------------------------------------------------------------------------------

##-----------------------------------------------------------------------------------------------------------
## Clear Global Env
##-----------------------------------------------------------------------------------------------------------

rm(list = ls())

##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

source("~/Research/BiolinksAnalysis/Scripts/functions.R")
source("~/Research/BiolinksAnalysis/Scripts/functions2.R")
##-----------------------------------------------------------------------------------------------------------
## Read in the data
##-----------------------------------------------------------------------------------------------------------
Data_Names <- c("BLCA", "HNSC")
Data_Names <- c("KICH", "LUAD")
Data_Names <- c("KIRC", "LUSC")
Data_Names <- c("PAAD", "STAD")
Data_Names <- c("LIHC", "KIRP")

##-----------------------------------------------------------------------------------------------------------
## The cancers were split into sets of 2 for processing reason. To download the data, select the 
## Data_Names vector and run the script below.
##-----------------------------------------------------------------------------------------------------------

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
    Set_1 <- df %>% filter(race == "black or african american", pipeline == pipelines[j]) %>% select(age_at_diagnosis, hugo_symbol)
    Set_2 <- df %>% filter(race == "white", pipeline == pipelines[j]) %>% select(age_at_diagnosis, hugo_symbol)
    
    age_1 <- Set_1 %>% 
      group_by(age_at_diagnosis) %>% 
      count() %>% 
      arrange(age_at_diagnosis) %>% 
      mutate(new_age = round(age_at_diagnosis/365))
    age_2 <- Set_2 %>% 
      group_by(age_at_diagnosis) %>%
      count() %>%
      arrange(age_at_diagnosis) %>% 
      mutate(new_age = round(age_at_diagnosis/365))
    
    age_1$status <- "AfrAmr"
    age_2$status <- "EurAmr"
    
    age <- rbind(age_1, age_2)
    
    age_data <- age %>% split(.$status) %>% map(~lm(new_age~n, data = .))
    
    r_squared <- age_data %>% map(summary) %>% map_dbl("r.squared")
    
    dat <- tribble(
      ~NS, ~S, ~info,
      r_squared[1], r_squared[2], "rsquared values comparing age vs mutational burden"
    )
    
    age_graph <- ggplot(age) + geom_point(aes(new_age, n)) + geom_smooth(method = 'lm', aes(new_age, n), se = F) +
      facet_wrap(~status)
    
    write_delim(dat, file.path("Output/Race/Age/files/", paste0(pipelines[j], "/", Data_Names[i], "_rsquared.tsv")),
                delim = '\t')
    
    ggsave(age_graph, file = file.path("Output/Race/Age/Graphs/", paste0(pipelines[j], "/", Data_Names[i], "_age_burden.jpg")),
           width = 6, height = 6, units = "in")
  }
}
