## Author: Chris Holt
## Purpose: Looks at mutational burden and age between Smokers and Non-smokers
## Date Created: 25/Mar/2019
## Date of Last Update: 9/Apr/2019

##-----------------------------------------------------------------------------------------------------------
## This file will compare smoker vs nonsmoker data from 10 cancer datasets
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
source("~/Research/BiolinksAnalysis/Scripts/functions1.R")

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
    Set_1 <- df %>% filter(status == cats[1], pipeline == pipelines[j]) %>% select(age_at_diagnosis, hugo_symbol)
    Set_2 <- df %>% filter(status == cats[2], pipeline == pipelines[j]) %>% select(age_at_diagnosis, hugo_symbol)
    
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
    
    age_1$status <- abbr[1]
    age_2$status <- abbr[2]
    
    age <- rbind(age_1, age_2)
    
    age_data <- age %>% split(.$status) %>% map(~lm(n~new_age, data = .))
    
    r_squared <- age_data %>% map(summary) %>% map_dbl("r.squared")
    
    dat <- tribble(
      ~NS, ~S, ~info,
      r_squared[1], r_squared[2], "rsquared values comparing age vs mutational burden"
    )
    
    age_graph <- ggplot(age) + geom_point(aes(new_age, n))+
      geom_smooth(method = 'lm', aes(new_age, n), se = F) +
      facet_wrap(~status)
    
    write_delim(dat, file.path(paste0("Output/", category, "/Age/files/", pipelines[j], "/", Data_Names[i], "_rsquared.tsv")),
                delim = '\t')
    
    ggsave(age_graph, file = file.path(paste0("Output/", category, "/Age/Graphs/", pipelines[j], "/", Data_Names[i], "_age_burden.jpg")),
           width = 6, height = 6, units = "in")
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




    