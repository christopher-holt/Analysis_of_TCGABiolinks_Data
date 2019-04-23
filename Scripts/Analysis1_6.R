## Author: Chris Holt
## Purpose: Looks at the variant_classification and determines signficance 
## Date Created: 23/Apr/2019
## Date of Last Update: 23/Apr/2019

##-----------------------------------------------------------------------------------------------------------
## Clear Global Env
##-----------------------------------------------------------------------------------------------------------

rm(list = ls())

##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

source("Scripts/functions.R")

##----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------
## Read in the data
##-----------------------------------------------------------------------------------------------------------
main <- function(Data_Names, category){
  pipelines <- c("muse", "mutect", "somaticsniper", "varscan2")
  for (i in 1:length(Data_Names)){
    setwd("~/Research/BiolinksAnalysis/Datasets")
    df <- read_delim(paste0(Data_Names[i], "_select.csv"), delim = "\t")
    setwd("~/Research/BiolinksAnalysis/")
    # df <- df %>% filter(mutation_status == "Somatic",
    #                     variant_classification %in% c("Missense_Mutation", 
    #                                                   "Nonsense_Mutations", "Silent",
    #                                                   "Frame_Shift_Del",
    #                                                   "Frame_Shift_Ins", "In_Frame_Del", 
    #                                                   "In_Frame_Ins", "Indel"))
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
      ## ## Set_1 are smokers/male/EurAmr, Set_2 are non-smokers/female/AfrAmr
      
      ## These dataframes will be used to graph
      
      Set_1 <- df %>% 
        filter(status == cats[1], pipeline == pipelines[j], variant_classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Silent", "5'UTR", "3'UTR")) %>% 
        group_by(variant_classification) %>% summarise(n = n()) %>% mutate(perc = n/sum(n)) %>%
        ungroup()
      Set_1$status <- cats[1]
      
      Set_2 <- df %>% 
        filter(status == cats[2], pipeline == pipelines[j], variant_classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Silent", "5'UTR", "3'UTR")) %>% 
        group_by(variant_classification) %>% summarise(n = n()) %>% mutate(perc = n/sum(n)) %>%
        ungroup()
      Set_2$status <- cats[2]
      
      combined_set <- rbind(Set_1, Set_2)
      
      
      
      ## This will be used to calculate significance
      range_vc <- df %>%
        filter(variant_classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Silent", "5'UTR", "3'UTR"), pipeline == pipelines[j]) %>% 
        group_by(status, tumor_barcode) %>% count(variant_classification) %>% ungroup()
      
      
      vars <- unique(range_vc$variant_classification)
      pVal_df <- tribble()
      for(k in 1:length(vars)){
        tryCatch({
          temp <- range_vc %>% filter(variant_classification == vars[k])
          pValue <- wilcox.test(n~status, data = temp, paired = F)$p.value
          temp_df <- tribble(~variant_classification, ~pValue, ~Cancer, ~Pipeline,
                             paste0(vars[k]), pValue, paste0(Data_Names[i]), paste0(pipelines[j]))
          pVal_df <- rbind(pVal_df, temp_df)
        }, error = function(e){})
      }
      
      pVal_df <- pVal_df %>% mutate(sig = ifelse(pVal_df$pValue < 0.05, "1", "0"))

      pVal_df$loc <- 1
      
      var_graph <- ggplot(combined_set) + geom_bar(aes(x = variant_classification, y = perc, fill = status), 
                                                   stat = "identity", position = "dodge") +
        coord_flip() + geom_point(data = pVal_df %>% filter(sig == 1), aes(x = variant_classification, y = loc),
                                  shape = "*", size=5, show.legend = FALSE,
                                  color = "black") 
      
      
      setwd(paste0("~/Research/BiolinksAnalysis/Output/", category, "/Graphs/var_class_Graph/", pipelines[j], "/"))
      ggsave(var_graph, file = file.path(paste0(pipelines[j], "_", Data_Names[i], "_var_classification.jpg"))
             , width = 6, height = 6, units = "in")
      
      
      setwd(paste0("~/Research/BiolinksAnalysis/Output/", category, "/pValues/var_class_pval/", pipelines[j], "/"))
      write_delim(pVal_df, paste0(pipelines[j], "_", Data_Names[i], "_var_class_pvalues", ".tsv"), delim = "\t")
      setwd("~/Research/BiolinksAnalysis/")
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

