##-----------------------------------------------------------------------------------------------------------
## This file will compare African-American vs EurAmr data from 10 cancer datasets
##-----------------------------------------------------------------------------------------------------------

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
##
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
    df_1 <- df %>% filter(pipeline == pipelines[j])
    valid_sites <- race_sites(df_1)
    if(length(valid_sites) > 0){
      for(k in 1:length(valid_sites)){
        ## Subset the data into AfrAmr and EurAmr. These will be in a loop that separates cancer, pipeline, and location
        ## Set_1 are AfrAmr, Set_2 are EurAmr
        Set_1 <- df %>% filter(race == "black or african american", pipeline == pipelines[j], tissue_or_organ_of_origin == valid_sites[k])
        Set_2 <- df %>% filter(race == "white", pipeline == pipelines[j], tissue_or_organ_of_origin == valid_sites[k])
        
        
        ##-----------------------------------------------------------------------------------------------------------
        ## For each Cancer, determine the valid sites, then for each location and each pipeline,
        ## calculate the pValue for the distribution of all somatic point mutations
        ##-----------------------------------------------------------------------------------------------------------
        
        ## For each dataset, find the total number of somatic point mutations per person 
        Set_2_Sum <- number_mut(Set_2, "white")
        Set_1_Sum <- number_mut(Set_1, "black or african american")
        Sum <- rbind(Set_2_Sum, Set_1_Sum)
        total_mut <- ggplot(Sum) + geom_boxplot(aes(x = status, y = total))
        
        ## Calculate the pValue and quartile data for both groups
        tryCatch({
          pVal_Sum <- wilcox.test(total ~ status, data = Sum, paired = F)$p.value
          Overview_Set_2 <- summary_data(Set_2_Sum)
          Overview_Set_1 <- summary_data(Set_1_Sum)
          
          ## Combine summary data
          Overview_Set_2$status <- "white"
          Overview_Set_1$status <- "black or african american"
          Overview <- rbind(Overview_Set_2, Overview_Set_1)
          Overview$pValue <- pVal_Sum
          
          ## Save summary data and plots
          write_delim(Overview, file.path("~/Research/BiolinksAnalysis/Output/Race/pValues/total_mut_pval/", paste0(pipelines[j], "/", 
                                                                                                                     Data_Names[i],"_", valid_sites[k], 
                                                                                                                     ".tsv")), delim  = "\t")
        }, error = function(e){})
        ggsave(total_mut, file = file.path("~/Research/BiolinksAnalysis/Output/Race/Graphs/total_mut_Graph/", paste0(pipelines[j], "/", 
                                                                                                                      Data_Names[i], "_", valid_sites[k], 
                                                                                                                      ".jpg")), width = 6,
               height = 6, units = "in")
        
        ##-----------------------------------------------------------------------------------------------------------
        ## For each Cancer, determine the valid sites, then for each location and each pipeline,
        ## calculate the pValue between the percentage of nucleotide changes
        ##-----------------------------------------------------------------------------------------------------------
        Set_1_nuc <- nucChange_Sum(Set_1, "black or african american", "AfrAmr")
        Set_2_nuc <- nucChange_Sum(Set_2, "white", "EurAmr")
        
        combined_nuc <- rbind(Set_1_nuc, Set_2_nuc)
        combined_nuc_1 <- combined_nuc %>%  group_by(status, Var1, Var2, abbr) %>% 
          summarise(total = sum(Freq)) %>% group_by(Var1) %>% mutate(perc = (total/sum(total))*100)
        
        
        combined_nuc_1 <- combined_nuc_1 %>% group_by(Var1) %>% mutate(perc = (total/sum(total))*100)
        nucChanges <- as.character(unique(combined_nuc_1$Var2))
        
        nuc_graph <- ggplot(combined_nuc_1) + geom_boxplot(aes(x = Var2, y = perc, fill = status)) + coord_flip()
        
        ggsave(nuc_graph, file = file.path("~/Research/BiolinksAnalysis/Output/Race/Graphs/nucChange_Graph/", paste0(pipelines[j], "/", 
                                                                                                                      Data_Names[i], "_", valid_sites[k], 
                                                                                                                      ".jpg")), width = 6,
               height = 6, units = "in")
        for(n in 1:length(nucChanges)){
          change <- combined_nuc_1 %>% filter(Var2 == nucChanges[n])
          tryCatch({
            pValue <- wilcox.test(total~status, data = change, paired = F)$p.value
            
            
            change_1_summary <- change %>% filter(abbr == "AfrAmr")
            change_1_summary_data <- summary_data(change_1_summary)
            change_1_summary_data$status <- "AfrAmr"
            
            change_2_summary <- change %>% filter(abbr == "EurAmr")
            change_2_summary_data <- summary_data(change_2_summary)
            change_2_summary_data$status <- "EurAmr"
            
            change_summary <- rbind(change_1_summary_data, change_2_summary_data)
            change_summary$pValue <- pValue
            write_delim(change_summary, file.path("~/Research/BiolinksAnalysis/Output/Race/pValues/nucChange_pVal/", paste0(pipelines[j], "/", 
                                                                                                                             Data_Names[i],"_", valid_sites[k],
                                                                                                                             "_", nucChanges[n],
                                                                                                                             ".tsv")), delim  = "\t")
            
          }, error = function(e){})
          
        }
        
        ##-----------------------------------------------------------------------------------------------------------
        ## For each Cancer, determine the valid sites, then for each location and each pipeline,
        ## calculate the pValue  for Ti vs Tv
        ##-----------------------------------------------------------------------------------------------------------
        
        Set_1_TiTv <- TiTv_count(Set_1, "black or african american", "AfrAmr")
        Set_2_TiTv <- TiTv_count(Set_2, "white", "EurAmr")
        
        combined_Ti_Tv <- rbind(Set_1_TiTv, Set_2_TiTv)
        
        TiTv <- ggplot(combined_Ti_Tv) + geom_boxplot(aes(x = status, y = pc, fill = Var2))
        ggsave(TiTv, file = file.path("~/Research/BiolinksAnalysis/Output/Race/Graphs/TiTv_Graph/", paste0(pipelines[j], "/", 
                                                                                                            Data_Names[i], "_", valid_sites[k], 
                                                                                                            ".jpg")), width = 6,
               height = 6, units = "in")
        
        
        combined_Ti_Tv <- combined_Ti_Tv %>% rename("total" = "pc")
        type <- unique(combined_Ti_Tv$Var2)
        
        for(l in 1:length(type)){
          type_data <- combined_Ti_Tv %>% filter(Var2 == type[l])
          tryCatch({
            type_pval <- wilcox.test(total ~ status, data = type_data, paired = F)$p.value
            
            type_data_1 <- type_data %>% filter(abbr == "AfrAmr")
            type_data_1_summary <- summary_data(type_data_1)
            type_data_1_summary$status <- "AfrAmr"
            
            type_data_2 <- type_data %>% filter(abbr == "EurAmr")
            type_data_2_summary <- summary_data(type_data_2)
            type_data_2_summary$status <- "EurAmr"
            
            type_data_combined <- rbind(type_data_1_summary,type_data_2_summary)
            
            type_data_combined$pvalue <- type_pval
            
            write_delim(change_summary, file.path("~/Research/BiolinksAnalysis/Output/Race/pValues/TiTv_pVal/", paste0(pipelines[j], "/", 
                                                                                                                        Data_Names[i],"_", valid_sites[k],
                                                                                                                        "_", type[l],
                                                                                                                        ".tsv")), delim  = "\t")
          }, error = function(e){})
        }
      }
    }else{
      
    }
  }
}















