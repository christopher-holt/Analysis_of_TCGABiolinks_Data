## Author: Chris Holt
## Purpose: create heatmap of genes that are significant in all cancers
## Date Created: 18/Apr/2019
## Date of last update: 18/Apr/2019

source("Scripts/functions.R")


rm(list = ls())
Data_Names <- c("BLCA", "HNSC", "KICH", "LUAD", "LUSC", "PAAD", "KIRP")
## Read in Data
load_data <- function(Data_Names, category){
  temp_df <- tibble()
  pipelines <- c("muse", "mutect", "somaticsniper", "varscan2")
  setwd(paste0("~/Research/BiolinksAnalysis/Output/",category,"/Genes_Pvalues/"))
  for (i in 1:length(Data_Names)){
    for (j in 1:length(pipelines)){
      df <- read_delim(paste0(Data_Names[i], "_", 
                            pipelines[j], "_FINAL_PVALUES.tsv"), delim = "\t")
      df$cohort <- Data_Names[i]
      df$pipeline <- pipelines[j]
      
      temp_df <- rbind(temp_df, df)
    }
  }
  return(temp_df)
}

main_df <- load_data(Data_Names, "Gender")

### Subset somatic sniper data
main_df_som <- main_df %>% filter(pipeline == "somaticsniper") %>% mutate(sig = ifelse(pvalue < 0.05, "1", "0"))

## Find the intersecting genes
## There is no overlap
cancers <- unique(main_df_som$cohort)
for (i in 1:length(cancers)){
  temp <- main_df_som %>% filter(cohort == cancers[i], sig == 1)
  if(nrow(temp) > 0){
    assign(cancers[i], temp)
  }
}
genes <- unique(main_df_som$gene)
sig_genes = c()
for (i in 1:length(genes)){
  if(genes[i] %in% BLCA$gene){
    if(genes[i] %in% HNSC$gene){
      if(genes[i] %in% KIRP$gene){
        if(genes[i] %in% LUAD$gene){
          if(genes[i] %in% LUSC$gene){
            sig_genes <- sig_genes %>% append(genes[i])
            
          }
        }
      }
    } 
  }
}



cancers <- unique(main_df_som$cohort)
for (i in 1:length(cancers)){
  temp <- main_df_som %>% filter(cohort == cancers[i], sig == 1)
  if(nrow(temp) > 0){
    assign(cancers[i], temp)
  }
}





ggplot(main_df_som %>% filter(sig == 1)) + geom_tile(aes(cohort, gene, fill = sig)) + theme(axis.text.y = element_blank())

