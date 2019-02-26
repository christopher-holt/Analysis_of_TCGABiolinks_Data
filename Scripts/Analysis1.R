##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

library(tidyverse)
library(data.table)
setwd("~/Research/BiolinksAnalysis/Datasets")
##-----------------------------------------------------------------------------------------------------------
## Read in the data
##-----------------------------------------------------------------------------------------------------------
Data_Names <- c("LAML", "HNSC")#, "BLCA", "ESCA", "KICH", "KIRC", "LUSC", "PAAD", "STAD", "LIHC", "KIRP", "LUAD")

## This will read in each Cancer and determine the important sites
for (i in 1:length(Data_Names)){
  df <- read_delim(paste0(Data_Names[i], "_select.csv"), delim = "\t")
  assign(paste(Data_Names[i]), df)
  rm(df)
}

##-----------------------------------------------------------------------------------------------------------
## Subset the data into smokers and non-smokers. These will be in a loop that separates cancer, pipeline, and 
## location
##-----------------------------------------------------------------------------------------------------------
not_smoke <- HNSC %>% filter(is.na(cigarettes_per_day))
smoke <- HNSC %>% filter(!(is.na(cigarettes_per_day)))


##-----------------------------------------------------------------------------------------------------------
## For each dataset, find the total number of somatic point mutations per person 
##-----------------------------------------------------------------------------------------------------------
smoke_Sum <- number_mut(smoke, "smoker")
not_smoke_Sum <- number_mut(not_smoke, "non-smoker")
Sum <- rbind(smoke_Sum, not_smoke_Sum)



ggplot(Sum) + geom_boxplot(aes(x = status, y = total)) + ylim(0,1000)

pVal_Sum <- wilcox.test(total ~ status, data = Sum, paired = F)$p.value
write.csv(pVal_Sum, file = file.path("/home/chris-holt/Research/BiolinksAnalysis/Output/Smoke/total_mut_pval/",
                                     paste(cancers[i], '_', pipelines[j],'_', site[n],'_pvalue', 
                                           '.csv')), row.names = F)

  
summary_smoke <- as.data.frame(do.call(cbind, lapply(smoke_Sum, summary)))
summary_not_smoke <- as.data.frame(do.call(cbind, lapply(not_smoke_Sum, summary)))
summary_smoke <- add_rownames(summary_smoke, "Values")
summary_not_smoke <- add_rownames(summary_not_smoke, "Values")
Overview <- summary_smoke %>% select(Values, total)
Overview$pValue <- pVal_Sum
  

























