library(tidyverse)
library(magrittr)
##-----------------------------------------------------------------------------------------------------------
## Does each cancer site have enough data? This will only consider somatic point mutations
##-----------------------------------------------------------------------------------------------------------

sites <- function(Data, category){
  pot_sites <- as.character(unique(Data$tissue_or_organ_of_origin))
  new_site <- character()
  
  if ("white" %in% Data$category){
    Data %<>% filter(race %in% c("white", "black or african american"))
    cats <- as.character(unique(Data$category))
  } else if(category == "gender"){
    cats <- as.character(unique(Data$category))
  }
  
  for (i in 1:length(pot_sites)){
    set1 <- Data %>% filter(mutation_status == "Somatic",
                            variant_classification %in% c("Missense_Mutation", 
                                                          "Nonsense_Mutations", "Silent",
                                                          "Frame_Shift_Del",
                                                          "Frame_Shift_Ins", "In_Frame_Del", 
                                                          "In_Frame_Ins", "Indel"), 
                            tissue_or_organ_of_origin == pot_sites[i])
    
    set1_a <- set1 %>% filter(category == cats[1])
    set1_b <- set1 %>% filter(category == cats[2])
    
    if((nrow(set1_b) > 0 & nrow(set1_a) > 0) & ((length(unique(set1_b$tumor_barcode)) >= 3 & length(unique(set1_a$tumor_barcode)) >= 3))){
      new_site <- c(new_site, pot_sites[i])
    } 
    
  }
  return(new_site)
}
##-----------------------------------------------------------------------------------------------------------
## Testing sites function
##-----------------------------------------------------------------------------------------------------------

a <- sites(HNSC, race)
b <- sites(HNSC, gender)


##-----------------------------------------------------------------------------------------------------------
## End of Script
##-----------------------------------------------------------------------------------------------------------
