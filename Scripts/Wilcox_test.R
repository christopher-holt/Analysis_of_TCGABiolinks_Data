library(tidyverse)

HNSC <- read.csv("~/Desktop/Biolinks/HNSC_biolinks.csv")
cancers <- as.character(unique(HNSC$disease))
site <- as.character(unique(HNSC$tissue_or_organ_of_origin))
## "Cheek mucosa" does not have enough AA for the t.test
## "Base of tongue, NOS" has no AA in set
## "Posterior wall of oropharynx" has no AA in set
## Supraglottis has no AA in set
## Hypopharynx, NOS does not have enough x values for t.test
site <- site[site != "Cheek mucosa"]
site <- site[site != "Base of tongue, NOS"]
site <- site[site != "Posterior wall of oropharynx"]
site <- site[site != "Supraglottis"]
site <- site[site != "Gum, NOS"]
site <- site[site != "Hypopharynx, NOS"]
site <- site[site != "Mouth, NOS"]
site <- site[site != "Hard palate"]
site <- site[site != "Palate, NOS"]
site <- site[site != "Anterior floor of mouth"]
site <- site[site != "Mandible"]
site <- site[site != "Border of tongue"]
site <- site[site != "Lip, NOS"]
site <- site[site != "Pharynx, NOS"]
site <- site[site != "Upper Gum"]
site <- site[site != "Lower gum" ]
site <- site[site != "Oropharynx, NOS"]
pipelines = as.character(unique(HNSC$pipeline))
i = 1
j = 1
n = 1
k = 1
for(i in 1:length(cancers)){
  for (j in 1:length(pipelines)){
    for (n in 1:length(site)){
      Dataset_1 <- HNSC %>% filter(disease == cancers[i], mutation_status == "Somatic",
                                   variant_classification %in% c("Missense_Mutation", 
                                                                 "Nonsense_Mutations", "Silent",
                                                                 "Frame_Shift_Del",
                                                                 "Frame_Shift_Ins", "In_Frame_Del", 
                                                                 "In_Frame_Ins", "Indel"), tissue_or_organ_of_origin == "Oropharynx, NOS", pipeline == pipelines[j])
      
      
      ####### Dataset 1 analysis
      ## Distribution of point mutations
      
      EurAmr <- Dataset_1 %>% filter(race == "white") %>% mutate_if(is.factor,as.character)
      AfrAmr <- Dataset_1 %>% filter(race == "black or african american") %>% mutate_if(is.factor,as.character)
      
      EA_nuc <- EurAmr %>% filter(!(nucChange %in% c("deletion", "insertion", "other")) )  %>% mutate_if(is.factor,as.character)
      AA_nuc <- AfrAmr %>% filter(!(nucChange %in% c("deletion", "insertion", "other")) )  %>% mutate_if(is.factor,as.character)
      
      EA_nuc_table <- as.data.frame(table(EA_nuc$tumor_barcode,EA_nuc$nucChange)) %>% mutate_if(is.factor,as.character)
      AA_nuc_table <- as.data.frame(table(AA_nuc$tumor_barcode,AA_nuc$nucChange)) %>% mutate_if(is.factor,as.character)
      
      
      EA_nuc_table$race <- paste0(length(unique(EA_nuc_table$Var1)), "-white")
      AA_nuc_table$race <- paste0(length(unique(AA_nuc_table$Var1)), "-african_american")
      
      combined_table <- rbind(EA_nuc_table, AA_nuc_table)
      
      combined_table_1 <- combined_table %>% group_by(race, Var1, Var2) %>% summarise(total = sum(Freq))
      combined_table_1 <- combined_table_1 %>% group_by(Var1) %>% mutate(perc = (total/sum(total))*100)
      nucChanges <- as.character(unique(combined_table_1$Var2))
      for(k in 1:length(nucChanges)){
        change <- combined_table_1 %>% filter(Var2 == nucChanges[k])
        pValue <- wilcox.test(total~race, data = change, paired = F)$p.value
        pVal_pc <- wilcox.test(perc~race, data = change, paired = F)$p.value
        write.csv(pValue, file = file.path("/home/chris-holt/Desktop/Biolinks/Wilcox/pValues/",
                                           paste(cancers[i], '_', pipelines[j],'_', site[n],'_', nucChanges[k],'_pvalue', 
                                                 '.csv')), row.names = F)
        write.csv(pVal_pc, file = file.path("/home/chris-holt/Desktop/Biolinks/Wilcox/pValues/",
                                            paste(cancers[i], '_', pipelines[j],'_', site[n],'_', nucChanges[k],'_pvalue_pc', 
                                                  '.csv')), row.names = F)
        
        assign(paste(cancers[i], '_', pipelines[j],'_', nucChanges[k],'_pvalue',sep=''),pValue)
        assign(paste(cancers[i], '_', pipelines[j],'_', nucChanges[k],'_pvalue_pc',sep=''),pVal_pc)
      }
    }
    
  }
}
## Now for total num of mutations

for(i in 1:length(cancers)){
  for (j in 1:length(pipelines)){
    for (n in 1:length(site)){
      Dataset_1 <- HNSC %>% filter(disease == cancers[i], mutation_status == "Somatic",
                                   variant_classification %in% c("Missense_Mutation", 
                                                                 "Nonsense_Mutations", "Silent",
                                                                 "Frame_Shift_Del",
                                                                 "Frame_Shift_Ins", "In_Frame_Del", 
                                                                 "In_Frame_Ins", "Indel"), tissue_or_organ_of_origin == "Larynx, NOS", pipeline == pipelines[j])
      
      
      ####### Dataset 1 analysis
      ## Distribution of point mutations
      
      EurAmr <- Dataset_1 %>% filter(race == "white") %>% mutate_if(is.factor,as.character)
      AfrAmr <- Dataset_1 %>% filter(race == "black or african american") %>% mutate_if(is.factor,as.character)
      
      
      
      # EurAmr_table <- EurAmr %>% select(tumor_barcode, variant_classification) %>% group_by(tumor_barcode) %>% count(variant_classification)
      # AfrAmr_table <- AfrAmr %>% select(tumor_barcode, variant_classification) %>% group_by(tumor_barcode) %>% count(variant_classification)
      # 
      # names(EurAmr_table) <- c("Var1", "Var2", "Freq")
      # names(AfrAmr_table) <- c("Var1", "Var2", "Freq")
      # 
      # 
      
      EurAmr_table <- as.data.frame(table(EurAmr$tumor_barcode, EurAmr$variant_classification))
      AfrAmr_table <- as.data.frame(table(AfrAmr$tumor_barcode, AfrAmr$variant_classification))
      
      EurAmr_table$Var1 <- type.convert(EurAmr_table$Var1, as.is = T)
      EurAmr_table$Var2 <- type.convert(EurAmr_table$Var2, as.is = T)
      AfrAmr_table$Var1 <- type.convert(AfrAmr_table$Var1, as.is = T)
      AfrAmr_table$Var2 <- type.convert(AfrAmr_table$Var2, as.is = T)
      
      
      EurAmr_Sum <- EurAmr_table %>% group_by(Var1) %>% summarise(total = sum(Freq))
      AfrAmr_Sum <- AfrAmr_table %>% group_by(Var1) %>% summarise(total = sum(Freq))
      
      EurAmr_Sum$race <- paste0(length(unique(EurAmr_table$Var1)), "-white")
      AfrAmr_Sum$race <- paste0(length(unique(AfrAmr_table$Var1)), "-african_american")  
      
      Sum <- rbind(EurAmr_Sum, AfrAmr_Sum)
      
      pVal_Sum <- wilcox.test(total ~ race, data = Sum, paired = F)$p.value
      write.csv(pVal_Sum, file = file.path("/home/chris-holt/Desktop/Biolinks/Wilcox/pValues/Distro",
                                           paste(cancers[i], '_', pipelines[j],'_', site[n],'_pvalue_mut_distro', 
                                                 '.csv')), row.names = F)
      assign(paste(cancers[i], '_', pipelines[j],'_pvalue_mut_distro',sep=''),pVal_Sum)
      
      
    }
  }
  
}

## Now Ti vs TV


for(i in 1:length(cancers)){
  for (j in 1:length(pipelines)){
    for (n in 1:length(site)){
      Dataset_1 <- HNSC %>% filter(disease == cancers[i], mutation_status == "Somatic",
                                   variant_classification %in% c("Missense_Mutation", 
                                                                 "Nonsense_Mutations", "Silent",
                                                                 "Frame_Shift_Del",
                                                                 "Frame_Shift_Ins", "In_Frame_Del", 
                                                                 "In_Frame_Ins", "Indel"), tissue_or_organ_of_origin == "Larynx, NOS", pipeline == pipelines[j])
      
      
      ####### Dataset 1 analysis
      ## Distribution of point mutations
      
      EurAmr <- Dataset_1 %>% filter(race == "white") %>% mutate_if(is.factor,as.character)
      AfrAmr <- Dataset_1 %>% filter(race == "black or african american") %>% mutate_if(is.factor,as.character)
      
      Ti_Tv_EA <- EurAmr %>% filter(Mut_Type %in% c("Ti", "Tv")) %>% mutate_if(is.factor,as.character)
      Ti_Tv_AA <- AfrAmr %>% filter(Mut_Type %in% c("Ti", "Tv")) %>% mutate_if(is.factor,as.character)
      
      EA_table <- as.data.frame(table(Ti_Tv_EA$tumor_barcode,Ti_Tv_EA$Mut_Type)) %>% mutate_if(is.factor,as.character)
      AA_table <- as.data.frame(table(Ti_Tv_AA$tumor_barcode,Ti_Tv_AA$Mut_Type)) %>% mutate_if(is.factor,as.character)
      
      EA_Sum <- EA_table %>% group_by(Var1) %>% mutate(pc = Freq/sum(Freq)*100)
      AA_Sum <- AA_table %>% group_by(Var1) %>% mutate(pc = Freq/sum(Freq)*100)
      
      EA_Sum$race <- paste0(length(unique(EA_table$Var1)), "-white")
      AA_Sum$race <- paste0(length(unique(AA_table$Var1)), "-african_american")
      
      combined_Ti_Tv <- rbind(EA_Sum, AA_Sum)
      type <- unique(combined_Ti_Tv$Var2)
      for(l in 1:length(type)){
        subset <- combined_Ti_Tv %>% filter(Var2 == type[l])
        pVal_Sum_ti_tv <- wilcox.test(Freq ~ race, data = subset, paired = F)$p.value
        pVal_Sum_ti_tv_pc <- wilcox.test(pc ~ race, data = subset, paired = F)$p.value
        write.csv(pVal_Sum_ti_tv, file = file.path("/home/chris-holt/Desktop/Biolinks/Wilcox/pValues/Ti_TV",
                                                   paste(cancers[i], '_', pipelines[j],'_', site[n],'_', type[k], '_pvalue_TiTv', 
                                                         '.csv')), row.names = F)
        write.csv(pVal_Sum_ti_tv_pc, file = file.path("/home/chris-holt/Desktop/Biolinks/Wilcox/pValues/Ti_TV",
                                                   paste(cancers[i], '_', pipelines[j],'_', site[n],'_', type[k], '_pvalue_TiTv_pc', 
                                                         '.csv')), row.names = F)
        assign(paste(cancers[i], '_', pipelines[j], '_', type[k], '_pvalue_mut_distro',sep=''),pVal_Sum_ti_tv)
        assign(paste(cancers[i], '_', pipelines[j], '_', type[k], '_pvalue_mut_distro',sep=''),pVal_Sum_ti_tv_pc)
      }
    }
  }
  
}
