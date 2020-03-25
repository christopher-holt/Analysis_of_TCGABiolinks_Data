##-----------------------------------------------------------------------------------------------------------
## This file will generate all general functions need throughout this project 
##-----------------------------------------------------------------------------------------------------------


##-----------------------------------------------------------------------------------------------------------
## imports
##-----------------------------------------------------------------------------------------------------------

library(tidyverse)



##-----------------------------------------------------------------------------------------------------------
## Downloads the clinical data
##-----------------------------------------------------------------------------------------------------------

download_clinical <- function(Data){
  clinical <- GDCquery_clinic(paste0("TCGA-",Data))
  clinical$tumor_barcode <- str_replace_all(clinical$submitter_id, pattern = '-',replacement = '.')
  clinical$tumor_barcode <- str_sub(clinical$tumor_barcode, start = 1, end = 12)
  clinical$tumor_barcode <- sapply(clinical$tumor_barcode, tolower)
  clinical$submitter_id <- NULL
  colnames(clinical) <- tolower(colnames(clinical))
  
  return(clinical)
}

##-----------------------------------------------------------------------------------------------------------
## Dowloads and merges the mutational data from each of the four pipelines
## as well as formatting two columns to remove characters
##-----------------------------------------------------------------------------------------------------------

download_mutational<- function(file){
  pipelines <- c("somaticsniper")
  combined <- data.frame()
  
  for (i in 1:length(pipelines)){
    Data <- GDCquery_Maf(paste0(file), pipelines = pipelines[i])
    Data$tumor_barcode <- str_replace_all(Data$Tumor_Sample_Barcode, pattern = '-',replacement = '.')
    Data$tumor_barcode <- str_sub(Data$tumor_barcode, start = 1, end = 12)
    Data$tumor_barcode <- sapply(Data$tumor_barcode, tolower)
    Data$pipeline = paste0(pipelines[i])
    colnames(Data) <- tolower(colnames(Data))


    hgvsc1 <- str_replace_all(Data$hgvsc, pattern = "^c.", replacement = "")
    c <- as.data.frame(hgvsc1)
    c <- type.convert(c, as.is = T)
    Data_1 <- cbind(Data, c) 
    
    combined <- rbind(combined, Data_1)
    
  }
  return(combined)
}



##-----------------------------------------------------------------------------------------------------------
## This will take two columns containing nucleotide changes and create their own column
## as well as classify the nucleotide changes ad transitions and transversion
##-----------------------------------------------------------------------------------------------------------

classify_Changes <- function(Data){
  ti_c = c("^.*G>A$", "^.*A>G$", "^.*C>T$", "^.*T>C$")
  tv_c = c("^.*C>A$", "^.*C>G$", "^.*A>C$", "^.*G>C$", "^.*G>T$", "^.*T>G$", "^.*T>A$","^.*A>T$")
  
  Data <- Data %>% mutate(mut_type = ifelse(grepl(paste(ti_c, collapse = '|'), Data$hgvsc1) , "Ti", 
                                                              ifelse(grepl(paste(tv_c, collapse = '|'), Data$hgvsc1), "Tv",
                                                                     ifelse(grepl("^.*ins.*$", Data$hgvsc1), "insertion",
                                                                            ifelse(grepl("^.*del.*$", Data$hgvsc1), "deletion", 
                                                                                   ifelse(grepl("^.$", Data$hgvsc1), "None",
                                                                                          ifelse(is.na(Data$hgvsc1), NA, "other")))))))
  
  
  Data <- Data %>% mutate(nucchange = ifelse(grepl("^.*G>A$",Data$hgvsc1) , "G > A", 
                                                               ifelse(grepl("^.*A>G$",Data$hgvsc1), "A > G", 
                                                                      ifelse(grepl("^.*C>T$",Data$hgvsc1), "C > T",
                                                                             ifelse(grepl("^.*T>C$",Data$hgvsc1), "T > C",
                                                                                    ifelse(grepl("^.*C>A$",Data$hgvsc1), "C > A" ,
                                                                                           ifelse(grepl("^.*C>G$",Data$hgvsc1), "C > G", 
                                                                                                  ifelse(grepl("^.*A>C$",Data$hgvsc1), "A > C",
                                                                                                         ifelse(grepl("^.*G>C$",Data$hgvsc1), "G > C",
                                                                                                                ifelse(grepl("^.*G>T$",Data$hgvsc1), "G > T",
                                                                                                                       ifelse(grepl("^.*T>G$",Data$hgvsc1), "T > G",
                                                                                                                              ifelse(grepl("^.*T>A$",Data$hgvsc1), "T > A",
                                                                                                                                     ifelse(grepl("^.*A>T$",Data$hgvsc1), "A > T",
                                                                                                                                            ifelse(grepl("deletion",Data$hgvsc1), "deletion",
                                                                                                                                                   ifelse(grepl("insertion", Data$hgvsc1), "insertion", "other" )))))))))))))))
  
  
  
  
}


##-----------------------------------------------------------------------------------------------------------
## This function will take a dataframe and return the pValues for variant_classification and 
## nucchange categories, and Ti vs Tv, chromosome
##-----------------------------------------------------------------------------------------------------------
calc_pValues <- function(dataframe, col){
  df1 <- dataframe %>% rename("type" = paste0(col))
  main_df <- tribble(~cohort, ~type, ~pValue)
  main_mut <- unique(df1$type)
  canc <- unique(df1$cohort)
  ## This loop takes a while to run
  for(k in 1:length(canc)){
    for (j in 1:length(main_mut)){
      tryCatch({
      temp <- df1 %>% filter(cohort == canc[k], type == main_mut[j])
      pValue <- wilcox.test(n~sex, data = temp,paired = FALSE)$p.value
      temp_df <- tribble(~cohort, ~type, ~pValue,
                         paste0(canc[k]), paste0(main_mut[j]), pValue)
      
      main_df <- rbind(main_df, temp_df)
      }, error = function(e){})
    }
  }
  
  main_df <- main_df %>% mutate(sig = ifelse(main_df$pValue < 0.05, 1, 0)) 
  if(col == "variant_classification"){
    main_df1 <- main_df %>% rename("variant_classification" = "type")
  } else if(col == "nucchange"){
    main_df1 <- main_df %>% rename("nucchange" = "type")
  } else if(col == "mut_type"){
    main_df1 <- main_df %>% rename("mut_type" = "type")
  } else if(col == "chromosome"){
    main_df1 <- main_df %>% rename("chromosome" = "type")
  }
  
  return(main_df1)
}





##-----------------------------------------------------------------------------------------------------------
## This function will take a dataset determine the gene frequencies
##-----------------------------------------------------------------------------------------------------------

gene_count <- function(df){
  new_df <- df %>% select(hugo_symbol) %>% 
    group_by(hugo_symbol) %>% count()
  
  return(new_df)
}

diff_genes <- function(df, df2){
  df3 <- as.data.frame(setdiff(df$hugo_symbol, df2$hugo_symbol)) %>%
    rename("hugo_symbol" = `setdiff(df$hugo_symbol, df2$hugo_symbol)`) %>% mutate_if(is.factor, as.character)
  
  return(df3)
}

gene_freq <- function(df){
  main_cancers <- tibble()
  temp_cancers <- unique(df$disease)
  for (i in 1:length(temp_cancers)){
    temp_df <- df %>% filter(disease == temp_cancers[i])

    ## Set 1 is male, set 2 is female
    Set_1 <- df %>% filter(sex == "male")
    Set_2 <- df %>% filter(sex == "female")
    
    
    
    ## for each person, count the number of times a gene is mutated per person
    genes_1 <- as.data.frame(gene_count(Set_1)) 
    genes_2 <- as.data.frame(gene_count(Set_2))
    
  
    ## Determine the genes in common
    same_genes <- as.data.frame(intersect(genes_1$hugo_symbol, genes_2$hugo_symbol)) %>% 
      rename("hugo_symbol" = `intersect(genes_1$hugo_symbol, genes_2$hugo_symbol)`) %>% mutate_if(is.factor, as.character)
    
    diff_1 <- diff_genes(genes_1, genes_2) ## Genes in genes_1 that are not in genes_2
    diff_2 <-diff_genes(genes_2, genes_1) ## Genes in genes_2 that are not in genes_1
    
    diff_1$n <- as.integer(0)
    diff_2$n <- as.integer(0)
    
    ## number of people who have a mutation in that gene

    genes_1 <- genes_1 %>% select(hugo_symbol, n)
    genes_2 <- genes_2 %>% select(hugo_symbol, n)
    
    genes_1 <- rbind(genes_1, diff_2)
    genes_2 <- rbind(genes_2, diff_1)
    
    
    genes_1$status <- "male"
    genes_2$status <- "female"
    

    
    ## n.x is the number of mutations a specific person has for gene
    ## n.y is the number of people in the whole population who have >=1 mutation in that gene

    genes <- rbind(genes_1, genes_2)

    
    genes$disease <- temp_cancers[i]
    
    main_cancers <- rbind(main_cancers, genes)
  }
  
  return(main_cancers)

  }


##-----------------------------------------------------------------------------------------------------------
## This function will sort the hclusted dendograms of a matrix row and column
##-----------------------------------------------------------------------------------------------------------

sort_hclust <- function(matrix){
  as.hclust(dendsort(as.dendrogram(matrix)))
}

##-----------------------------------------------------------------------------------------------------------
## This function will save the pheatmap as a png
##-----------------------------------------------------------------------------------------------------------


pheatsave <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  
  dev.off()
}


##-----------------------------------------------------------------------------------------------------------
## End of Script
##-----------------------------------------------------------------------------------------------------------
















