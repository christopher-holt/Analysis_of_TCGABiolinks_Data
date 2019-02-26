#!/usr/bin/env python
# coding: utf-8
#!/usr/bin/python3
import pandas
import os
import glob
os.chdir("../Datasets")
filenames = (glob.glob("*[!_select]*.csv")) 

## For system memory purposes, it may be easier to run the script in jupyter and rerun the last cell as range(0,1) then range(1,2) etc through range(11,12) instead
## of one big loop
for file in range(len(filenames)):
    df = pandas.read_csv(filenames[file], sep="\t", low_memory = False)
    df1 = df[["tumor_barcode", "primary_diagnosis", "tumor_stage", "gender", "ethnicity", "race", "age_at_diagnosis", 
                       "pipeline", "Mut_Type", "nucChange", "hugo_symbol", "disease", "years_smoked", "cigarettes_per_day",
                       "chromosome", "consequence", "one_consequence", "gene", "exon_number", "tissue_or_organ_of_origin",
                       "variant_classification", "mutation_status"]]
    
    file_title = str(df.disease[2]+"_select.csv")
    df1.to_csv(file_title, sep = "\t")

