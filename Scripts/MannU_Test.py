#!/usr/bin/python3
# coding: utf-8

"""
@author: Chris Holt

This file will take flat files containing frequencies of genes and calculate the pValues using mannwhitneyu comparison 

Date Created: 15/Mar/2019

Date of Last Update: 26/Mar/2019

"""

## Imports
import pandas
from scipy import stats
import os
import glob


## Function that will calculate pValues 
def pvalues(string, disease):
    os.chdir("/home/chris-holt/Research/BiolinksAnalysis/Output/%s/Genes_Pvalues" %(str(string)))
    filenames = glob.glob("*_%s*[!_FINAL_PVALUES]*.tsv" % disease)
    for file in range(len(filenames)):
        df = pandas.read_table(filenames[file], sep = "\t", low_memory = False)
        genes = df["hugo_symbol"].unique()
        final_df = pandas.DataFrame(columns = ["gene", "pvalue"])
        for k in range(len(genes)):
            try:
                df1 = df.loc[df.hugo_symbol == genes[k]]
                df1_1 = df1.loc[df1.status == "NS"]
                df1_2 = df1.loc[df1.status == "S"]
                pval = stats.mannwhitneyu(df1_1["perc"], df1_2["perc"])
                data = [[genes[k], pval.pvalue]]
                b = pandas.DataFrame(data, columns = ["gene","pvalue"])
                final_df = pandas.concat([final_df, b])
            except ValueError:
                pass
    
        name = str(filenames[file].split("_")[1]+ "_" + filenames[file].split("_")[0] + "_FINAL_PVALUES")
        final_df.to_csv(name, sep = "\t")
        print(name + " download complete")
    
 
def main():
    
    ## I have broken down the process to be easier on the computer processor
    ## Enter in a valid directory and cancer cohort
    group = str(input("Please enter a directory (Smoke/Race/Gender): " ))
    Cancer = str(input("Please enter a Cancer: "))
    pvalues(group, Cancer)
    
    ## Valid directories:
        ## Smoke
        ## Race
        ## Gender
    ## Valid Cohorts:
        ## BLCA 
        ## HNSC
        ## KICH
        ## KIRP
        ## LUAD
        ## LUSC
        ## PAAD
    ## (Race/Gender Only)
        ## KIRC, LIHC, STAD

## Takes ~15min to run for one combination

## Completed combos(Smoke:BLCA)


main()