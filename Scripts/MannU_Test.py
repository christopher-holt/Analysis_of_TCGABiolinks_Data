#!/usr/bin/python3
# coding: utf-8

"""
@author: Chris Holt

This file will take flat files containing frequencies of genes and calculate the pValues using mannwhitneyu comparison 



Date Created: 15/Mar/2019

Date of Last Update: """ + now.strftime("%d/%B/%y") + """

"""

## Imports
import pandas
import time
from scipy import stats
import os
import glob
import datetime


## Function that will calculate pValues 
def pvalues(string, pipeline):
    print(time.ctime())
    os.chdir("/home/chris-holt/Research/BiolinksAnalysis/Output/%s/Genes_Pvalues" %(str(string)))
    filenames = glob.glob("%s_*[!_FINAL_PVALUES]*.tsv" % pipeline)
    for file in range(len(filenames)):
        df = pandas.read_table(filenames[file], sep = "\t", low_memory = False)
        genes = df["hugo_symbol"].unique()
        final_df = pandas.DataFrame(columns = ["gene", "pvalue"])
        stat = df.status.unique()
        for k in range(len(genes)):
            try:
                df1 = df.loc[df.hugo_symbol == genes[k]]
                df1_1 = df1.loc[df1.status == stat[0]]
                df1_2 = df1.loc[df1.status == stat[1]]
                pval = stats.mannwhitneyu(df1_1["n.x"], df1_2["n.x"])
                data = [[genes[k], pval.pvalue]]
                b = pandas.DataFrame(data, columns = ["gene","pvalue"])
                final_df = pandas.concat([final_df, b])
            except ValueError:
                pass
    
        name = str(filenames[file].split("_")[1]+ "_" + filenames[file].split("_")[0] + "_FINAL_PVALUES")
        final_df.to_csv(name, sep = "\t")
        print(name + " download complete")
        print("Files left to be analysed: " + str(len(filenames) - file - 1))
        print("-----------------------------------------------------------")
        print(time.ctime())
 
def main():
    
    ## I have broken down the process to be easier on the computer processor

## somatic sniper pipelines
 #   pvalues("Smoke", "somaticsniper")
 #   pvalues("Race", "somaticsniper")
 #   pvalues("Gender", "somaticsniper")
    
    ## Muse pipeline
    
   # pvalues("Smoke", "muse")
    #pvalues("Race", "muse")
   # pvalues("Gender", "muse")
  
  ## Mutect
  
    #pvalues("Smoke", "mutect")
  #  pvalues("Race", "mutect")
  #  pvalues("Gender", "mutect")
  
  ## varscan2
    #pvalues("Smoke", "varscan2")
  #  pvalues("Race", "varscan2")
  
  ## Just need to run this command
  #  pvalues("Gender", "varscan2")
    

main()