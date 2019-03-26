#!/usr/bin/env python
# coding: utf-8

import pandas
from scipy import stats
import os
import glob



def pvalues(string):
    os.chdir("/home/chris-holt/Research/BiolinksAnalysis/Output/%s/Genes_Pvalues" %(str(string)))
    filenames = glob.glob("*.tsv")
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
                continue
    
    name = str(filenames[file].split("_")[1]+ "_" + filenames[file].split("_")[0] + "FINAL_PVALUES")
    final_df.to_csv(name, sep = "\t")
    print(name + " download complete")
    
 
def main():
    pvalues("Smoke")
    pvalues("Race")
    pvalues("Gender")

main()