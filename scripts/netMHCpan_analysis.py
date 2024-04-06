#! /python3

# 2023-11-08
# cd Desktop/msc_thesis/task1_predict_binding_to_HLA
# input: out/xls file (from netMHCpan_specify_hla.sh)
# output: a csv file which has binding affinity scores for each (specified) HLA for a list of desired variants
# author: Barbara Walkowiak bw450

# %% 
# to run: 
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/xls_sub
# python3 ../../scripts/netMHCpan_analysis.py /Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/xls_sub

# %%
import xlrd 
import sys
import pandas as pd
import os

import time
timestr = time.strftime("%Y%m%d") # get current date 
# print (timestr) 

directory = sys.argv[1] # supply the directory when you are running the script

dfs_all = pd.DataFrame() # create an empty dataframe to append to at the end of the loop

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".xlsx"): 

        # read in the file 
        print(os.path.join(directory, filename))
        xlsx = pd.read_excel(file, header=[0, 1])

        # for each HLA, calculate:
        # sum of peptides with score < 0.5
        # sum of peptides with 0.5 =< score < 2
        # extract best (ie lowest) Rank EL score # this is what I understood is what I want 

        # select EL rank
        EL_rank = xlsx.xs('EL_Rank', level=1, axis=1) # select EL rank scores for all HLAs
        
        min_EL_rank = EL_rank.min(numeric_only=True) # get the minimum score for a variant 
        # print("the min BA score is", min_BA_score, "for each HLA")

        sum_EL_rank_05 = (EL_rank.iloc[:,:] < 0.5).sum() # get the sum of peptides w/ score < 0.5 for each HLA
        # print("The sum of peptides with the BA score < 0.5 is:", sum_BA_score_05, "for each HLA")

        sum_EL_rank_05_2 = (EL_rank.iloc[:,:] < 2).sum() - (EL_rank.iloc[:,:] < 0.5).sum()
        # get the sum of peptides with scores below 2 - scores below 0.5 (> this gives scores between 0.5 and 2)
        # print("The sum of peptides with the BA score between 0.5 and 2 is:", sum_BA_score_05_2, "for each HLA")

        sum_EL_rank_2 =(EL_rank.iloc[:,:] < 2).sum()
        # I am adding this column bc more bio relevant to ask about all binders than just weak ones

        # okay, so now I have a df that has the peptide, its BA score for each HLAs examined, and indication of the gene / varinat this peptide is derived from
        # now, I want to get out summary stats (nr peptides w/ scores < 0.5, scores between 0.5 and 2, and best rank)

        to_bind = min_EL_rank, sum_EL_rank_05, sum_EL_rank_2 # including this not weak binders bc all binders makes more sense to me?
        df = pd.concat(to_bind, axis=1)
        
        # add data on the gene and variant analysed 
        gene = filename.split('_')[2]
        variant = filename.split('_')[3]
        genotype = filename.split('_')[4]

        df["gene"] = gene # add column with gene 
        df["variant"] = variant # add column with variant 
        df['genotype'] = genotype

        # make sure that allele goes to a column not rownames 
        df.index.name = 'allele'
        df.reset_index(inplace=True)

        df.columns = ['allele', 'min_rank', 'sum_peptides_below_05', 'sum_peptides_below_2', 'gene', 'variant', 'genotype']
    
        dfs_all = pd.concat([dfs_all, df])
    
    else:
         continue

# once this is done for every file, concatenate (rbind) files to get one dataframe

dfs_all.to_csv("../../scores/" + timestr + "_percent_ranks_for_each_variant_by_HLA.csv", index = False)

# desired output:
    # for each variant (e.g., DNMT3A R882H) we want to have a single score for a given HLA allele

# %%
