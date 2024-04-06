#! /python3

# 2023-11-14 (modified 23-11-23)
# input: out/csv file (from PRIME_specify_hla.sh)
# output: a csv file which has binding affinity scores for each (specified) HLA for a list of desired variants
# author: Barbara Walkowiak bw450

# %% 
# to run: 
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/PRIME_out/csv
# python3 ../../scripts/PRIME_analysis1.py /Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/PRIME_out/csv/csv_20240220

# %%
import sys
import pandas as pd
import os

import time
timestr = time.strftime("%Y%m%d") # get current date 
# print (timestr) 

directory = sys.argv[1] # supply the directory when you are running the script

df_all = pd.DataFrame() # create an empty dataframe to append to at the end of the loop

# quick note on what I am selecting
# for each HLA allele, PRIME outputs 3 parameters: %rank, score, %rank_binding 
# %Rank indicates the fraction of random peptides from the human proteome (length 8 to 14) that would have a score higher or equal to the peptide given in input (does this mean I should extend my k-mer script to 14?)

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".csv"): 

        # read in the file 
        # print(os.path.join(directory, filename)) # sanity check
        csv = pd.read_csv(file)
        
        # get the gene and variant data first 
        gene = filename.split('_')[2]
        variant = filename.split('_')[3] # what HOTSPOT we are looking at
        genotype = filename.split('_')[4] 
        # do we have the CH variant or the reference seq corresponding to this variant (eg what is normally there in the wt)
        # values of genotype are either "refseq" or "ch"

        # Looking at the %Rank for EL 
        rank = '%Rank_'
        ranks = csv.filter(regex=rank, axis=1)
        stacked_ranks = ranks.stack().reset_index()
        stacked_ranks.columns = ['index', 'allele', '%Rank']
        stacked_ranks['allele'] = stacked_ranks['allele'].str.split('_').str[1]

        # remove the "bestAllele" thing, it is really annoying
        mask = stacked_ranks['allele'] != 'bestAllele'
        stacked_ranks = stacked_ranks[mask]
    
        min_values_rank = stacked_ranks.groupby('allele')['%Rank'].min().reset_index()
        sum_05_rank = stacked_ranks.groupby('allele')['%Rank'].apply(lambda x: (x < 0.5).sum()).reset_index(name='sum_peptides_below_05')
        sum_2_rank = stacked_ranks.groupby('allele')['%Rank'].apply(lambda x: (x < 2).sum()).reset_index(name='sum_peptides_below_2')

        # merge dfs by shared column (allele)
        merged_ranks = pd.merge(sum_05_rank, sum_2_rank, on='allele', how="inner")
        merged_ranks = pd.merge(min_values_rank, merged_ranks, on='allele', how="inner")

        merged_ranks['gene'] = gene
        merged_ranks['variant'] = variant 
        merged_ranks['genotype'] = genotype
        
        # make allele names consistent across different methdos 
        string_to_add = 'HLA-' # add "HLA-" to make consistent with other lists  
        merged_ranks['allele'] = string_to_add + merged_ranks['allele']
        character_to_add = ':'
        merged_ranks['allele'] = merged_ranks['allele'].apply(lambda x: x[:-2] + character_to_add + x[-2:])

        # print(merged_ranks)
        merged_ranks.columns = ['allele', 'min_rank', 'sum_peptides_below_05', 'sum_peptides_below_2', 'gene', 'variant', 'genotype']
        
        df_all = pd.concat([df_all, merged_ranks])

    else:
         continue

# # write the complete file to csv
df_all.to_csv("../../scores/" + timestr + "_percent_ranks_for_each_variant_by_HLA.csv", index = False)
