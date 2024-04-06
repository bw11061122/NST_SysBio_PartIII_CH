#! /python3

# 2023-11-20
# input: data/UKBB/ukb_typed_hla_with_ids_07_threshold_processed.txt
# new scores dataset which compares mutant and refseq and takes different measures of difference 
# this is bc there are many reports that distance b/n mut and wt (differential )

# %% 
# to run: 
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA
# python3 scripts/add_dia_indexes_to_scores.py /Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/PRIME_out/scores/20231116_percent_ranks_for_each_variant_by_HLA.csv
import pandas as pd
import sys 
import time

timestr = time.strftime("%Y%m%d") # get current date 

# scores_file = "/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/PRIME_out/scores/20240103_percent_ranks_for_each_variant_by_HLA.csv"
scores_file = sys.argv[1]

path = scores_file.split('/')[:8]
path = '/'.join(map(str, path))
print(path)
filename = scores_file.split('/')[8] # get file name but without the csv extension
filename = filename.split('.')[0]
print("The scores file examined is:", filename, "in", path)

scores = pd.read_csv(scores_file)
# add expression data 

# %%
# add a new column based on comparison of variant / ref-seq for the same gene
# I am creating two DAI indexes which were discussed in https://www.biorxiv.org/content/10.1101/2022.03.14.484285v1.full 
# DAI 1 = wt - mutant (HIGHER score indicates high immunogenicity: means that wt > mutant so mutant is more immunogenic)
# DAI 2 = mutant / wt (LOWER score indicates higher immunogenicity)
# note that you can try doing this for peptide sum but most peptides will have 0 so that's gonna be very fun to divide by

def wt_minus_mut(row, param, df):
    if row['variant'] == "W288fs":
        pass
    else:
        if row['genotype'] == 'refseq':
            return 0
        else:
            filter_condition = (df['genotype'] == 'refseq') & (df['gene'] == row['gene']) & (df['variant'] == row['variant']) & (df['allele'] == row['allele'])
            return df.loc[filter_condition, param].values[0] - row[param]

def ratio_mut_wt(row, param, df):
    if row['variant'] == "W288fs":
        pass
    else:
        if row['genotype'] == 'refseq':
            return 1
        else:
            filter_condition = (df['genotype'] == 'refseq') & (df['gene'] == row['gene']) & (df['variant'] == row['variant']) & (df['allele'] == row['allele'])
            return row[param] / (df.loc[filter_condition, param].values[0])

scores['min_rank_wt_minus_mut'] = scores.apply(wt_minus_mut, param = "min_rank", df = scores, axis=1)
scores['min_rank_ratio_mut_wt'] = scores.apply(ratio_mut_wt, param = "min_rank", df = scores, axis=1)

# %%
# take log of minimum rank here to make life easier for yourself 
# THIS IS -log10(minimum EL %rank)
import numpy as np
scores["min_rank_log"] = -1 * np.log10(scores['min_rank'])

# now doing the log difference
# log(x / y) = log(x) - log(y)
# so you want to have raw mutant over wt ratio and then log this
# e.g., mut = 0.0001, wt = 0.1, -log(mut) = 3, -log(wt) = 1, mut/wt = 0.01, -log(mut/wt) = 2
# name the column like this bc it is log()
scores['min_rank_ratio_mut_wt_log'] = -1 * np.log10(scores['min_rank_ratio_mut_wt'])

# %%
# # add expression data from GTEx 
# # not considering expression levels because I assume they are the same for everyone - I have no way of knowing if they are not
# tpm_expression = '~/Desktop/msc_thesis/task1_predict_binding_to_HLA/data/ch_var_test_gtex_tpm.csv'
# tpm_expression = pd.read_csv(tpm_expression)
# tpm_var = tpm_expression[['gene', 'GTEX_medianTPM_whole_blood']].drop_duplicates()

# scores = pd.merge(scores,tpm_var, on='gene')
# print(scores.head)

# scores['consider_based_on_tpm_gtex'] = [True if i > 1 else False for i in scores['GTEX_medianTPM_whole_blood']]

# # U2AF1 is expressed at VERY low levels 
# # Interestingly, DNMT3A is ALSO expressed at very low levels
# # Note sure what to make of this at the moment

# %%
# # add expression data from CD34+ HSC bulk RNA-seq data
# cd34_expression = '~/Desktop/msc_thesis/task1_predict_binding_to_HLA/data/ch_var_test_cd34_hsc.csv'
# cd34_expression = pd.read_csv(cd34_expression)

# scores = pd.merge(scores, cd34_expression, on='gene')
# print(scores.head)

# scores['consider_based_on_counts_cd34'] = [True if i > 5 else False for i in scores['mean_logcounts_cd34_hsc']]

# %%
# save the file with added columns to a new file 
scores.to_csv(path + '/' + filename + "_dia_added.csv") 
