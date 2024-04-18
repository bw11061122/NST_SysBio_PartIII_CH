# python file because jupyter cannot handle the large stuff I am doing 

# %%

# imported packages
import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import matplotlib.ticker as plticker
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.patches import Polygon
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import cm
import scipy.special
import scipy.integrate as it
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.stats import kde
import copy
import glob, os
import re
# from sklearn import datasets, linear_model
import pandas as pd
from decimal import *
from operator import itemgetter    
from collections import OrderedDict
import timeit
import time 
import csv
import seaborn as sns 
import scipy as sp
from sklearn.preprocessing import LabelEncoder
from matplotlib.ticker import LogLocator
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D

from scipy.stats import mannwhitneyu
from scipy.stats import kruskal
from myfunctions import transform_format


axisfont=11
titlefont=20
subtitlefont = 15
axislabelfont=12
legendfont = 11
tpfont = 12
plt.rcParams.update({'font.sans-serif':'Verdana'})

timestr = time.strftime("%Y%m%d") # get current date 

# %%

folder_path = '/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/data/ukb_hotspot_calls_annotated'  # Path to ukb files from Hamish 

# Get all files in the folder 
files = glob.glob(os.path.join(folder_path, '*.txt'))  
# print('Files examined are:', files)

# Initialize an empty dictionary 
dataframes = {}

# Iterate through the CSV files and read each one with pandas
for csv_file in files:
    
    df = pd.read_csv(csv_file, sep = '\t')
    dataframes[csv_file] = df

# Extract dataframes 
for file_name, df in dataframes.items():
    
    variable_name = file_name.split('/')[8].split('.')[0] + '_data'  # Removing the file extension
    print('Examined file:', variable_name)
    
    df['batch'] = variable_name # add column to indicate source 
    globals()[variable_name] = df  # assign dataframe to the variable 


# Concat all into one df
dfs_to_concat = [v for k, v in globals().items() if k.endswith('_data') and isinstance(v, pd.DataFrame)]

# Row bind all batch dataframes into a single one 
batch_all = pd.concat(dfs_to_concat, ignore_index=True)
batch_all = batch_all.dropna(subset=['batch']) # remove rows which are read incorrectly (w/o batch number)

print('Number of samples with variants examined:', batch_all.shape[0])


# Change numerical variables to integers

batch_all['end_position'] = batch_all['end_position'].astype(int)
batch_all['position'] = batch_all['position'].astype(int)
batch_all['sample_ID'] = batch_all['sample_ID'].astype(int)

# Subset and create new useful columns

batch_all = batch_all[['sample_ID', 'chromosome', 'end_position', 'VAF', 'var_depth', 'depth', 'Amino_acids', 'SYMBOL', 'Codons', 'batch']]
batch_all['alt_variant'] = batch_all['Amino_acids'].str.split('/', expand = True)[1] # alternative (CH) variant
batch_all['ref_variant'] = batch_all['Amino_acids'].str.split('/', expand = True)[0] # reference variant 

# There are some cases where there is no change in amino acids, for now save as NaN 

batch_all['alt_variant'].fillna(batch_all['Amino_acids'], inplace=True)
batch_all['ref_variant'].fillna(batch_all['Amino_acids'], inplace=True)

# Exclude data with singletons (likely errors)

batch_all_ns = batch_all[batch_all['var_depth'] >= 2]

# Number of samples with two variant reads or more 
print('Number of samples carrying more than a single read with the variant sequence:', batch_all_ns.shape[0])
# okay so we have 2823 of these but then note that only 2249 are with CH variants 
# so like 20% of >= 2 reads are not v likely to be real? should conslut with Jamie 

batch_all_ns.head(n = 10)

# Import indexes tested in each batch 
folder_path = '/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/data/ukb_hotspot_calls/batch_ids'  # Path to ukb files from Hamish 

# Get all files in the folder 
files_ids = glob.glob(os.path.join(folder_path, '*.tsv'))  

# Read each file one by one 

indexes = {}

# Iterate through the CSV files and read each one with pandas
for file in files_ids:
    
    id = pd.read_csv(file, sep = '\t')
    id = id.rename(columns={'batch ID': 'sample_ID'})
    id['sample_ID'] = id['sample_ID'].str.split('_', n = 1).str[0]
    indexes[file] = id

for file_name, df in indexes.items():
    
    variable_name = file_name.split('/')[9].split('.')[0]    # Remove file extension
    globals()[variable_name] = df  # Assign the DataFrame to a variable with the file name

# %%
    
# Find out how many gene_variants were called:

# don't show warnings 
import warnings
warnings.filterwarnings("ignore")

# identify variants called 
batch_all_ns['variant_coord'] = batch_all_ns['chromosome'].astype(str) + "_" + batch_all_ns['end_position'].astype(str) # specific position in the genome 
batch_all_ns['variant_coord'] = batch_all_ns['variant_coord'].astype('category')
batch_all_ns['variant_coord_pos'] = batch_all_ns['variant_coord'].astype(str) +  "_" + batch_all_ns['SYMBOL'].astype(str) + "_" + batch_all_ns['ref_variant'].astype(str) + "_" + batch_all_ns['alt_variant'].astype(str) # change to a specific aa
batch_all_ns['variant_coord_pos'] = batch_all_ns['variant_coord_pos'].astype('category')

# remove samples that have not been annotated (you can tell from the coordinate what is likely but these could be different mutations)
batch_all_ns = batch_all_ns.dropna(subset=['SYMBOL']) # remove column where gene is not known 
print('Number of samples which have been correctly annotated:', batch_all_ns.shape[0]) # but at each site, you are getting reads modified to sth else 
# ok so sth worked wrong with annotation in only 5 cases > that looks good 

# identify the number of variants in a specific position
num_variants = pd.DataFrame(batch_all_ns['variant_coord_pos'].value_counts())
num_variants = num_variants[num_variants['count']!=0]
num_variants = num_variants.sort_values(by = 'count')
num_variants['variant_coord_pos'] = num_variants.index
num_variants = num_variants.reset_index(drop=True)
num_variants = num_variants.sort_values(by = 'variant_coord_pos')
print('Number of variants identified in batches analysed:', num_variants.shape[0]) # but at each site, you are getting reads modiifed to sth else 

# identify the number of positions we looked at
num_sites = pd.DataFrame(batch_all_ns['variant_coord'].value_counts())
num_sites = num_sites[num_sites['count']!=0]
num_sites = num_sites.sort_values(by = 'count')
num_sites['variant_coord'] = num_sites.index
num_sites = num_sites.reset_index(drop=True)
num_sites = num_sites.sort_values(by = 'variant_coord')
print('Number of sites identified in batches analysed:', num_sites.shape[0]) # okay see so you only found 37 sites 

# we have more variants than sites because for each site, we look at any change from the reference sequence (any possible variant)
# at the same time, because we have this data, we can check if we have more variants that "random" mutations

# %%

# Annotation (using df with variant names and genomic coordinates)

# read in the df with coordinates
coord_gene_var = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/data/ch_variants_coordinates_tp53_added_nts.csv')

# all coordinates identified in the batches  
coord_out = num_sites['variant_coord'].tolist()

# intersection (annotate)
coord_gene_var['variant_coord'] = coord_gene_var['chromosome'] + "_" + coord_gene_var['end_position'].astype(str) # find cariant coordinates
coord_gene_var['SYMBOL'] = coord_gene_var['gene_var'].str.split('_').str[0] # find gene analysed
coord_gene_var['ref_variant'] = coord_gene_var['gene_var'].str.split('_').str[1].str[0] # reference sequence variant 
coord_gene_var['alt_variant'] = coord_gene_var['gene_var'].str.split('_').str[1].str[-1] # alternative sequence (CH / mutation) variant
coord_gene_var['variant_coord_pos'] = coord_gene_var['variant_coord'].astype(str) + "_" + coord_gene_var['SYMBOL'].astype(str) + "_" + coord_gene_var['ref_variant'].astype(str) + "_" + coord_gene_var['alt_variant'].astype(str) # specific mutation 

print('Number of variants which have been investigated:', len(coord_gene_var.variant_coord_pos.unique()))
print('Number of sites which have been investigated:',len(coord_gene_var.variant_coord.unique()))

# NB we removed one variant bc the coordinates were incorrect

# NOTE: 
# TP53_R249S was searched for but not found in any of the batches so far 
# there have been 5 samples identified (in DNMT3A, at 2 different sites) with variants that I did not originally search for 
# in addition, some of the variants were offset by one base (could have been a deletion / insertion)
# this is why we have differences in the number of sites in the two datasets 

# %%

# identify CH variants typed in the dataset 
gene_vars = coord_gene_var[['variant_coord_pos', 'gene_var']]

# subset batch df to only include variants which were successfully identified 
batch_gene_vars = pd.merge(batch_all_ns, gene_vars, on = 'variant_coord_pos', how = 'inner')
batch_gene_vars['gene_var'] = batch_gene_vars['gene_var'].astype('category')
batch_gene_vars['gene_var'].value_counts()
batch_gene_vars.rename(columns = {'sample_ID' : 'Person_ID'}, inplace = True)

gene_vars_count = pd.DataFrame(batch_gene_vars['gene_var'].value_counts())
gene_vars_count['gene_var'] = gene_vars_count.index
gene_vars_count = gene_vars_count.reset_index(drop=True)
gene_vars_sorted = gene_vars_count.sort_values(by = 'count', ascending=False)

print('Number of variants identified with annotations:', gene_vars_sorted.shape[0])
# so we will have 30 variants to look at bc the rest was not mapped correctly 

print('Number of samples with a mutation in a CH-relevant position (with >2 reads):', batch_all_ns.shape[0])
print('Number of patients with a mutation in a CH-relevant position (with >2 reads):', len(batch_all_ns.sample_ID.unique()))
print('Number of samples with a CH hotspot mutation (with > 2 reads):', batch_gene_vars.shape[0])
print('Number of patients with a CH hotspot mutation (with > 2 reads):', len(batch_gene_vars.Person_ID.unique()))


# Filtering 2: we will only be looking at variants for which we have more than 10 samples 
gene_vars_count = pd.DataFrame(batch_gene_vars['gene_var'].value_counts())
gene_vars_count['gene_var'] = gene_vars_count.index
gene_vars_count = gene_vars_count.reset_index(drop=True)
gene_vars_sorted = gene_vars_count.sort_values(by = 'count', ascending=False)

variants_to_examine = gene_vars_sorted[gene_vars_sorted['count'] >= 15].gene_var.tolist()
print('Number of variants to examine (min 10 samples with variant identified):', len(variants_to_examine))

# filter the dataframe to only include CH variants with minimum 10 samples 
# we will be using this dataset going forward 
batch_gene_vars_10 = batch_gene_vars[batch_gene_vars['gene_var'].isin(variants_to_examine)]
print('Number of samples with a CH hotspot mutation (with > 2 reads):', batch_gene_vars.shape[0])
print('Number of patients with a CH hotspot mutation (with > 2 reads):', len(batch_gene_vars.Person_ID.unique()))
print('Number of samples with a CH hotspot mutation (with > 2 reads) (common variants: >= 10 samples):', batch_gene_vars_10.shape[0])
print('Number of patients with a CH hotspot mutation (with > 2 reads) (common variants: >= 10 samples):', len(batch_gene_vars_10.Person_ID.unique()))

# %%

# load the age dataset 
age_data = pd.read_csv('~/Desktop/msc_thesis/task1_predict_binding_to_HLA/data/2022-05-20_healthy_pheno.tsv', sep = '\t')

age_df = age_data[['ID_v0', 'Age.when.attended.assessment.centre_v0']]
age_df.columns.values[0] = 'Person_ID'
age_df.columns.values[1] = 'age'
print('Number of individuals for whom age data is available:', age_df.shape[0])

# add age data to the rest of the data 
batch_gene_vars_10.gene_var = batch_gene_vars_10.gene_var.astype(str)
batch_gene_vars_10.gene_var = batch_gene_vars_10.gene_var.astype('category') # remove previous categories we will not be looking at
batch_gene_age = pd.merge(batch_gene_vars_10, age_df, on = 'Person_ID') # add age to everyone who has a variant
print("Number of sample with available CH variant:", batch_gene_vars_10.shape[0]) # note we are using the filtered dataframe here 
print("Number of sample with available CH variant and age:", batch_gene_age.shape[0]) # these are healthy / non-cancer cases (can think if we want to filter out the cancer ones, but tbf probably yes)
# note: this is throwing away patients who have been diagnosed with cancer (this dataset is for healthy individuals only)

# %%

# add MHC II genotype data to CH cases

# path to file
file_hla = "/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/data/UKBB/ukb_typed_hla_with_ids_07_threshold_processed.txt"

# get the header 
header = pd.read_csv(file_hla, sep='\t', nrows=1, header=None).values.tolist()
head = [item for sublist in header for item in sublist]

# get the actual dataframe 
df = pd.read_csv(file_hla, skiprows = 1, sep = ' ', header = None)

# add columns 
df.columns = head 

# subset data for HLA-I and HLA-II class alleles 
# we are only interested in HLA-I for the momnet 
df_hla1 = df.filter(regex='^(Person_|A_|B_|C_)') # 488377 cases 
df_hla2 = df.filter(regex='^(Person_|D)') # 488377 cases 

df_clean_hla2 = df_hla2[~df_hla2.isin([0.5]).any(axis=1)] # exclude people for whom we lack full genotype data
df_clean_hla2 = df_clean_hla2[~df_clean_hla2.isin([1.5]).any(axis=1)] # exclude people for whom we lack full genotype data

# we will analyse DR and DP / DQ alleles separately because the "rules" are different (ie alpha / beta chain combinations etc.)
# NOTE: I am removing samples for which we don't have correct allele annotation
df_hla2_dp = df_clean_hla2.filter(regex='^(Person_|DP)')
df_hla2_dq = df_clean_hla2.filter(regex='^(Person_|DQ)')
df_hla2_dr = df_clean_hla2.filter(regex='^(Person_|DR)')

print("Number of samples with MHC II genotype imputation data", df_hla2.shape[0])
print("Number of samples with MHC II genotype confidently imputed", df_clean_hla2.shape[0]) # not ideal as this really throws out a lot of the data we have actually 

# %%

# okay, a problem with DP / DQ is that we have predictions for COMBINATIONS
# therefore, we need to create columns with combinations and use them for when we make the prediction

# what we can do is create additional columns and set the value to 1 if someone has both alleles (who knows if they will have the combination
# but if predictions are available for it in NetMHCII we can assume that maybe it's legit)

from itertools import combinations

# combinations for DQ
cc = list(combinations(df_hla2_dq.columns[1:],2)) # possible combinations of these columns (first is Person_ID)
df_hla2_dq_comb = pd.concat([df_hla2_dq[c[0]].multiply(df_hla2_dq[c[1]]) for c in cc], axis=1, keys=cc) 
# we are multiplying
# in this way, if someone has >= 1 in both alleles they will get a positive number
# even if they have 2 alleles in one allele but 0 in the other, they will get a zero

df_hla2_dq_comb.columns = df_hla2_dq_comb.columns.map('-'.join) # join with a '-'

# repeat the same for DP 
cc = list(combinations(df_hla2_dp.columns[1:],2)) # possible combinations of these columns (first is Person_ID)
df_hla2_dp_comb = pd.concat([df_hla2_dp[c[0]].multiply(df_hla2_dp[c[1]]) for c in cc], axis=1, keys=cc) 
# we are multiplying
# in this way, if someone has >= 1 in both alleles they will get a positive number
# even if they have 2 alleles in one allele but 0 in the other, they will get a zero

df_hla2_dp_comb.columns = df_hla2_dp_comb.columns.map('-'.join) # join with a '-'

# %%
# okay, now: we do not want to predict stuff for 2 alpha or 2 beta chains because these are not created so we can remove these columns
# we also want to rename columns from XXB-XXA to XXA-XXB
def is_valid_column_name(col):
    parts = col.split('-') # two part of the name of the column
    return 'A' in parts[0] and 'B' in parts[1] or 'B' in parts[0] and 'A' in parts[1]

def switch_parts(col):
    parts = col.split('-') 
    new_name = parts[1] + '-' + parts[0]
    return new_name

# %%

# remove columns with BB or AA combinations 
filtered_columns_dp = [col for col in df_hla2_dp_comb.columns if is_valid_column_name(col)]
filtered_columns_dq = [col for col in df_hla2_dq_comb.columns if is_valid_column_name(col)]

# Create a new DataFrame with the filtered and reordered columns
new_df_dp = df_hla2_dp_comb[filtered_columns_dp]
new_df_dq = df_hla2_dq_comb[filtered_columns_dq]

renamed_columns_dp = [switch_parts(col) for col in new_df_dp.columns]
renamed_columns_dq = [switch_parts(col) for col in new_df_dq.columns]

new_df_dp.columns = renamed_columns_dp
new_df_dq.columns = renamed_columns_dq

# %%
# once we have selected the correct columns, combine with other data

# DP
df_hla2_dp_all = pd.concat([df_hla2_dp, new_df_dp], axis = 1)
df_hla2_dp_all.replace(2, 1, inplace=True) # if someone got a 2 (2 * 1), replace to 1
df_hla2_dp_all.replace(4, 1, inplace=True) # if someone got a 2 (2 * 2), replace to 1

# let's remove columns where noone has the combination (it's possible that some combinations are never seen in our dataset)
colsums = pd.DataFrame(df_hla2_dp_all.sum()).reset_index()
colsums.rename(columns={'index': 'col_name', 0:'col_sum'}, inplace=True)
col_to_retain = colsums[colsums['col_sum']>=1]['col_name'] # these are the columns to retain 
df_hla2_dp_all = df_hla2_dp_all[col_to_retain]
df_hla2_dp_all 

# DQ
df_hla2_dq_all = pd.concat([df_hla2_dq, new_df_dq], axis = 1)
df_hla2_dq_all.replace(2, 1, inplace=True) # if someone got a 2 (2 * 1), replace to 1
df_hla2_dq_all.replace(4, 1, inplace=True) # if someone got a 2 (2 * 1), replace to 1

# let's remove columns where noone has the combination (it's possible that some combinations are never seen in our dataset)
colsums = pd.DataFrame(df_hla2_dq_all.sum()).reset_index()
colsums.rename(columns={'index': 'col_name', 0:'col_sum'}, inplace=True)
col_to_retain = colsums[colsums['col_sum']>=1]['col_name'] # these are the columns to retain 
df_hla2_dq_all = df_hla2_dq_all[col_to_retain]
df_hla2_dq_all 

# %%

# parameters
param = '%Rank_EL' # let's do it on EL ranks first 

# NOTE: In the original file (sub), I was missing a score for DNMT3A_R736C_ch so I generated it separately with the same code and will concatenate the files here 
# load the file with binding predictions across variants
pred_file_dpq_sub = "/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/NetMHCII_out/scores/20240211_NetMHC_HLA_UKBB_with_affinities_DP_DQ_bestscores.csv"
pred_file_dpq_r736c = "/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/NetMHCII_out/scores/20240215_NetMHC_HLA_UKBB_with_affinities_DP_DQ_R736C_bestscores.csv"
pred_file_dpq_stop = "/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/NetMHCII_out/scores/20240214_NetMHC_HLA_UKBB_with_affinities_DP_DQ_stopcodons_bestscores.csv"
pred_method = pred_file_dpq_sub.split('_out')[0] # all are with the same method so does not matter which file you end up using 
 
# organize file with prediction scores for DP and DQ alleles 
pred_filename_dpq_sub = pred_file_dpq_sub.split('/')[2].split('.')[0]
pred_filename_dpq_stop = pred_file_dpq_stop.split('/')[2].split('.')[0]
pred_df_dpq_sub = pd.read_csv(pred_file_dpq_sub)
pred_df_dpq_r736c = pd.read_csv(pred_file_dpq_r736c)
pred_df_dpq_stop = pd.read_csv(pred_file_dpq_stop)

# In the sub file, I had some predictions for STOP-codon-containing variants, I don't want these so remove them
pred_df_dpq_sub = pred_df_dpq_sub[~pred_df_dpq_sub['gene_var'].str.contains('\*')]

# concat the files 
pred_df_dpq_stop.drop(columns=['Affinity (nM)'], inplace=True)
pred_df_dpq_r736c.drop(columns=['Affinity (nM)'], inplace=True)
pred_df_dpq = pd.concat([pred_df_dpq_sub, pred_df_dpq_stop, pred_df_dpq_r736c], ignore_index = True)
pred_df_dpq['gene_var_gt'] = pred_df_dpq['gene'] + '_' + pred_df_dpq['variant'] + '_' + pred_df_dpq['genotype'] # create a column that includes all genotype data

# organize file with prediction scores for DR alleles 
pred_file_dr_sub = "/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/NetMHCII_out/scores/20240214_NetMHC_HLA_UKBB_with_affinities_DR_bestscores.csv"
pred_file_dr_r736c = "/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/NetMHCII_out/scores/20240215_NetMHC_HLA_UKBB_with_affinities_DR_R736C_bestscores.csv"
pred_file_dr_stop = "/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/NetMHCII_out/scores/20240214_NetMHC_HLA_UKBB_with_affinities_DR_stopcodons_bestscores.csv"

# organize file with prediction scores for DR alleles 
pred_filename_dr_sub = pred_file_dr_sub.split('/')[2].split('.')[0]
pred_filename_dr_stop = pred_file_dr_stop.split('/')[2].split('.')[0]
pred_df_dr_sub = pd.read_csv(pred_file_dr_sub)
pred_df_dr_r736c = pd.read_csv(pred_file_dr_r736c)
pred_df_dr_stop = pd.read_csv(pred_file_dr_stop)

# In the sub file, I had some predictions for STOP-codon-containing variants, I don't want these so remove them
pred_df_dr_sub = pred_df_dr_sub[~pred_df_dr_sub['gene_var'].str.contains('\*')]

# remove affinity columns
pred_df_dr_stop.drop(columns=['Affinity (nM)'], inplace=True)
pred_df_dr_r736c.drop(columns=['Affinity(nM)'], inplace=True)

# concat the files 
pred_df_dr = pd.concat([pred_df_dr_sub, pred_df_dr_stop, pred_df_dr_r736c], ignore_index = True)
pred_df_dr['gene_var_gt'] = pred_df_dr['gene'] + '_' + pred_df_dr['variant'] + '_' + pred_df_dr['genotype'] # create a column that includes all genotype data
pred_df_dr.head(n = 5)

# %%

# we need to move the format of the DP/DQ predictions to something that matches our UKBB genotyping and everything else 
def transform_format_DPQ(input_string):
    # Define a regular expression pattern to match the input format
    pattern = re.compile(r'HLA-(\w{3})(\d{1})(\d{4})-(\w{3})(\d{1})(\d{4})') # okay so this is the pattern we are trying to match

    # check if there is a match
    match = pattern.match(input_string)

    # if match, apply transformation
    if match:
        group1 = match.group(1) # we are not including HLA annotations, this is the name of the first allele   
        group2 = int(match.group(2)) # there will be no zeroes, leave as it is 
        group3 = int(match.group(3)) # remove zeros at the start 
        group4 = match.group(4) # name of the second allele in the combination
        group5 = int(match.group(5)) # there will be no zeroes, leave as it is 
        group6 = int(match.group(6)) # remove zeros at the start 

        # Format the output string
        output_string = f'{group1}{group2}_{group3}-{group4}{group5}_{group6}' # stitch back 

        return output_string # return transformed string 

    # if no much, return original string 
    return 0

def transform_format_DR(input_string):
    # Define a regular expression pattern to match the input format
    pattern = re.compile(r'(\w{3})(\d{1})_(\d+)') # okay so this is the pattern we are trying to match

    # check if there is a match
    match = pattern.match(input_string)

    # if match, apply transformation
    if match:
        group1 = match.group(1) # we are not including HLA annotations, this is the name of the first allele   
        group2 = int(match.group(2)) # there will be no zeroes, leave as it is 
        group3 = int(match.group(3)) # remove zeros at the start 

        # Format the output string
        output_string = f'{group1}{group2}_{group3}' # stitch back 

        return output_string # return transformed string 

    # if no much, return original string 
    return 0

# %%
# change format of HLA allele naming 
pred_df_dpq['HLA_formatted'] = pred_df_dpq['HLA'].map(transform_format_DPQ)
pred_df_dr['HLA_formatted'] = pred_df_dr['HLA'].map(transform_format_DR)

# combine all predictions
pred_df_all = pd.concat([pred_df_dpq, pred_df_dr], axis = 0)
pred_df_all.head()


# now you need to merge the data that you want to use with the allele data to get predictions
# most comprehensive dataset we have is batch_gene_age (includes genetic variant someone carries and age)
# Find IDs of CH-affected persons and match to their HLA genotype 

batch_gene_vars.rename(columns={'sample_ID': 'Person_ID'}, inplace=True)
batch_gene_age.rename(columns={'sample_ID': 'Person_ID'}, inplace=True)
ids_batch = batch_gene_age['Person_ID']
hla_batch_ids_dr = df_hla2_dr[df_hla2_dr['Person_ID'].isin(ids_batch)] 

# add HLA genotype data to samples with annotated CH variant and age
batch_gene_age_hla_dr = pd.merge(batch_gene_age, df_hla2_dr, on='Person_ID')
batch_gene_age_hla_dp = pd.merge(batch_gene_age, df_hla2_dp_all, on='Person_ID')
batch_gene_age_hla_dq = pd.merge(batch_gene_age, df_hla2_dq_all, on='Person_ID')

batch_gene_hla_dr = pd.merge(batch_gene_vars, df_hla2_dr, on='Person_ID')
batch_gene_hla_dp = pd.merge(batch_gene_vars, df_hla2_dp_all, on='Person_ID')
batch_gene_hla_dq = pd.merge(batch_gene_vars, df_hla2_dq_all, on='Person_ID')

print("Number of samples with available CH variant:", batch_gene_vars.shape[0])
print('Number of samples with annotated variants CH variant and HLA genotype (DR):', batch_gene_hla_dr.shape[0]) # 1 HLA genotype missing 
print('Number of samples with annotated variants CH variant and HLA genotype (DP):', batch_gene_hla_dp.shape[0]) # 1 HLA genotype missing 
print('Number of samples with annotated variants CH variant and HLA genotype (DQ):', batch_gene_hla_dq.shape[0]) # 1 HLA genotype missing 
print("Number of samples with available CH variant and age:", batch_gene_age.shape[0])
print('Number of samples with annotated variants CH variant and age and HLA genotype (DR):', batch_gene_age_hla_dr.shape[0]) # 1 HLA genotype missing 
print('Number of samples with annotated variants CH variant and age and HLA genotype (DP):', batch_gene_age_hla_dp.shape[0]) # 1 HLA genotype missing 
print('Number of samples with annotated variants CH variant and age and HLA genotype (DQ):', batch_gene_age_hla_dq.shape[0]) # 1 HLA genotype missing 

# %%

# define function to find the best score for the variant that is carried (ie present in the person)
def find_best_score_for_variant_carried(row, df, param):

    '''
    This functiion is applied to the CH dataset
    The df is the dataset with predictions for given allele and genetic variant 
    parameter is what to base this prediction on (here will only be using %Rank_EL)
    '''
    
    # these are the values of rows to seach for 
    row_values = pd.to_numeric(row[1:-1], errors='coerce')
    
    # find alleles which are present 
    hla = row.index[1:-1][row_values >= 1]  
    
    # find values corresponding to these alleles 
    vals = df.loc[df['gene_var'] == row['gene_var'], hla].values.flatten() 
    
    # it can be that nothing was found e.g., bc there were no predictions made for this variant
    if vals.size == 0: 
        value = None
    
    else: 
    
        # note here we will only be using '%Rank_EL'
        if param == "Aff_nM":
            value = min(vals) # highest affinity corresponds to lowest value (see plots above)
        elif param == "%Rank_EL":
            value = min(vals) # the best rank is the lowest one (indicates peptide in top x % of binders)
        elif param == "Score_EL":
            value = max(vals) # the best score is the highest one 
        elif param == "%Rank_BA":
            value = min(vals) # the best rank is the lowest one (indicates peptide in top x % of binders)
        elif param == "Score_BA":
            value = max(vals) # the best score is the highest one
        else:
            print('Incorrect parameter provided') 
    
    return value  

# %%

# do this for DR alleles 
param = '%Rank_EL'

# DR
pred_sub_dr = pred_df_dr[['HLA_formatted', 'gene_var_gt', param]]
pred_sub_dr_wide = pd.pivot(pred_sub_dr, index='gene_var_gt', columns='HLA_formatted', values=param)
pred_sub_dr_wide = pred_sub_dr_wide.reset_index() # this is to make sure that you have the gene_var column in there too

hla_ukbb_dr = batch_gene_age_hla_dr.filter(regex='\d').columns.tolist() # relevant HLAs
hla_intersect_dr = pred_sub_dr_wide.columns[pred_sub_dr_wide.columns.isin(hla_ukbb_dr)] # HLA in the UKBB which I have predictions for 
hla_intersect_dr_list = hla_intersect_dr.tolist() 

# prepare gene variants names to match names in the Patient file 
pred_sub_dr = pred_sub_dr_wide[hla_intersect_dr_list + pred_sub_dr_wide.columns[pred_sub_dr_wide.columns.str.contains('gene_var')].tolist()]
pred_sub_dr = pred_sub_dr[pred_sub_dr['gene_var_gt'].str.contains('_ch', regex=True)] # retain CH scores only 
pred_sub_dr['gene_var'] = pred_sub_dr['gene_var_gt'].str.replace('_ch', '') # remove the ch / refseq annotation
pred_sub_dr['gene_var'] = pred_sub_dr['gene_var'].str.replace('_refseq', '') # remove refseq if present 
pred_sub_dr['gene_var'] = pred_sub_dr['gene_var'].astype(str)

# subset batch_gene_age_hla file 
ch_hla_sub_dr = batch_gene_age_hla_dr[hla_intersect_dr_list + batch_gene_age_hla_dr.columns[batch_gene_age_hla_dr.columns.str.contains('gene_var')].tolist()]
ch_hla_sub_dr = pd.concat([batch_gene_age_hla_dr["Person_ID"], ch_hla_sub_dr], axis=1) # add CH cases 
ch_hla_sub_dr['score'] = ch_hla_sub_dr.apply(find_best_score_for_variant_carried, df=pred_sub_dr, param=param, axis=1) # add score for the parameter
ch_hla_scores_dr = ch_hla_sub_dr.dropna() # remove NA (incorrectly annotated cases)

# merge scores with VAF and age
age_vaf_dr = batch_gene_age_hla_dr[['Person_ID', 'VAF', 'var_depth', 'age', 'gene_var']]
ch_hla_merge_dr = pd.merge(ch_hla_scores_dr, age_vaf_dr, on = ['Person_ID', 'gene_var'])

# now add the columns with VAF and age 
col_to_select = ['Person_ID', 'gene_var', 'score', 'age', 'var_depth', 'VAF'] # subset the data 
ch_hla_merge_sub_dr = ch_hla_merge_dr[col_to_select]
ch_hla_merge_sub_dr['log_score'] = -1*np.log10(ch_hla_merge_sub_dr['score']) # convert score to -log10(score)
ch_hla_merge_sub_dr['allele_type'] = 'DR'

# %%

# DP
pred_sub_dpq = pred_df_dpq[['HLA_formatted', 'gene_var_gt', param]]
pred_sub_dpq = pred_sub_dpq[~pred_sub_dpq.duplicated()] # remove duplicate columns 
pred_sub_dpq_reset = pred_sub_dpq.reset_index()
pred_sub_dpq_wide = pd.pivot(pred_sub_dpq_reset, index='gene_var_gt', columns='HLA_formatted', values=param)
pred_sub_dpq_wide = pred_sub_dpq_wide.reset_index() # this is to make sure that you have the gene_var column in there too

hla_ukbb_dp = batch_gene_age_hla_dp.filter(regex='\d').columns.tolist() # relevant HLAs
hla_intersect_dp = pred_sub_dpq_wide.columns[pred_sub_dpq_wide.columns.isin(hla_ukbb_dp)] # HLA in the UKBB which I have predictions for 
hla_intersect_dp_list = hla_intersect_dp.tolist() 

# prepare gene variants names to match names in the Patient file 
pred_sub_dp = pred_sub_dpq_wide[hla_intersect_dp_list + pred_sub_dpq_wide.columns[pred_sub_dpq_wide.columns.str.contains('gene_var')].tolist()]
pred_sub_dp = pred_sub_dp[pred_sub_dp['gene_var_gt'].str.contains('_ch', regex=True)] # retain CH scores only 
pred_sub_dp['gene_var'] = pred_sub_dp['gene_var_gt'].str.replace('_ch', '') # remove the ch / refseq annotation
pred_sub_dp['gene_var'] = pred_sub_dp['gene_var'].str.replace('_refseq', '') # remove refseq if present 
pred_sub_dp['gene_var'] = pred_sub_dp['gene_var'].astype(str)

# subset batch_gene_age_hla file 
ch_hla_sub_dp = batch_gene_age_hla_dp[hla_intersect_dp_list + batch_gene_age_hla_dp.columns[batch_gene_age_hla_dp.columns.str.contains('gene_var')].tolist()]
ch_hla_sub_dp = pd.concat([batch_gene_age_hla_dp["Person_ID"], ch_hla_sub_dp], axis=1) # add CH cases 
ch_hla_sub_dp['score'] = ch_hla_sub_dp.apply(find_best_score_for_variant_carried, df=pred_sub_dp, param=param, axis=1) # add score for the parameter
ch_hla_scores_dp = ch_hla_sub_dp.dropna() # remove NA (incorrectly annotated cases)

# merge scores with VAF and age
age_vaf_dp = batch_gene_age_hla_dp[['Person_ID', 'VAF', 'var_depth', 'age', 'gene_var']]
ch_hla_merge_dp = pd.merge(ch_hla_scores_dp, age_vaf_dp, on = ['Person_ID', 'gene_var'])

# now add the columns with VAF and age 
col_to_select = ['Person_ID', 'gene_var', 'score', 'age', 'var_depth', 'VAF'] # subset the data 
ch_hla_merge_sub_dp = ch_hla_merge_dp[col_to_select]
ch_hla_merge_sub_dp['log_score'] = -1*np.log10(ch_hla_merge_sub_dp['score']) # convert score to -log10(score)
ch_hla_merge_sub_dp['allele_type'] = 'DP'

# DQ
hla_ukbb_dq = batch_gene_age_hla_dq.filter(regex='\d').columns.tolist() # relevant HLAs
hla_intersect_dq = pred_sub_dpq_wide.columns[pred_sub_dpq_wide.columns.isin(hla_ukbb_dq)] # HLA in the UKBB which I have predictions for 
hla_intersect_dq_list = hla_intersect_dq.tolist() 

# prepare gene variants names to match names in the Patient file 
pred_sub_dq = pred_sub_dpq_wide[hla_intersect_dq_list + pred_sub_dpq_wide.columns[pred_sub_dpq_wide.columns.str.contains('gene_var')].tolist()]
pred_sub_dq = pred_sub_dq[pred_sub_dq['gene_var_gt'].str.contains('_ch', regex=True)] # retain CH scores only 
pred_sub_dq['gene_var'] = pred_sub_dq['gene_var_gt'].str.replace('_ch', '') # remove the ch / refseq annotation
pred_sub_dq['gene_var'] = pred_sub_dq['gene_var'].str.replace('_refseq', '') # remove refseq if present 
pred_sub_dq['gene_var'] = pred_sub_dq['gene_var'].astype(str)

# subset batch_gene_age_hla file 
ch_hla_sub_dq = batch_gene_age_hla_dq[hla_intersect_dq_list + batch_gene_age_hla_dq.columns[batch_gene_age_hla_dq.columns.str.contains('gene_var')].tolist()]
ch_hla_sub_dq = pd.concat([batch_gene_age_hla_dq["Person_ID"], ch_hla_sub_dq], axis=1) # add CH cases 
ch_hla_sub_dq['score'] = ch_hla_sub_dq.apply(find_best_score_for_variant_carried, df=pred_sub_dq, param=param, axis=1) # add score for the parameter
ch_hla_scores_dq = ch_hla_sub_dq.dropna() # remove NA (incorrectly annotated cases)

# merge scores with VAF and age
age_vaf_dq = batch_gene_age_hla_dq[['Person_ID', 'VAF', 'var_depth', 'age', 'gene_var']]
ch_hla_merge_dq = pd.merge(ch_hla_scores_dq, age_vaf_dq, on = ['Person_ID', 'gene_var'])

# now add the columns with VAF and age 
col_to_select = ['Person_ID', 'gene_var', 'score', 'age', 'var_depth', 'VAF'] # subset the data 
ch_hla_merge_sub_dq = ch_hla_merge_dq[col_to_select]
ch_hla_merge_sub_dq['log_score'] = -1*np.log10(ch_hla_merge_sub_dq['score']) # convert score to -log10(score)
ch_hla_merge_sub_dq['allele_type'] = 'DQ'

# we can actually combine all these scores now 
ch_hla_merge_sub_all = pd.concat([ch_hla_merge_sub_dr, ch_hla_merge_sub_dp, ch_hla_merge_sub_dq], axis = 0)

# filter so this only includes variants with >=10 samples present
gene_counts = ch_hla_merge_sub_all.gene_var.value_counts().reset_index()
variants_to_examine = gene_counts[gene_counts['count'] >= 30].gene_var.tolist()

ch_hla_merge_sub_all = ch_hla_merge_sub_all[ch_hla_merge_sub_all['gene_var'].isin(variants_to_examine)]
ch_hla_merge_sub_dr = ch_hla_merge_sub_dp[ch_hla_merge_sub_dp['gene_var'].isin(variants_to_examine)]
ch_hla_merge_sub_dp = ch_hla_merge_sub_dq[ch_hla_merge_sub_dq['gene_var'].isin(variants_to_examine)]
ch_hla_merge_sub_dq = ch_hla_merge_sub_dr[ch_hla_merge_sub_dr['gene_var'].isin(variants_to_examine)]
# NB in the cell paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6482006/) they don't consider stuff separately
# I would say let's do the analysis for 3 allele classes separately for now and we can then just pick the best if we want to 

# create a dataframe with all MHC class II genotype for each case
batch_gene_age_hla_all = pd.merge(batch_gene_age_hla_dq, batch_gene_age_hla_dp)
batch_gene_age_hla_all = pd.merge(batch_gene_age_hla_all, batch_gene_age_hla_dr)
batch_gene_age_hla_all.head()

# %%

# define function to find the best score for the variant that is carried (ie present in the person)
def find_best_score_for_variant_carried_all_hla(row, df, param):

    '''
    This functiion is applied to the CH dataset
    The df is the dataset with predictions for given allele and genetic variant 
    parameter is what to base this prediction on (here will only be using %Rank_EL)
    '''
    
    # these are the values of rows to seach for 
    row_values = pd.to_numeric(row[1:-1], errors='coerce')
    
    # find alleles which are present (NOTE this is changed from the previous function bc formatting is different)
    hla = row.index[1:-1][row_values >= 1]  
    
    # find values corresponding to these alleles 
    vals = df.loc[df['gene_var'] == row['gene_var'], hla].values.flatten() 
    
    # it can be that nothing was found e.g., bc there were no predictions made for this variant
    if vals.size == 0: 
        value = None
    
    else: 
    
        # note here we will only be using '%Rank_EL'
        if param == "Aff_nM":
            value = min(vals) # highest affinity corresponds to lowest value (see plots above)
        elif param == "%Rank_EL":
            value = min(vals) # the best rank is the lowest one (indicates peptide in top x % of binders)
        elif param == "Score_EL":
            value = max(vals) # the best score is the highest one 
        elif param == "%Rank_BA":
            value = min(vals) # the best rank is the lowest one (indicates peptide in top x % of binders)
        elif param == "Score_BA":
            value = max(vals) # the best score is the highest one
        else:
            print('Incorrect parameter provided') 
    
    return value  

# %%

# add a single score based on all of the MHC class II alleles someone has (not breaking down into DR, DP, DQ)
pred_sub_all = pred_df_all[['HLA_formatted', 'gene_var_gt', param]]
pred_sub_all = pred_sub_all[~pred_sub_all.duplicated()] # remove duplicate columns 
pred_sub_all_reset = pred_sub_all.reset_index()
pred_sub_all_wide = pd.pivot(pred_sub_all_reset, index='gene_var_gt', columns='HLA_formatted', values=param).reset_index()

hla_ukbb_all = batch_gene_age_hla_all.columns[16:].tolist() # this includes both separate DP / DQ and in combinations
hla_intersect_all = pred_sub_all_wide.columns[pred_sub_all_wide.columns.isin(hla_ukbb_all)] # HLA in the UKBB which I have predictions for 
hla_intersect_all_list = hla_intersect_all.tolist() 

# prepare gene variants names to match names in the Patient file 
pred_sub_all = pred_sub_all_wide[hla_intersect_all_list + pred_sub_all_wide.columns[pred_sub_all_wide.columns.str.contains('gene_var')].tolist()]
pred_sub_all = pred_sub_all[pred_sub_all['gene_var_gt'].str.contains('_ch', regex=True)] # retain CH scores only 
pred_sub_all['gene_var'] = pred_sub_all['gene_var_gt'].str.replace('_ch', '') # remove the ch / refseq annotation
pred_sub_all['gene_var'] = pred_sub_all['gene_var'].str.replace('_refseq', '') # remove refseq if present 
pred_sub_all['gene_var'] = pred_sub_all['gene_var'].astype(str)

# subset batch_gene_age_hla file 
ch_hla_sub_all = batch_gene_age_hla_all[hla_intersect_all_list + batch_gene_age_hla_all.columns[batch_gene_age_hla_all.columns.str.contains('gene_var')].tolist()]
ch_hla_sub_all = pd.concat([batch_gene_age_hla_all["Person_ID"], ch_hla_sub_all], axis=1) # add CH cases 
ch_hla_sub_all['score'] = ch_hla_sub_all.apply(find_best_score_for_variant_carried_all_hla, df=pred_sub_all, param=param, axis=1) # add score for the parameter
ch_hla_scores_all = ch_hla_sub_all.dropna() # remove NA (incorrectly annotated cases)

# merge scores with VAF and age
age_vaf_all = batch_gene_age_hla_all[['Person_ID', 'VAF', 'var_depth', 'age', 'gene_var']]
ch_hla_merge_all = pd.merge(ch_hla_scores_all, age_vaf_all, on = ['Person_ID', 'gene_var'])

# now add the columns with VAF and age 
col_to_select = ['Person_ID', 'gene_var', 'score', 'age', 'var_depth', 'VAF'] # subset the data 
ch_hla_merge_sub_all = ch_hla_merge_all[col_to_select]
ch_hla_merge_sub_all['log_score'] = -1*np.log10(ch_hla_merge_sub_all['score']) # convert score to -log10(score)

# filter out rare variants
ch_hla_merge_sub_all = ch_hla_merge_sub_all[ch_hla_merge_sub_all['gene_var'].isin(variants_to_examine)]

# %%

# get IDs of healthy individuals in the UKBB 
ch_sampleids = batch_gene_vars_10['Person_ID']
ukbb_no_ch = df_clean_hla2[~df_clean_hla2['Person_ID'].isin(ch_sampleids)]

df_hla2_dp_all_combinations = pd.concat([df_hla2_dp_all.iloc[:,0], df_hla2_dp_all.iloc[:,42:]], axis = 1)
df_hla2_dq_all_combinations = pd.concat([df_hla2_dq_all.iloc[:,0], df_hla2_dq_all.iloc[:,33:]], axis = 1)

# add their DP and DQ combinations
ukbb_no_ch_dp = pd.merge(ukbb_no_ch.iloc[:, :73], df_hla2_dp_all_combinations, on = 'Person_ID')
ukbb_no_ch_dp_dq = pd.merge(ukbb_no_ch_dp, df_hla2_dq_all_combinations, on = 'Person_ID')

# identify HLA alleles carried by non-CH-individuals in the UKBB data
hla_ukbb = ukbb_no_ch_dp_dq.filter(regex='\d').columns 

# prepare dataset with predictions (now just from NetMHC)
hla_ukbb_all = batch_gene_age_hla_all.columns[16:].tolist() # this includes both separate DP / DQ and in combinations
hla_intersect_all = pred_sub_all_wide.columns[pred_sub_all_wide.columns.isin(hla_ukbb_all)] # HLA in the UKBB which I have predictions for 
hla_intersect_all_list = hla_intersect_all.tolist() 
pred_sub_all_wide = pred_sub_all_wide[hla_intersect_all_list + pred_sub_all_wide.columns[pred_sub_all_wide.columns.str.contains('gene_var_gt')].tolist()]

ukbb_no_ch_hla = ukbb_no_ch_dp_dq[hla_intersect_all]
ukbb_no_ch_dp_dq = pd.concat([ukbb_no_ch_dp_dq['Person_ID'], ukbb_no_ch_hla], axis = 1)
ukbb_no_ch_dp_dq.head()

ids_non_carriers_healthy = ukbb_no_ch_dp_dq.Person_ID.unique().tolist()
print('Number of healthy (non-CH) individuals examined:', len(ids_non_carriers_healthy))


# select only HLAs for which we have predictions
df_clean_hla2_hla = pd.merge(df_clean_hla2.iloc[:,:73], df_hla2_dp_all_combinations, on = 'Person_ID')
df_clean_hla2_hla = pd.merge(df_clean_hla2_hla, df_hla2_dq_all_combinations, on = 'Person_ID')
df_clean_hla2_hla = df_clean_hla2_hla.drop('DRB3_9901', axis = 1)
df_clean_hla2_hla = df_clean_hla2_hla.drop('DRB4_9901', axis = 1)
df_clean_hla2_hla = df_clean_hla2_hla.drop('DRB5_9901', axis = 1)

# %%

# define new functions in case you want to get predictions for both CH and refseq

def find_best_score_all_variants_log_ref_ch(row, df, param):
    
    '''
    The same function but log scores 
    the only allowed parameters are %Rank_EL and %Rank_BA
    we want -1 * log(score) so the highest score if the most "immunogenic" (best binding) one
    '''
    
    hlas = row.index[1:-1][row[1:-1] >= 1] # select alleles which each Person (row) carries

    variants = df['gene_var_gt']
    
    scores = {} # initialise empty dictionaries

    if param == "%Rank_EL":
        for var in variants:

            # take maximum of the negative log score 
            best_value = max(-1 * np.log10(df.loc[df['gene_var_gt'] == var, hlas].values[0]))
            scores[f'score_{var}'] = best_value

        return pd.Series(scores)

    elif param == "%Rank_BA":
        for var in variants:

            # take maxium of the negative log score 
            best_value = max(-1 * np.log10(df.loc[df['gene_var_gt'] == var, hlas].values[0]))
            scores[f'score_{var}'] = best_value

        return pd.Series(scores)

    elif param == "Aff_nM":
        for var in variants:

            # take maxium of the negative log score (low aff = more immunogenic)
            best_value = max(-1 * np.log10(df.loc[df['gene_var_gt'] == var, hlas].values[0]))
            scores[f'score_{var}'] = best_value

        return pd.Series(scores)

# %%
    
# prepare to get the reference dataset 
np.random.seed(1) # ensure this is all reproducible 

# specify list of variants to look at 
variants = batch_gene_vars_10.gene_var.unique()
colors = ['#f00071', '#0422ed']
cols_ch = ['Person_ID', 'gene_var', 'score']
param = '%Rank_EL'

# initialize empty dataframe to store results in 
df_compare_to_ref = pd.DataFrame()

# %%
# make sure you only sample from people you screened and you know they don't have CH
batch_all['batch_number'] = batch_all['batch'].apply(lambda x: [int(num) for num in re.findall(r'\d+', x)])
batch_numbers_examined = set(batch_all['batch_number'].sum()) 
print('Batches that were examined for variants:', batch_numbers_examined)

all_dataframes = [f'batch_{i}_ids' for i in range(11, 61)]
selected_ids = []

# Iterate through all df name 
for df_name in all_dataframes:
    number = int(re.search(r'\d+', df_name).group())
    # check if number matches numbers of batches examined 
    if number in batch_numbers_examined:        
        selected_ids.append(globals()[df_name])

ids_examined = pd.concat(selected_ids, ignore_index=True)
print('Number of samples examined for variants:', ids_examined.shape[0])

ids_examined_list = ids_examined['sample_ID'].tolist()
ids_examined_list = [int(id) for id in ids_examined_list]
df_clean_hla2_hla = df_clean_hla2_hla[df_clean_hla2_hla['Person_ID'].isin(ids_examined_list)]

# %%
# we use 2k for the reference to have a really good idea of how the distribution looks like 
# we also do the plot for all variants together 

# full script to look across variants  

for var in variants:

    var_name = var.replace('_', ' ')
    
    # we will compare to %Rank EL (log) derived from NetMHC now
    ch_scores = ch_hla_merge_sub_all[ch_hla_merge_sub_all['gene_var'] == var]
    ch_scores = ch_scores[['Person_ID', 'gene_var', 'log_score']]
    ch_scores.rename(columns = {'log_score': 'score'}, inplace = True)

    # identify number of variant carriers 
    n =  ch_scores.shape[0]

    if n > 10: # we have filtered the dataset before so all variants should fulfill this condition 

        # find all the dataset that you can sample from
        # subset the dataframe to only include individuals who were actually screened for CH
        ukbb_not_carrier = df_clean_hla2_hla[df_clean_hla2_hla.Person_ID.isin(ids_examined_list)]
        ref = ukbb_not_carrier.sample(n=2000, replace=False) # no duplicates
            
        # specify parameter 
        param = '%Rank_EL'

        # add scores 
        ref_scores = pd.concat([ref, ref.apply(find_best_score_all_variants_log_ref_ch, df=pred_sub_all_wide, param=param, axis=1)], axis=1)
        ref_scores['gene_var'] = 'reference set'
        ref_scores.rename(columns={f'score_{var}_ch': 'score'}, inplace=True)
        ref_scores_var = ref_scores[cols_ch]
        
        # concatenate all dataframes
        scores_compare_all = pd.concat([ch_scores, ref_scores_var], axis=0)
        scores_compare_all["CH_variant_carrier"] = scores_compare_all['gene_var'] == var 
        scores_compare_all["CH_variant_carrier"] = scores_compare_all["CH_variant_carrier"].map({True: 'Carrier', False: 'Non-carrier'})
        
        # add rows to dataframe with data for all variants
        df_compare_to_ref = pd.concat([df_compare_to_ref, scores_compare_all], axis = 0)

# %%
# okay we need to add a colum with the variant for everyone (because you have reference sets to a variant but the gene_var column does not say which variant this is reference for)
labels = []
current_label = None

# Iterate through the DataFrame
for idx, row in df_compare_to_ref.iterrows():
    # If the category is 'ch', update the current label
    if row['CH_variant_carrier'] == 'Carrier':
        current_label = row['gene_var']
    # Append the current label to the list
    labels.append(current_label)

# Add the labels list as a new column to the DataFrame
df_compare_to_ref['CH variant'] = labels

# figure out desired order of variants on the plot 
df_compare_to_ref['CH variant2'] = df_compare_to_ref['CH variant'].str.replace('_', '\n')
df_compare_to_ref['gene_var2'] = df_compare_to_ref['gene_var'].str.replace('_', '\n')

# ordering based on scores in carriers only 
df_compare_to_ref_carr = df_compare_to_ref[df_compare_to_ref['CH_variant_carrier']=='Carrier']
df_compare_to_ref_carr['median_score'] = df_compare_to_ref_carr.groupby('CH variant2')['score'].transform('median')
df_compare_to_ref_sort = df_compare_to_ref_carr.sort_values(by='median_score', ascending = False)
order = df_compare_to_ref_sort.gene_var2.unique()

# %%

# plot for all variants together 

# jitterplot 
# facet by mutation type 

# will be using the same colors because otherwise one starts getting a rainbow tbf 
red0 = '#FC2A00'
blue0 = '#0077bb'
colors = [red0, blue0]

plt.figure(figsize=(16,4)) # set figure size

sns.stripplot(y=f'score', x='CH variant2', hue = 'CH_variant_carrier', data=df_compare_to_ref, dodge = True, jitter = True, palette = colors, size = 4, alpha = 0.3, legend = True, order = order)

plt.xlabel(f'CH hotspot variant', fontsize = 14)
plt.ylabel(r'$\log_{10}$(min %EL rank)', fontsize = 14)
plt.title('MHC II', fontsize = 16)
plt.xticks(fontsize = 12, rotation = 90)
plt.yticks(fontsize = 13)
plt.ylim(-2.5, 2.5)

# add annotation to indicate median

for i, category in enumerate(order):
            
    median_carrier = df_compare_to_ref[(df_compare_to_ref['CH variant2'] == f'{category}') & (df_compare_to_ref['CH_variant_carrier'] == 'Carrier')].score.median()
    median_noncarrier = df_compare_to_ref[(df_compare_to_ref['CH variant2'] == f'{category}') & (df_compare_to_ref['CH_variant_carrier'] == 'Non-carrier')].score.median()

    # Plot text for each hue group
    plt.text(i, median_carrier, '—', ha='right', va='center', fontsize=16, fontweight='bold', color = '#FC2A00')
    plt.text(i, median_noncarrier, '—', ha='left', va='center', fontsize=16, fontweight='bold', color = '#0077bb')

# add Mann Whitney wallis test between groups 
for i, category in enumerate(order):
    
    category_data = df_compare_to_ref[df_compare_to_ref['CH variant2'] == f'{category}']
    max_value = category_data['score'].max()

    scores_carrier = df_compare_to_ref[(df_compare_to_ref['CH variant2'] == f'{category}') & (df_compare_to_ref['CH_variant_carrier'] == 'Carrier')].score.tolist()
    scores_noncarrier = df_compare_to_ref[(df_compare_to_ref['CH variant2'] == f'{category}') & (df_compare_to_ref['CH_variant_carrier'] == 'Non-carrier')].score.tolist()

    statistic, p_value = mannwhitneyu(scores_carrier, scores_noncarrier)
    significance = ''
    if p_value > 0.05 / len(order): # adjust to the number of tests performed 
        significance = 'ns'
    elif p_value < 0.01 / len(order):
        significance = '**'
    else:
        significance = '*'
    plt.text(i, 0.2+max_value, significance, ha='center', va='bottom', fontsize=11)

plt.legend(title = 'CH status', markerscale = 2, loc = 'upper right', fontsize = 13, title_fontsize = 14)
legend = plt.gca().get_legend()
for lh in legend.legendHandles:
    lh.set_alpha(1)

plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/project_report/folders_for_git/figure 4/2b_{param}_for_all_compare_to_ref_set_2k_jitter_mhc2.pdf', bbox_inches='tight')

