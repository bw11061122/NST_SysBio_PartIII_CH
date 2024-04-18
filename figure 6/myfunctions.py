# user-defined functions
# Author: Barbara Walkowiak (bw450) + a lot of random Stack Overflow 

import re 
import pandas as pd 
import numpy as np

# for a given fasta seq, create kmers of kize 
# takes 21-aa long seq and looks at the central position where the CH hotspot residue is
def build_kmers(fasta, ksize):
    kmers = [] # initialise list 
    for i in range(ksize): 
        kmer = fasta[11-ksize+i : 11+i] 
        kmers.append(kmer)

    return kmers

# this function takes HLA alleles written in the format of HLA-A24:51 and transforms them into format A2451
# input: HLA-A01:01; output: A_101
# input: HLA-B12:04; output: B_1204

def transform_format(input_string):
    # Define a regular expression pattern to match the input format
    pattern = re.compile(r'HLA-(\w)(\d{2}):(\d+)')

    # check if there is a match
    match = pattern.match(input_string)

    # if match, apply transformation
    if match:
        group1 = match.group(1) # leave as it is 
        group2 = int(match.group(2)) # remove zeros at the start 
        group3 = match.group(3) # leave as it is 

        # Format the output string
        output_string = f'{group1}_{group2}{group3}' # stitch back 

        return output_string # return transformed string 

    # if no much, return original string 
    return input_string

def transform_format_from_netmhc(input_string):
    # Define a regular expression pattern to match the input format
    # The input is in a format of HLA-A*01:03, HLA-C*15:27
    pattern = re.compile(r'HLA-([ABC])\*(\d{2}):(\d+)')

    # check if there is a match
    match = pattern.match(input_string)
    
    # if match, apply transformation
    if match:
        group1 = match.group(1) # leave as it is 
        group2 = int(match.group(2)) # remove zeros at the start 
        group3 = match.group(3) # leave as it is 

        # Format the output string
        output_string = f'{group1}_{group2}{group3}' # stitch back 

        return output_string # return transformed string 

    # if no much, return original string 
    return input_string

def transform_format_for_netmhc(input_string):
    # Define a regular expression pattern to match the input format
    pattern = re.compile(r'(\w_)(\d{2}):(\d+)') # expression I am looking for 
    match = pattern.match(input_string) # check if there is a match

   # if match, apply transformation
    if match:
        group1 = match.group(1) # leave as it is 
        group2 = match.group(2) # remove zeros at the start 
        group3 = match.group(3) # leave as it is 

        # Format the output string
        output_string = f'HLA-{group1}{group2}:{group3}' # stitch back 

        return output_string # return transformed string 

    # if no much, return original string 
    return input_string

def find_min_value_for_present_var(row, df):
    row_values = pd.to_numeric(row[1:-2], errors='coerce')
    
    # find alleles which are present 
    hla = row.index[1:-2][row_values >= 1]  
    
    # find values corresponding to these alleles 
    vals = df.loc[df['gene_var'] == row['gene_var'], hla].values.flatten() 
    # print(vals)
    # print(type(vals))
    # it can be that nothing was found e.g., bc there were no predictions made for this variant
    if vals.size == 0: # if found some values, take the minimum
        value = None
    else: # if no values, get NA 
        value = min(vals)
    return value

def find_max_value_for_present_var(row, df):
    row_values = pd.to_numeric(row[1:-2], errors='coerce')
    
    # find alleles which are present 
    hla = row.index[1:-2][row_values >= 1]  
    
    # find values corresponding to these alleles 
    vals = df.loc[df['gene_var'] == row['gene_var'], hla].values.flatten() 
    # print(vals)
    # print(type(vals))
    # it can be that nothing was found e.g., bc there were no predictions made for this variant
    if vals.size == 0: # if found some values, take the minimum
        value = None
    else: # if no values, get NA 
        value = max(vals)
    return value


def find_best_score_for_present_var(row, df, param):
    row_values = pd.to_numeric(row[1:-2], errors='coerce')
    
    # find alleles which are present 
    hla = row.index[1:-2][row_values >= 1]  
    
    # find values corresponding to these alleles 
    vals = df.loc[df['gene_var'] == row['gene_var'], hla].values.flatten() 
    # print(vals)
    # print(type(vals))
    # it can be that nothing was found e.g., bc there were no predictions made for this variant
    if vals.size == 0: 
        value = None
    
    else: 
        if param == "min_rank":
            value = min(vals)
        elif param == "sum_peptides_below_05":
            value = max(vals)
        elif param == "sum_peptides_below_2":
            value = max(vals)
        elif param == "min_rank_wt_minus_mut":
            value = max(vals)
        elif param == "min_rank_ratio_mut_wt":
            value = min(vals)
        elif param == "min_rank_log":
            value = max(vals)
        elif param == "min_rank_ratio_mut_wt_log":
            value = max(vals)
    return value   

def find_min_values(row, df):
    hlas = row.index[1:-3][row[1:-3] >= 1] # select alleles which each Person (row) carries
    variants = df['gene_var'] # get out variants which are present 
    min_values = {} # initialise empy dictionaries

    for var in variants:
        # Find the minimum value for each variant in the category that is present
        min_value = min(df.loc[df['gene_var'] == var, hlas].values[0])
        # Update the dictionary with the minimum value for the corresponding variant
        min_values[f'min_rank_{var}'] = min_value

    return pd.Series(min_values)

def find_best_score(row, df, param):
    hlas = row.index[1:-1][row[1:-1] >= 1] # select alleles which each Person (row) carries
    variants = df['gene_var'] # get out variants which are present 
    scores = {} # initialise empy dictionaries
    # note: param is how you can evaluate binding: min_rank, sum_peptides_below_05, sum_peptides_below_2 

    if param == "min_rank":
        for var in variants:
            # Find the minimum value for each variant in the category that is present
            min_value = min(df.loc[df['gene_var'] == var, hlas].values[0])
            # Update the dictionary with the minimum value for the corresponding variant
            scores[f'score_{var}'] = min_value
        return pd.Series(scores)

    elif param == "sum_peptides_below_2":
        for var in variants:
            # Find the minimum value for each variant in the category that is present
            min_value = max(df.loc[df['gene_var'] == var, hlas].values[0])
            # Update the dictionary with the minimum value for the corresponding variant
            scores[f'score_{var}'] = min_value
        return pd.Series(scores)

    elif param == "sum_peptides_below_05":
        for var in variants:
            # Find the minimum value for each variant in the category that is present
            min_value = max(df.loc[df['gene_var'] == var, hlas].values[0])
            # Update the dictionary with the minimum value for the corresponding variant
            scores[f'score_{var}'] = min_value
        return pd.Series(scores)

    elif param == "min_rank_wt_minus_mut":
        for var in variants:
            # Find the minimum value for each variant in the category that is present
            min_value = max(df.loc[df['gene_var'] == var, hlas].values[0])
            # Update the dictionary with the minimum value for the corresponding variant
            scores[f'score_{var}'] = min_value
        return pd.Series(scores)

    elif param == "min_rank_ratio_mut_wt":
        for var in variants:
            # Find the minimum value for each variant in the category that is present
            min_value = max(df.loc[df['gene_var'] == var, hlas].values[0])
            # Update the dictionary with the minimum value for the corresponding variant
            scores[f'score_{var}'] = min_value
        return pd.Series(scores)

    elif param == "min_rank_log":
        for var in variants:
            # we took the -10log(%rank) therefore we are looking for max value
            # Find the minimum value for each variant in the category that is present
            min_value = max(df.loc[df['gene_var'] == var, hlas].values[0])
            # Update the dictionary with the minimum value for the corresponding variant
            scores[f'score_{var}'] = min_value
        return pd.Series(scores)

    elif param == "min_rank_ratio_mut_wt_log":
        for var in variants:
            # we took -log(x/y) so we want to find the max value in the end
            # Find the minimum value for each variant in the category that is present
            min_value = max(df.loc[df['gene_var'] == var, hlas].values[0])
            # Update the dictionary with the minimum value for the corresponding variant
            scores[f'score_{var}'] = min_value
        return pd.Series(scores)

def find_best_score_log(row, df, param):
    hlas = row.index[1:-1][row[1:-1] >= 1] # select alleles which each Person (row) carries
    variants = df['gene_var'] # get out variants which are present 
    scores = {} # initialise empy dictionaries
    # note: param is how you can evaluate binding: min_rank, sum_peptides_below_05, sum_peptides_below_2 

    if param == "min_rank":
        for var in variants:
            # Find the minimum value for each variant in the category that is present
            min_value = min(np.log10(df.loc[df['gene_var'] == var, hlas].values[0]))
            # Update the dictionary with the minimum value for the corresponding variant
            scores[f'score_{var}'] = min_value
        return pd.Series(scores)

    # not sure if log10 of this makes sense but it is 11 pm and IDK sorry
    elif param == "min_rank_wt_minus_mut":
        for var in variants:
            # Find the minimum value for each variant in the category that is present
            min_value = max(np.log10(df.loc[df['gene_var'] == var, hlas].values[0]))
            # Update the dictionary with the minimum value for the corresponding variant
            scores[f'score_{var}'] = min_value
        return pd.Series(scores)

    elif param == "min_rank_ratio_mut_wt":
        for var in variants:
            # Find the minimum value for each variant in the category that is present
            min_value = min(np.log10(df.loc[df['gene_var'] == var, hlas].values[0]))
            # Update the dictionary with the minimum value for the corresponding variant
            scores[f'score_{var}'] = min_value
        return pd.Series(scores)

# Two functions defined in order to do extra calculations scores dataframe
def wt_minus_mut(row, param, df):
    if row['genotype'] == 'refseq':
        return 0
    else:
        filter_condition = (df['genotype'] == 'refseq') & (df['gene'] == row['gene']) & (df['variant'] == row['variant']) & (df['allele'] == row['allele'])
        return df.loc[filter_condition, param].values[0] - row[param]

def ratio_mut_wt(row, param, df):
    if row['genotype'] == 'refseq':
        return 1
    else:
        filter_condition = (df['genotype'] == 'refseq') & (df['gene'] == row['gene']) & (df['variant'] == row['variant']) & (df['allele'] == row['allele'])
        return row[param] / (df.loc[filter_condition, param].values[0])

# Determine if a person presents a variant based on presence of a peptide derived from the CH-variant carrying protein with a score < 0.5
def presents_variant(row, threshold):
    scores = {}
    if row['score'] < threshold: # indicates strong binding
       scores[f'presents_variant_{threshold}'] = True
    else:
        scores[f'presents_variant_{threshold}'] = False
    return pd.Series(scores)

# this is to get fasta formatted sequences 
def create_combined_column(df):
    combined_values = []

    for index, row in df.iterrows():
        sample = row['Sample']
        wt_id = row['MUTATION_ID']
        mt_id = row['MUTATION_ID']
        wt_seq = row['WT.peptide']
        mt_seq = row['MT.peptide']

        combined_values.append(f'>{sample}|{wt_id}|WT / {wt_seq} / >{sample}|{mt_id}|MT / {mt_seq}')

    combined_df = pd.DataFrame({'combined_column': combined_values})
    expanded_df = combined_df['combined_column'].str.split('/', expand=True)
    expanded_df = expanded_df.apply(lambda x: x.str.strip())  # Remove leading/trailing whitespaces
    expanded_df = expanded_df.stack().reset_index(level=1, drop=True).to_frame('combined_column')

    return expanded_df
