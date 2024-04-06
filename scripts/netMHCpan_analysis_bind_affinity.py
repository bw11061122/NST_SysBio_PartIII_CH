#! /python3

# 2023-11-26
# cd Desktop/msc_thesis/task1_predict_binding_to_HLA
# input: out/txt file (from netMHCpan_specify_hla_kmers_add_affinity.sh)
# output: a csv file which has binding affinity scores for each (specified) HLA for a list of desired variants
# author: Barbara Walkowiak bw450

# %% 
# to run: 
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/txt
# python3 ../../scripts/netMHCpan_analysis_bind_affinity.py /Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/txt

# %%

import sys
import pandas as pd
import os
import time
timestr = time.strftime("%Y%m%d") # get current date 
print (timestr) 

# %%

# identify the directory to do this in 
import csv

directory = sys.argv[1] # supply the directory when you are running the script
dfs_all = pd.DataFrame()

for file in os.listdir(directory):

    filename = os.fsdecode(file)
    if filename.endswith(".txt"): 

        input_file_path = os.path.join(directory, filename)
        

        out_filename = filename.split('.')[0]
        out_filename = out_filename + "_scores.csv"
        output_csv_path = os.path.join(directory, out_filename)

        # Initialize variables to store header and data
        header = None
        data = []

        # extract gene and variant data from file name
        print(filename)

        gene = filename.split('_')[2]
        variant = filename.split('_')[3]
        genotype = filename.split('_')[4]

        # read in the input file 
        with open(input_file_path, 'r') as file:
            lines = file.readlines()

        # Process lines to extract header and data
        for line in lines:
            if line.startswith('#'):
                # ignore comment lines 
                continue
            elif '---' in line:
                # if there is a ---, the next line is going to be the header
                header = lines[lines.index(line) + 1].strip().split()
            elif line.strip():
                # if there is sth, treat as data 
                data.append(line.strip().split())

        # check these columns for a specific value 
        column_to_check1 = 1  
        column_to_check2 = 9 

        # Specify the value to exclude rows with in the specified column
        value_to_exclude1 = "PEPLIST."  # we don't want to write rows where there is "PEPLIST."
        value_to_exclude2 = "neighbor" # we don't want to write HLA data so can exclude by looking for 'neighbour

        # Remove rows with the specified value in the specified column
        data = [row for row in data if row[column_to_check1] != value_to_exclude1]
        data = [row for row in data if row[column_to_check2] != value_to_exclude2]

        # Find and keep only one instance of each duplicate row as the header
        seen_rows = set()
        unique_rows = []
        for row in data:
            row_tuple = tuple(row)
            if row_tuple not in seen_rows:
                seen_rows.add(row_tuple)
                unique_rows.append(row)

        # columns_to_retain = [0, 1, 2, 11, 12, 13, 14, 15]

        # Extract the specified columns
        # unique_rows = [[row[i] for i in columns_to_retain] for row in data]

        df = pd.DataFrame(unique_rows)
        df.columns = df.iloc[0]
        df = df[1:].reset_index(drop=True)

        df["gene"] = gene
        df["variant"] = variant 
        df["genotype"] = genotype
        df["gene_var"] = df["gene"] + "_" + df["variant"] # add gene_var column

        dfs_all = pd.concat([dfs_all, df])

    else:
        continue

# remove wrong headers 
dfs_all = dfs_all[dfs_all["Pos"] != "Pos"] # clean up headers
dfs_all.rename(columns = {'MHC':'HLA'}, inplace = True) # rename column to match desired name going forwad

dfs_all.to_csv("../scores/" + timestr + '_NetMHC_HLA_UKBB_with_affinities.csv', index=False)
