#!/usr/bin/python3

# 2024-02-08 
# script to automate generation of k-mers from fasta sequences 

# how to run this:
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA
# python3 scripts/get_kmers_from_csv_all_updated_class2.py data/all_variants_url.csv
# test_ch_var.csv is a dataframe that contains three columns
# the first column = gene name, the second column = variant (e.g., R882H), the third column = URL to UniProt fasta seq
# example URL to UniProt fasta seq: https://rest.uniprot.org/uniprotkb/P22681.fasta where P22681 is a UniProt ID for the proteoform of interest

# MAIN UPDATE FOR MHC CLASS II ALLELES: predictions are made only for 15 amino acid-long peptides
# "for class II only one length is admitted with 15 being the default value"
# see source publication for NetMHCII_pan4.0 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7319546/ 

# download the sequence of a given protein from NCBI
from urllib.request import urlopen
import sys
import pandas as pd

file = sys.argv[1]
data = pd.read_csv(file, sep=",", chunksize=1) # read the file
# print(data.head())

import time
timestr = time.strftime("%Y%m%d") # get current date 

for line in data:

# link = str(sys.argv[1]) # use if you want to provide this as command line 
    url = line.iloc[:,2].values[0]
    print(url)
    url = urlopen(url)
    fasta = url.read().decode("utf-8", "ignore") # import the url sequence 

    # get what gene you are looking at 
    import re 
    desc=fasta.split("\n",1)[0] # get the description of the sequence 
    gene = re.sub(".+GN=", "", desc)
    gene = re.sub(r'\s* PE.*', '', gene)
    print("The name of the gene being investigated is", gene)

    # get the protein sequence 
    fasta = fasta.split("\n",1)[1] # get only the protein sequence 
    fasta = "".join(fasta.split()) # remove spaces and lines 
    print("The sequence of the", gene, "protein is", fasta, "The sequence is", len(fasta), "amino acids long.") 

    # identify the specific residue
    # mut = sys.argv[2] # specific variant provided as second argument (as a string, " ")
    mut = line.iloc[:,1].values[0] # you need to provide the ID of the mutation
    nr_sub = int(mut[1:-1]) # determine the nr of the residue the mutation occurs in
    nr_sub_ref = int(nr_sub) - 1 # python starts indexing at 0, so your n turns into n-1 
    original = mut[0]
    ch_variant = mut[-1:] # what is the sequence at that site in the CH variant?
    print("At position:", nr_sub, "the reference residue is", original, "In CH, this residue is mutated to", ch_variant)

    # is R at 882 position?
    if (fasta[nr_sub_ref] == original) == True:
        print("All good, can proceed!")
    else: 
        print("Error! There is an issue in how the fasta sequence was processed")
        exit() # stop executing the script if the residue does not match what should be in the variant 

    # find aa +/- 15 from that residue

    fasta_sub = fasta[nr_sub-min(15, nr_sub_ref) : nr_sub+min(14, (len(fasta) - nr_sub))]
    print("The sequence containing the original variant and the residues 15 aa away from it is:", fasta_sub)

    # %%
    # set path to current directory 
    import os
    path = "/Users/barbarawalkowiak/Desktop/msc_thesis"
    os.chdir(path)

    # %% 
    
    from myfunctions import build_kmers 

    def build_kmers(fasta, ksize):
        for i in range(ksize):
            kmer = fasta[15-ksize+i : 15+i] 
            kmers.append(kmer)
        return kmers

    # create new folder to save k-mers into if it does not exist yet
    if not os.path.exists("/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/kmers/kmers15_" + timestr): 
        os.mkdir("/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/kmers/kmers15_" + timestr)
    
    kmers = []
    
    ki = build_kmers(fasta_sub, 15) # only build kmers of length 15
    print(ki)

    with open('task1_predict_binding_to_HLA/kmers/kmers15_' + timestr + '/' + timestr + '_kmers_' + gene + '_' + mut + '_refseq.txt', 'w') as outfile: # refseq indicates I am using the "wild-type" sequence of the protein 
        for i, item in enumerate(ki, start=1):
            if len(item) >= 15:
                outfile.write(f">seq_{i}\n")
                outfile.write(f"{item}\n")
    
    print("Completed writing reference k-mer sequences from", gene, "to a file")

    # %% 
    # modify the fasta_sub seqeunce so that it carries the MUTANT CH variant 
    fasta_ch = list(fasta_sub)

    # sanity check 
    if (fasta_sub[min(14, nr_sub_ref)] == original) == True: # note that indexing starts from 0 so if you want aa 15 you need position 14
        print("All good, can proceed!")
    else: 
        print("Error! There is an issue in subsetting the fasta sequence")
        exit() # stop executing the script if the residue does not match what should be in the variant 

    if ch_variant  == '*':

        # find the last 
        fasta_sub_ch = fasta[nr_sub-16:nr_sub-1] # you just take the last 15 amino acids 
        
        with open('task1_predict_binding_to_HLA/kmers/kmers15_' + timestr + '/' + timestr + '_kmers_' + gene + '_' + mut +'_ch_variant.txt', 'w') as outfile: # refseq indicates I am using the "wild-type" sequence of the protein 
            if len(fasta_sub_ch) >= 15: # only write sequences longer than 15 aa
                outfile.write(f">seq_1\n")
                outfile.write(f"{fasta_sub_ch}\n")

    else:
        fasta_ch[min(14, nr_sub_ref)] = str(ch_variant) # replace with variant identified in CH
        fasta_ch = ''.join(fasta_ch)
        print("The sequence containing the ch variant and the residues 15 aa away from it is:", fasta_ch)
        kmers = []

        ki = build_kmers(fasta_ch, 15) # only build kmers of length 15

        with open('task1_predict_binding_to_HLA/kmers/kmers15_' + timestr + '/' + timestr + '_kmers_' + gene + '_' + mut +'_ch_variant.txt', 'w') as outfile: # refseq indicates I am using the "wild-type" sequence of the protein 
            for i, item in enumerate(ki, start=1):
                if len(item) >= 15: # only write sequences longer than 15 aa
                    outfile.write(f">seq_{i}\n")
                    outfile.write(f"{item}\n")
        
        print("Completed writing ch-variant", mut, "k-mer sequences from", gene, "to a file")

# %%
