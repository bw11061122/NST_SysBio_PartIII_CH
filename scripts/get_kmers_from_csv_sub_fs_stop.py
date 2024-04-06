#!/usr/bin/python3

# 20230-11-02 / 2023-11-03
# script to automate generation of k-mers from fasta sequences 

# how to run this:
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA
# python3 scripts/get_kmers_from_csv_sub_fs_stop.py data/ch_var_test.csv
# test_ch_var.csv is a dataframe that contains three columns
# the first column = gene name, the second column = variant (e.g., R882H), the third column = URL to UniProt fasta seq
# example URL to UniProt fasta seq: https://rest.uniprot.org/uniprotkb/P22681.fasta where P22681 is a UniProt ID for the proteoform of interest

# download the sequence of a given protein from NCBI
from urllib.request import urlopen
import sys
import pandas as pd
import re 

file = sys.argv[1]
data = pd.read_csv(file, sep=",", chunksize=1) # read the file
# print(data.head())

import time
timestr = time.strftime("%Y%m%d") # get current date 

# I will use this function to pull out numbers from the CH variants later on
def extract_numbers(input_string):
    # Define a regular expression pattern to match numbers in the specified patterns
    pattern = re.compile(r'\d+')
    # Use findall to extract all occurrences of the pattern in the input string
    numbers = pattern.findall(input_string)
    numbers_string = ''.join(numbers)
    return numbers_string

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
    print("The sequence of the", gene, "protein is", fasta,". The sequence is", len(fasta), "amino acids long.") 

    # identify the specific residue
    # mut = sys.argv[2] # specific variant provided as second argument (as a string, " ")
    mut = line.iloc[:,1].values[0] # you need to provide the ID of the mutation
    nr_sub = extract_numbers(mut) # determine the nr of the residue the mutation occurs in
    nr_sub = int(nr_sub) - 1 # python starts indexing at 0, so your n turns into n-1 
    original = mut[0]

    if mut[-2:] == "fs":
        print("there is a frameshift mutation")
        ch_variant = "fs"
    elif mut[-1:] == "*":
        print("there is a premature STOP codon")
        ch_variant = "*"
    else:
        ch_variant = mut[-1:] # what is the sequence at that site in the CH variant?
        print("At position:", nr_sub+1, "the reference residue is", original,". In CH, this residue is mutated to", ch_variant)

    # is R at 882 position?
    if (fasta[nr_sub] == original) == True:
        print("All good, can proceed!")
    else: 
        print("Error! There is an issue in how the fasta sequence was processed")
        exit() # stop executing the script if the residue does not match what should be in the variant 

    # find aa +/- 10 from that residue
    fasta_sub = fasta[nr_sub-10 : nr_sub+11]
    print("The sequence containing the original variant and the residues 10 aa away from it is:", fasta_sub)

    # %%
    # set path to current directory 
    import os
    path = "/Users/barbarawalkowiak/Desktop/msc_thesis"
    os.chdir(path)

    # %% 
    
    from myfunctions import build_kmers 

    def build_kmers(fasta, ksize):
        kmers = [] # initialise list 
        for i in range(ksize): 
            kmer = fasta[11-ksize+i : 11+i] 
            kmers.append(kmer)

        return kmers

    # create new folder to save k-mers into if it does not exist yet
    if not os.path.exists("/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/kmers/kmers_" + timestr): 
        os.mkdir("/Users/barbarawalkowiak/Desktop/msc_thesis/task1_predict_binding_to_HLA/kmers/kmers_" + timestr)
    
    kmers = []
    for i in range(8,12):
        ki = build_kmers(fasta_sub, i)
        kmers.append(ki)
    
    unlisted_kmers = [item for sublist in kmers for item in sublist]
    unlisted_kmers = set(unlisted_kmers)
    print(unlisted_kmers, len(unlisted_kmers))
    with open('task1_predict_binding_to_HLA/kmers/kmers_' + timestr + '/' + timestr + '_kmers_' + gene + '_' + mut + '_refseq.txt', 'w') as outfile: # refseq indicates I am using the "wild-type" sequence of the protein 
        outfile.write('\n'.join(unlisted_kmers))
    
    print("Completed writing reference k-mer sequences from", gene, "to a file")

    # %% 
    # modify the fasta_sub seqeunce so that it carries the MUTANT CH variant 
    fasta_ch = list(fasta_sub)

    # sanity check 
    if (fasta_sub[10] == original) == True:
        print("All good, can proceed!")
        
    else: 
        print("Error! There is an issue in subsetting the fasta sequence")
        exit() # stop executing the script if the residue does not match what should be in the variant 

    # the first option is that there is a frameshift mutation 
    if ch_variant == "fs":

        if gene == "NPM1" and mut == "W288fs": # this is the only fs variant I have
            print("Using sequence generated by the most common insertion TCTG") 
            npm1_path = "task1_predict_binding_to_HLA/data/npm1_w288fs_aa_sequence.txt"
            with open(npm1_path, "r") as file:
                npm1 = file.read().strip()
            print(npm1)
            npm1_sub = npm1[nr_sub-10 : nr_sub+11]
            print(npm1_sub)

            kmers = []
            
            for i in range(8,12):
                ki = build_kmers(npm1_sub, i)
                kmers.append(ki)
            
            unlisted_kmers = [item for sublist in kmers for item in sublist]
            unlisted_kmers = set(unlisted_kmers)
            print(unlisted_kmers, len(unlisted_kmers))
            with open('task1_predict_binding_to_HLA/kmers/kmers_' + timestr + '/' + timestr + '_kmers_' + gene + '_' + mut +'_ch_variant.txt', 'w') as outfile: # refseq indicates I am using the "wild-type" sequence of the protein 
                outfile.write('\n'.join(unlisted_kmers))

        else: 
            print("I don't know this mutation, you need to look it up")

    # the second option is that there is a stop codon introduced - but this does not generate any new neo-antigens 
    elif ch_variant == "*":
        print("no new neoantigens will be generated: premature stop codon")
        fasta_short = fasta_sub[:10] # the sequence ends at this residue 
        
        kmers = []
        for i in range(8,12):
            ki = build_kmers(fasta_short, i)
            kmers.append(ki)
        
        unlisted_kmers = [item for sublist in kmers for item in sublist]
        unlisted_kmers = set(unlisted_kmers)
        print(unlisted_kmers, len(unlisted_kmers))
        with open('task1_predict_binding_to_HLA/kmers/kmers_' + timestr + '/' + timestr + '_kmers_' + gene + '_' + mut + '_refseq.txt', 'w') as outfile: # refseq indicates I am using the "wild-type" sequence of the protein 
            outfile.write('\n'.join(unlisted_kmers))

    # the last option is getting a substitution we know how to deal with 
    else:       
        fasta_ch[10] = str(ch_variant) # replace with variant identified in CH
        fasta_ch = ''.join(fasta_ch)
        print("The sequence containing the ch variant and the residues 10 aa away from it is:", fasta_ch)
        kmers = []
        
        for i in range(8,12):
            ki = build_kmers(fasta_ch, i)
            kmers.append(ki)
        
        unlisted_kmers = [item for sublist in kmers for item in sublist]
        unlisted_kmers = set(unlisted_kmers)
        print(unlisted_kmers, len(unlisted_kmers))
        with open('task1_predict_binding_to_HLA/kmers/kmers_' + timestr + '/' + timestr + '_kmers_' + gene + '_' + mut +'_ch_variant.txt', 'w') as outfile: # refseq indicates I am using the "wild-type" sequence of the protein 
            outfile.write('\n'.join(unlisted_kmers))
        
    print("Completed writing ch-variant", mut, "k-mer sequences from", gene, "to a file")
