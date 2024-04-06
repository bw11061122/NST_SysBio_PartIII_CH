#!/bin/bash

# 2023-11-14 
# This is a test script to automate analysis of peptide affinity using PRIME (see https://github.com/GfellerLab/PRIME)
# author: Barbara Walkowiak bw450

# specify the fasta file and HLA allele background 
# fasta file = (contains sequences of peptides to examine - kmers generated with get_kmers_from_csv.py or get_kmers_from_arg.py)

# how to run this script 
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA
# example usage: bash scripts/PRIME_specify_hla_kmers.sh data/hla_15_most_common.csv kmers/kmers_20231116
# example usage: bash scripts/PRIME_specify_hla_kmers.sh data/HLA_list_1.txt kmers/kmers_20231116

mkdir -p ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/PRIME_out/csv
mkdir -p ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/PRIME_out/txt
mkdir -p ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/log

# the kmers directory has sequences of peptides (8- to 11-mers) which can be analysed for presentation 
directory=$2

for fasta in "$directory"/*
do
    base=`basename $fasta .txt`
    hla=$1
    base_hla=`basename $hla .txt`
    hla_list=$(<$hla)

    log=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/log/${base}_${base_hla}_PRIME_log.txt
    out_txt=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/PRIME_out/txt/${base}_${base_hla}_PRIME_out.txt
    out_txt2=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/PRIME_out/txt/${base}_${base_hla}_PRIME_out_data.txt
    out_csv=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/PRIME_out/csv/${base}_${base_hla}_PRIME_out.csv
    # exec >$log 2>&1

    echo "Sample name is $base, HLA alleles examined are $hla_list"
    
    PRIME -i $fasta -o $out_txt -a $hla_list 

    # I can't be bothered to convert each of these .txt files to .csv
    grep -v '^#' $out_txt > $out_txt2
    tr -s " " < $out_txt2 | sed 's/\t/,/g' > $out_csv

    # # remove the out_data files (no need for these and can always re-create)
    # rm $out_txt2

done
