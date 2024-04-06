#!/bin/bash

# 2023-11-24
# This is a test script to automate analysis of peptide affinity using netMHCpan-4.1
# author: Barbara Walkowiak bw450

# specify the fasta file and HLA allele background 
# fasta file = (contains sequences of peptides to examine - kmers generated with get_kmers_from_csv.py or get_kmers_from_arg.py)
# I modified this script because I need to get binding affinity values for the Luksza 2017 model and these are not saved in the CSV/XLS format

# how to run this script 
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA
# example usage: bash scripts/netMHCpan_specify_hla_kmers_add_affinity.sh data/20231108_list_of_HLA_alleles.txt kmers/kmers_231116

mkdir -p ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/xls
mkdir -p ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/txt
mkdir -p ~/Desktop/msc_thesis/task1_predict_binding_to_HLA/log

# read the file with HLA alleles
directory=$2

for fasta in "$directory"/*
do

    base=`basename $fasta .txt`
    hla=$1
    base_hla=`basename $hla .txt`
    hla_list=$(<$hla)

    log=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/log/${base}_${base_hla}_netMHC4.1_log.txt
    out=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/txt/${base}_${base_hla}_netMHC_out_affinities.txt
    out_xls=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/xls/${base}_${base_hla}_netMHC_out_affinities.xls
    out_xlsx=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/xls/${base}_${base_hla}_netMHC_out.xlsx
    # exec >$log 2>&1
    echo "Sample name is $base, HLA alleles examined are $hla_list"
    
    ../netMHCpan-4.1/netMHCpan -a $hla_list -BA -f $fasta > $out # run prediction 
    ../netMHCpan-4.1/netMHCpan -p $fasta -BA -xls -a $hla_list -xlsfile $out_xls # save as .xls file 
    ssconvert $out_xls $out_xlsx

done


