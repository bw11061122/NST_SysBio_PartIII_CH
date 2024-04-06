#!/bin/bash

# 2023-11-02 
# This is a test script to automate analysis of peptide affinity using netMHCpan-4.1
# author: Barbara Walkowiak bw450

# specify the fasta file and HLA allele background 
# fasta file = (contains sequences of peptides to examine - kmers generated with get_kmers_from_csv.py or get_kmers_from_arg.py)

# how to run this script 
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA
# example usage: bash scripts/netMHCpan_specify_hla.sh HLA-A01:01 
# example usage: bash scripts/netMHCpan_specify_hla.sh data/hla_list.txt 
# example usage: bash scripts/netMHCpan_specify_hla_kmers.sh data/20231108_list_of_HLA_alleles.txt kmers/kmers_231116

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
    out=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/txt/${base}_${base_hla}_netMHC_out.txt
    out_xls=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/xls/${base}_${base_hla}_netMHC_out.xls
    out_xlsx=~/Desktop/msc_thesis/task1_predict_binding_to_HLA/netMHC_out/xls/${base}_${base_hla}_netMHC_out.xlsx
    # exec >$log 2>&1
    echo "Sample name is $base, HLA alleles examined are $hla_list"
    
    ../netMHCpan-4.1/netMHCpan -p -a $hla_list $fasta > $out # run prediction 
    ../netMHCpan-4.1/netMHCpan -p $fasta -BA -xls -a $hla_list -xlsfile $out_xls # save as .xls file 

    ssconvert $out_xls $out_xlsx

done

