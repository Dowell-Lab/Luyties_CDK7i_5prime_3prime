#!/usr/bin/env python
# coding: utf-8

# pauseratio.py
# This script takes in featurecounts results from two bins to find ratio

# Import libraries

import csv
import math
import numpy as np
import pandas as pd
import sys


# Define variables
basedir = sys.argv[1]
filelist = sys.argv[2]
norms = sys.argv[3]
pause_window = sys.argv[4]
if len(sys.argv) > 5:
    gene_list = sys.argv[5]
else:
    gene_list = None
pr_file = basedir + "/results/pause_region_counts.txt"
rog_file = basedir + "/results/restofgene_counts.txt"
input_file = basedir + "/" + filelist
if gene_list:
    genes_path = basedir + "/" + gene_list
cols = [
    "gene",
    "transcript",
    "length",
]
full_counts = []

def readin_filelist(inputf,inlist=False):
    with open(inputf) as f:
        data = []
        for line in f:
            files = tuple(line.strip().split("\t"))
            if inlist:
                data.append(files[1])
            else:
                data.append(files[0])
    return data

cols = cols + readin_filelist(input_file,True)
factors = readin_filelist(norms,True)
factors = [float(i) for i in factors]
norm_factors = pd.Series(factors,index=cols[3:])
if gene_list:
    genes = readin_filelist(gene_list)

def readin_countfile(inputf,col,prev_counts):
    with open(inputf) as f:
        count_input = []
        accum = []
        i=0
        for line in f:
            i += 1
            if i == 1:
                continue
            counts = tuple(line.strip().split("\t"))
            count_input.append(counts)
            counts_to_add = [float(x) for x in counts[3:]]
            if prev_counts:
                pr = prev_counts[(i-2)]
                prev_counts[(i-2)] = [sum(x) for x in zip(pr, counts_to_add)]
            else:
                accum.append(counts_to_add)
        data_table = pd.DataFrame(count_input, columns=col)
    return data_table, accum

pr_table, full_counts = readin_countfile(pr_file, cols, [])
rog_table, full_counts = readin_countfile(rog_file, cols, full_counts)

full_table = pd.DataFrame(full_counts, columns=cols[3:])
full_table = full_table.div(norm_factors,axis='columns')
keep_conds = ['total_IFN','total_SY5609']
full_table = full_table[keep_conds]

full_table['rsums'] = full_table.sum(axis=1)
full_table['gene'] = rog_table['gene']
full_table['transcript'] = rog_table['transcript']
full_table['length'] = rog_table['length'].astype(float)
full_table['normcounts'] = full_table['rsums']/(full_table['length'] + int(pause_window))

# Filter to and log max transcript per gene
tx_filt = full_table.sort_values(by=['gene','normcounts'], ascending=[True,False]).drop_duplicates(['gene'])

# Find top 8000 transcript sums and associated ids
top_genes = tx_filt.sort_values(by='normcounts', ascending=False).iloc[0:6000,:]
tx_filt = pr_table.sort_values(by=['transcript']).drop_duplicates(['transcript'])
counts_bin1 = tx_filt[tx_filt['transcript'].isin(top_genes['transcript'])]
tx_filt = rog_table.sort_values(by=['transcript']).drop_duplicates(['transcript'])
counts_bin2 = tx_filt[tx_filt['transcript'].isin(top_genes['transcript'])]

# Filter to gene list
if gene_list:
  transcripts = tx_filt.loc[tx_filt['gene'].isin(genes)]['transcript']
  counts_bin1 = pr_table.loc[pr_table['transcript'].isin(transcripts)]
  counts_bin2 = rog_table.loc[rog_table['transcript'].isin(transcripts)]
#else:
#  counts_bin1 = pr_table.copy()
#  counts_bin2 = rog_table.copy()

# Iterate over bins and find averages
norm_counts1 = counts_bin1.iloc[:,3:len(counts_bin1.columns)].astype(float).div(counts_bin1['length'].astype(float),axis=0)
norm_counts1.index = counts_bin1['transcript']
norm_counts2 = counts_bin2.iloc[:,3:len(counts_bin2.columns)].astype(float).div(counts_bin2['length'].astype(float),axis=0)
norm_counts2.index = counts_bin2['transcript']

# For total polII analysis
norm_counts1 = norm_counts1.loc[norm_counts1['total_SY5609'] > 0]
norm_counts1 = norm_counts1.loc[norm_counts1['total_IFN'] > 0]
norm_counts2 = norm_counts2.loc[norm_counts1.index]
total_counts1 = norm_counts1[keep_conds]
total_counts2 = norm_counts2[keep_conds]
pauseratio_all = total_counts1/total_counts2

# Export data

indir = pr_file.split("/")
outdir = "/".join(indir[0:-1])
outfile = "".join([outdir, "/pause_ratios_bygene.txt"])

pauseratio_all.to_csv(
    outfile,
    sep='\t',
    columns=list(pauseratio_all.columns),
    header=True,
    index=True,
    quoting=csv.QUOTE_NONE,
    escapechar='\\'
)
