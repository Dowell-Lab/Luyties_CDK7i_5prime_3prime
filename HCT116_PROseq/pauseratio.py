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

# Define refseq input, output basename

basedir = sys.argv[1]
filelist = sys.argv[2]
pr_file = basedir + "/results/pause_region_counts.txt"
rog_file = basedir + "/results/restofgene_counts.txt"
input_file = basedir + "/" + filelist

# Dump input files into variables
cols = [
    "gene",
    "transcript",
    "length",
]

with open(input_file) as f:
    for line in f:
        files = tuple(line.strip().split("\t"))
        cols.append(files[1])

full_counts = []
with open(pr_file) as f:
    count_input = []
    i=0
    for line in f:
        i += 1
        if i == 1:
            continue
        counts = tuple(line.strip().split("\t"))
        count_input.append(counts)
        counts_to_add = [float(x) for x in counts[3:]]
        full_counts.append(counts_to_add)

pr_table = pd.DataFrame(count_input, columns=cols)

with open(rog_file) as f:
    count_input = []
    i=0
    for line in f:
        i += 1
        if i == 1:
            continue
        counts = tuple(line.strip().split("\t"))
        count_input.append(counts)
        pr = full_counts[(i-2)]
        rog = [float(x) for x in counts[3:]]
        full_counts[(i-2)] = [sum(x) for x in zip(pr, rog)]

rog_table = pd.DataFrame(count_input, columns=cols)
full_table = pd.DataFrame(full_counts, columns=cols[3:])
full_table['rsums'] = full_table.sum(axis=1)
full_table['gene'] = rog_table['gene']
full_table['transcript'] = rog_table['transcript']
full_table['length'] = rog_table['length'].astype(float)
full_table['normcounts'] = full_table['rsums']/(full_table['length'] + 330)

# Filter to and log max transcript per gene
tx_filt = full_table.sort_values(by=['gene','rsums'], ascending=[True,False]).drop_duplicates(['gene'])

# Find top 6000 transcript sums and associated ids
top_genes = tx_filt.sort_values(by='normcounts', ascending=False).iloc[0:6000,:]
tx_filt = pr_table.sort_values(by=['transcript']).drop_duplicates(['transcript'])
counts_bin1 = tx_filt[tx_filt['transcript'].isin(top_genes['transcript'])]
tx_filt = rog_table.sort_values(by=['transcript']).drop_duplicates(['transcript'])
counts_bin2 = tx_filt[tx_filt['transcript'].isin(top_genes['transcript'])]

# Iterate over bins and find averages
norm_counts1 = counts_bin1.iloc[:,3:len(counts_bin1.columns)].astype(float).div(counts_bin1['length'].astype(float),axis=0)
norm_counts2 = counts_bin2.iloc[:,3:len(counts_bin2.columns)].astype(float).div(counts_bin2['length'].astype(float),axis=0)

pauseratio_all = norm_counts1/norm_counts2

mean_counts1 = norm_counts1.mean(axis=0)
mean_counts2 = norm_counts2.mean(axis=0)
pauseratio = mean_counts1/mean_counts2

outdata = pd.DataFrame(columns=cols[3:len(cols)])
outdata.loc[0] = list(pauseratio)

# Export data

indir = pr_file.split("/")
outdir = "/".join(indir[0:-1])
outfile = "".join([outdir, "/pause_ratio_values.txt"])

outdata.to_csv(
    outfile,
    sep='\t',
    columns=list(outdata.columns),
    header=True,
    index=False,
    quoting=csv.QUOTE_NONE,
    escapechar='\\'
)

outfile = "".join([outdir, "/pause_ratios_bygene.txt"])

pauseratio_all.to_csv(
    outfile,
    sep='\t',
    columns=list(pauseratio_all.columns),
    header=True,
    index=False,
    quoting=csv.QUOTE_NONE,
    escapechar='\\'
)
