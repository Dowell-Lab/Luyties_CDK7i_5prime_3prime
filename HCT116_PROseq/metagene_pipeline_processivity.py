# metagene_pipeline_processivity.py
# for parsing counts across genes into binned values

### Import packages
import csv
import math
import numpy as np
import os.path
import pandas as pd
import sys
from pathlib import Path

### Import system arguments
# (input dir and endtype (processivity here))
indir = sys.argv[1]
size_factor_file = sys.argv[2]
endtype = sys.argv[3]

### Initialize global variables
pathlist = list(Path(indir).glob('**/*counts.bed'))
#numbins_list = [10,100]
numbins_list = [10]
size_factors = {}

### Define functions
def readin_counts(countdict,repdict,filename,sample,sample_base,numbins):
  if sample not in countdict.keys():
    countdict[sample] = {}
  if sample_base not in repdict.keys():
    repdict[sample_base] = {}
  with open(filename) as f:
    counts = [0]
    gene = 'A'
    strand = '+'
    genelength = 1
    firstline = True
    for line in f:
      countin = tuple(line.strip().split("\t"))
      if countin[7] != gene:
        if not firstline:
          if strand == '-':
            counts.reverse()
          counts = list(np.array(counts)/size_factors[sample])
          counts = [int(i) if int(i) == i else i for i in counts]
          countdict = add_counts(countdict,sample,counts,gene,genelength,numbins)
          repdict = add_replicate_counts(repdict,sample_base,counts,gene,genelength,numbins)
        firstline = False
        genelength = (int(countin[6]) - int(countin[5]))
        counts = [0 for i in range(0,genelength)]
        gene = countin[7]
        start_coord = int(countin[5])
        strand = countin[9]
      coord = int(countin[1]) - start_coord
      counts[coord] = int(countin[3])
  if strand == '-':
    counts.reverse()
  counts = list(np.array(counts)/size_factors[sample])
  counts = [int(i) if int(i) == i else i for i in counts]
  countdict = add_counts(countdict,sample,counts,gene,genelength,numbins)
  repdict = add_replicate_counts(repdict,sample_base,counts,gene,genelength,numbins)
  return countdict,repdict

def add_counts(countdict,sample,rawcounts,gene,genelength,numbins):
  if np.sum(rawcounts) > 100: # Only include genes that have some counts
    binsize = int(genelength/numbins)
    binned = []
    for i in range(0,(numbins-1)):
      binsum = sum(rawcounts[i*binsize:(i+1)*binsize])
      binned.append(binsum/binsize)
    last_binsize = len(rawcounts[(numbins-1)*binsize:])
    binsum = sum(rawcounts[(numbins-1)*binsize:])
    binned.append(binsum/last_binsize)
    countdict[sample][gene] = binned
  return countdict

def add_replicate_counts(repdict,sample_base,repcounts,gene,genelength,numbins):
  if gene in repdict[sample_base].keys():
    repdict[sample_base][gene] = [
      ((repdict[sample_base][gene][i] + repcounts[i])/2) for i in range(0,len(repcounts))
    ]
    if np.sum(repdict[sample_base][gene]) < 101: # Only include genes that have some counts
      repdict[sample_base].pop(gene)
    else:
      repcounts = repdict[sample_base][gene]
      binsize = int(genelength/numbins)
      binned = []
      for i in range(0,(numbins-1)):
        binsum = sum(repcounts[i*binsize:(i+1)*binsize])
        binned.append(binsum/binsize)
      last_binsize = len(repcounts[(numbins-1)*binsize:])
      binsum = sum(repcounts[(numbins-1)*binsize:])
      binned.append(binsum/last_binsize)
      repdict[sample_base][gene] = binned
  else:
    repdict[sample_base][gene] = repcounts
  return repdict

def find_genes_to_keep(countdict,sample):
  current_counts = pd.DataFrame(countdict[sample])
  # Find total reads per gene
  colsums = current_counts.sum(axis=0)
  colsums_sorted = colsums.sort_values(ascending = False)
  fivepct = int(len(colsums_sorted)/20)
  # Find middle 90% of genes
  keep_genes = colsums_sorted[fivepct:(len(colsums_sorted) - fivepct + 1)]
  keep_genes = list(keep_genes.index)
  return keep_genes

def sample_calc_output(countdict,sample,genes):
  current_counts = pd.DataFrame(countdict[sample])
  # Extract common genes across samples
  current_counts = current_counts[genes]
  # Calculate metagene values
  calculations = {
    'allgenes': current_counts,
    'med': current_counts.median(axis=1),
    'mean': current_counts.mean(axis=1),
  }
  if numbin == 10:
    # Do last bin over 2nd bin calculation
    binratio = current_counts.iloc[9,:]/current_counts.iloc[1,:]
    calculations['binratio'] = binratio
  output_sample_data(calculations,sample)
  return calculations

def output_sample_data(metacalcs,sample):
  # Make output filenames
  outfiles = {
    'allgenes': os.path.join(indir,sample + suffixes['allgenes']),
    'binratio': os.path.join(indir,sample + suffixes['binratio']),
    'med': os.path.join(indir,sample + suffixes['median']),
    'mean': os.path.join(indir,sample + suffixes['mean']),
  }
  # Output data
  for key in metacalcs.keys():
    metacalcs[key].to_csv(
      outfiles[key],
      sep='\t',
      header=False,
      index=True,
      quoting=csv.QUOTE_NONE,
    )

### Execute script for each binnum

# 0) Read in size factors file
with open(size_factor_file) as f:
  for line in f:
    inline = tuple(line.strip().split("\t"))
    size_factors[inline[0]] = float(inline[1])

for numbin in numbins_list:

  count_store = {}
  rep_store = {}
  sample_names = []
  sample_basenames = []
  # 1) Grab count numbers by looping through intersect files
  for path in pathlist:
    sample_name = str(path).split('.')[0].split('/')[-1]
    sample_names.append(sample_name)
    sample_basename = sample_name.split('_')
    sample_basename = '_'.join(sample_basename[0:(len(sample_basename)-1)])
    if sample_basename not in sample_basenames:
      sample_basenames.append(sample_basename)
    count_store,rep_store = readin_counts(count_store,rep_store,path,sample_name,sample_basename,numbin)
  for sample_basename in sample_basenames:
    genes_to_scan = list(rep_store[sample_basename].keys())
    for gid in genes_to_scan:
      if len(rep_store[sample_basename][gid]) != numbin:
        rep_store[sample_basename].pop(gid)

  # 2) Find middle 90% of genes in each sample, then find intersection
  all_genes = []
  all_rep_genes = []
  for sample_name in sample_names:
    sample_genes = find_genes_to_keep(count_store,sample_name)
    all_genes.append(sample_genes)
  for sample_basename in sample_basenames:
    sample_rep_genes = find_genes_to_keep(rep_store,sample_basename)
    all_rep_genes.append(sample_rep_genes)
  final_gene_list = list(set.intersection(*map(set,all_genes)))
  final_rep_gene_list = list(set.intersection(*map(set,all_rep_genes)))

  suffixes = {
    'allgenes': ('.' + str(numbin) + 'bins.' + endtype + '.allgenes.txt'),
    'binratio': ('.' + str(numbin) + 'bins.' + endtype + '.binratio.txt'),
    'median': ('.' + str(numbin) + 'bins.' + endtype + '.median.txt'),
    'mean': ('.' + str(numbin) + 'bins.' + endtype + '.mean.txt'),
  }

  # 3) For each sample, extract genes, output counts and calcs
  all_sample_means = pd.DataFrame()
  all_rep_means = pd.DataFrame()
  all_sample_binratios = pd.DataFrame()
  all_rep_binratios = pd.DataFrame()
  for sample_name in sample_names:
    sample_calcs = sample_calc_output(
      count_store,
      sample_name,
      final_gene_list,
    )
    all_sample_means[sample_name] = sample_calcs['mean']
    all_sample_binratios[sample_name] = sample_calcs['binratio']

  for sample_basename in sample_basenames:
    rep_calcs = sample_calc_output(
      rep_store,
      sample_basename,
      final_rep_gene_list,
    )
    all_rep_means[sample_basename] = rep_calcs['mean']
    all_rep_binratios[sample_basename] = rep_calcs['binratio']

  colnames = list(all_sample_means.columns)
  colnames.sort()
  all_sample_means = all_sample_means[colnames]
  all_sample_binratios = all_sample_binratios[colnames]
  colnames = list(all_rep_means.columns)
  colnames.sort()
  all_rep_means = all_rep_means[colnames]
  all_rep_binratios = all_rep_binratios[colnames]
  # 4) Output all sample data
  outfile = os.path.join(indir,'all_samples.' + str(numbin) + 'bins.' + endtype + '.mean.txt')
  all_sample_means.to_csv(
    outfile,
    sep='\t',
    header=True,
    index=True,
    quoting=csv.QUOTE_NONE,
  )
  outfile = os.path.join(indir,'all_samples.' + str(numbin) + 'bins.' + endtype + '.binratio.txt')
  all_sample_binratios.to_csv(
    outfile,
    sep='\t',
    header=True,
    index=True,
    quoting=csv.QUOTE_NONE,
  )
  outfile = os.path.join(indir,'all_samples_combined_replicates.' + str(numbin) + 'bins.' + endtype + '.mean.txt')
  all_rep_means.to_csv(
    outfile,
    sep='\t',
    header=True,
    index=True,
    quoting=csv.QUOTE_NONE,
  )
  outfile = os.path.join(indir,'all_samples_combined_replicates.' + str(numbin) + 'bins.' + endtype + '.binratio.txt')
  all_rep_binratios.to_csv(
    outfile,
    sep='\t',
    header=True,
    index=True,
    quoting=csv.QUOTE_NONE,
  )

# metagene_pipeline_processivity.py ends here
