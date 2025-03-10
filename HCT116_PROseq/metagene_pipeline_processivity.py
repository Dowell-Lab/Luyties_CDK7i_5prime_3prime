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
count_store,rep_store,reptimes,size_factors = {},{},{},{}
sample_names,sample_basenames = [],[]
pathlist = list(Path(indir).glob('**/*counts.bed'))
#numbins_list = [10,100]
conds = ['DMSO','IFN','SY5609','KB0742','SY_KB'] # File condition basenames
numbins_list = [10]
nrep = 2

### Define functions
def readin_counts(sampdict,repdict,iters,filename,sample,sample_base,nbin):
  if sample not in sampdict.keys():
    sampdict[sample] = {}
  if sample_base not in repdict.keys():
    repdict[sample_base] = {}
  with open(filename) as f:
    region = {}
    firstline = True
    counts,region = set_region_info(region,False,firstline)
    for readline in f:
      line = tuple(readline.strip().split("\t"))
      if int(line[1]) == region['coord']:
        continue
      if (line[7] == region['name']) & (int(line[5]) != region['start_coord']):
        continue
      if (line[7] == region['name']) & (int(line[6]) != region['end_coord']):
        continue
      if (line[7] != region['name']) & (int(line[5]) in range((region['start_coord']-1),(region['end_coord']+1))):
        continue
      if (line[7] != region['name']) & (int(line[6]) in range((region['start_coord']-1),(region['end_coord']+1))):
        continue
      elif line[7] != region['name']:
        if not firstline:
          sampdict,repdict = store_region(sampdict,repdict,counts,region,sample,sample_base,iters,nbin)
        firstline = False
        counts,region = set_region_info(region,line,firstline)
      region['coord'] = int(line[1])
      curr_coord = region['coord'] - region['start_coord']
      counts[curr_coord] = int(line[3])
  sampdict,repdict = store_region(sampdict,repdict,counts,region,sample,sample_base,iters,nbin)
  return sampdict,repdict

def set_region_info(reg,ln,first):
  if first:
    reg['name'] = 'A'
    reg['coord'],reg['start_coord'],reg['end_coord'], reg['length'] = 0,0,0,0
    reg['strand'] = '+'
    cts = [0]
  else:
    reg['name'] = ln[7]
    reg['start_coord'] = int(ln[5])
    reg['end_coord'] = int(ln[6])
    reg['strand'] = ln[9]
    reg['length'] = reg['end_coord'] - reg['start_coord']
    cts = [0 for i in range(0,(reg['length']))]
  return cts,reg

def store_region(sampcts,repcts,cts,reg,samp,samp_base,niter,numbins):
  if reg['strand'] == '-':
    cts.reverse()
  cts = list(np.array(cts)/size_factors[samp])
  cts = [int(i) if int(i) == i else i for i in cts]
  sampcts = add_counts(sampcts,cts,reg['name'],samp,reg['length'],numbins)
  repcts = add_replicate_counts(repcts,cts,reg['name'],samp_base,reg['length'],niter,numbins)
  return sampcts,repcts

def add_counts(store,ct,name,smp,reglen,nbins):
  if np.sum(ct) > 100: # Only include genes that have some counts
    binsize = int(reglen/nbins)
    binned = []
    for i in range(0,(nbins-1)):
      binsum = sum(ct[i*binsize:(i+1)*binsize])
      binned.append(binsum/binsize)
    last_binsize = len(ct[(nbins-1)*binsize:])
    binsum = sum(ct[(nbins-1)*binsize:])
    binned.append(binsum/last_binsize)
    store[smp][name] = binned
  return store

def add_replicate_counts(store,ct,name,smp,reglen,n,nbins):
  if name in store[smp].keys():
    store[smp][name] = [
      (store[smp][name][i] + ct[i]) for i in range(0,reglen)
    ]
    if ((n[smp] > (nrep-1)) & (np.sum(store[smp][name]) < 101)): # Only include genes that have some counts
      store[smp].pop(name)
    else:
      ct = store[smp][name]
      binsize = int(reglen/nbins)
      binned = []
      for i in range(0,(nbins-1)):
        binsum = sum(ct[i*binsize:(i+1)*binsize])
        binned.append(binsum/binsize)
      last_binsize = len(ct[(nbins-1)*binsize:])
      binsum = sum(ct[(nbins-1)*binsize:])
      binned.append(binsum/last_binsize)
      store[smp][name] = binned
  else:
    store[smp][name] = ct
  return store

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
    'mean': pd.DataFrame(current_counts.mean(axis=1),columns=['mean_counts'])
  }
  if numbin == 10:
    # Do last bin over 2nd bin calculation
    binratio = current_counts.iloc[9,:]/current_counts.iloc[1,:]
    calculations['binratio'] = pd.DataFrame(binratio,columns=['10th_over_2nd_binratio'])
  output_sample_data(calculations,sample)
  return calculations

def output_sample_data(metacalcs,sample):
  # Make output filenames
  outfiles = {
    'allgenes': os.path.join(indir,sample + suffixes['allgenes']),
    'binratio': os.path.join(indir,sample + suffixes['binratio']),
    'mean': os.path.join(indir,sample + suffixes['mean']),
  }
  # Output data
  for key in metacalcs.keys():
    metacalcs[key].to_csv(
      outfiles[key],
      sep='\t',
      header=True,
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

  # 1) Grab count numbers by looping through intersect files
  for path in pathlist:
    sample_name = str(path).split('.')[0].split('/')[-1]
    sample_names.append(sample_name)
    sample_basename = sample_name.split('_')
    sample_basename = '_'.join(sample_basename[0:(len(sample_basename)-1)])
    if sample_basename not in conds:
      continue
    if sample_basename not in sample_basenames:
      sample_basenames.append(sample_basename)
      reptimes[sample_basename] = 1
    else:
      reptimes[sample_basename] += 1
    count_store,rep_store = readin_counts(count_store,rep_store,reptimes,path,sample_name,sample_basename,numbin)
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
    all_sample_means[sample_name] = pd.Series(sample_calcs['mean']['mean_counts'])
    all_sample_binratios[sample_name] = pd.Series(sample_calcs['binratio']['10th_over_2nd_binratio'])

  for sample_basename in sample_basenames:
    rep_calcs = sample_calc_output(
      rep_store,
      sample_basename,
      final_rep_gene_list,
    )
    all_rep_means[sample_basename] = pd.Series(rep_calcs['mean']['mean_counts'])
    all_rep_binratios[sample_basename] = pd.Series(rep_calcs['binratio']['10th_over_2nd_binratio'])

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
