# metagene_pipeline_geneends_v2.py
# for parsing base-wise gene counts to binned metagene values

### Import packages
import csv
import numpy as np
import os.path
import pandas as pd
import sys
from pathlib import Path

### Import system arguments
# (input dir, a window size, and {3prime,5prime}))
indir = sys.argv[1]
size_factor_file = sys.argv[2]
window = int(sys.argv[3])
endtype = sys.argv[4]

### Initialize global variables
count_store,rep_store,reptimes,size_factors = {},{},{},{}
sample_names,sample_basenames = [],[]
senselist = list(Path(indir).glob('**/*.counts.bed'))
nrep = 2
conds = ['IFN','SY5609']
if endtype == '5prime':
  anticount_store,antirep_store = {},{}
  antisenselist = list(Path(indir).glob('**/*antisense_counts.bed'))
  binsizes = [20,50,100,150]
else:
  binsizes = [500,800,1000]

### Define functions
def readin_counts(sampdict,repdict,iters,filename,sample,sample_base):
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
      if ((line[7] == region['name']) & (int(line[5]) != region['start_coord'])):
        continue
      if line[7] != region['name']:
        if not firstline:
          sampdict,repdict = store_region(sampdict,repdict,counts,region,sample,sample_base,iters)
        firstline = False
        counts,region = set_region_info(region,line,firstline)
      region['coord'] = int(line[1])
      curr_coord = region['coord'] - region['start_coord']
      counts[curr_coord] = int(line[3])
  sampdict,repdict = store_region(sampdict,repdict,counts,region,sample,sample_base,iters)
  return sampdict,repdict

def set_region_info(reg,ln,first):
  cts = [0 for i in range(0,window)]
  if first:
    reg['name'] = 'A'
    reg['coord'],reg['start_coord'] = 0,0
    reg['strand'] = '+'
  else:
    reg['name'] = ln[7]
    reg['start_coord'] = int(ln[5])
    reg['strand'] = ln[9]
  return cts,reg

def store_region(sampcts,repcts,cts,reg,samp,samp_base,niter):
  if reg['strand'] == '-':
    cts.reverse()
  cts = list(np.array(cts)/size_factors[samp])
  cts = [int(i) if int(i) == i else i for i in cts]
  sampcts = add_counts(sampcts,cts,reg['name'],samp)
  repcts = add_replicate_counts(repcts,cts,reg['name'],samp_base,niter)
  return sampcts,repcts

def add_counts(store,ct,name,smp):
  if np.sum(ct) > 50: # Only include regions that have some counts
    store[smp][name] = ct
  return store

def add_replicate_counts(store,ct,name,smp,n):
  if name in store[smp].keys():
    store[smp][name] = [
      (store[smp][name][i] + ct[i]) for i in range(0,len(ct))
    ]
    if ((n[smp] > (nrep-1)) & (np.sum(store[smp][name]) < 51)): # Only include regions that have some counts
      store[smp].pop(name)
  else:
    store[smp][name] = ct
  return store

def find_regions_to_keep(sampcts,samp):
  curr_cts = pd.DataFrame(sampcts[samp])
  # Find total reads per region
  colsums = curr_cts.sum(axis=0)
  colsums = colsums.sort_values(ascending = False)
  shavepct = int(len(colsums)/20)
  # Find middle 90% of regions
  keep_regs = colsums[shavepct:(len(colsums) - shavepct + 1)]
  keep_regs = list(keep_regs.index)
  return keep_regs

def bin_and_calculate(ctdf,bsize):
  bct = []
  for j in range(0,(window-bsize-1)):
    bsum = ctdf.iloc[j:(j+bsize),:].sum(axis=0)
    bct.append(list(bsum/bsize))
  bdf = pd.DataFrame(bct,columns=bsum.index)
  # Calculate values
  metacalcs = {
    'allregs': bdf,
    'med': bdf.median(axis=1),
    'mean': bdf.mean(axis=1),
  }
  return metacalcs

def sample_calc_output(countdict,samp,regs):
  rep_df = pd.DataFrame(countdict[samp])
  rep_df = rep_df.loc[:,regs]
  # Calculate metagene values
  calcs = bin_and_calculate(rep_df,binsize)
  output_sample_data(calcs,samp)
  return calcs['mean']

def output_sample_data(metacalcs,smp):
  # Make output filenames
  outfiles = {
    'allregs': os.path.join(indir,smp + suffixes['allregs']),
    'med': os.path.join(indir,smp + suffixes['median']),
    'mean': os.path.join(indir,smp + suffixes['mean']),
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

### Execute script

# 0) Read in size factors file
with open(size_factor_file) as f:
  for readline in f:
    line = tuple(readline.strip().split("\t"))
    size_factors[line[0]] = float(line[1])

# 1) Grab count numbers by looping through intersect files
for path in senselist:
  sample_name = str(path).split('.')[0].split('/')[-1]
  sample_basename = sample_name.split('_')
  sample_basename = '_'.join(sample_basename[0:(len(sample_basename)-1)])
  if sample_basename not in conds:
    continue
  sample_names.append(sample_name)
  if sample_basename not in sample_basenames:
    sample_basenames.append(sample_basename)
    reptimes[sample_basename] = 1
  else:
    reptimes[sample_basename] += 1
  count_store,rep_store = readin_counts(count_store,rep_store,reptimes,path,sample_name,sample_basename)

for cond in conds:
  reptimes[cond] = 0

if endtype == '5prime':
  for path in antisenselist:
    sample_name = str(path).split('.')[0].split('/')[-1]
    sample_basename = sample_name.split('_')
    sample_basename = '_'.join(sample_basename[0:(len(sample_basename)-1)])
    if sample_basename not in conds:
      continue
    reptimes[sample_basename] += 1
    anticount_store,antirep_store = readin_counts(anticount_store,antirep_store,reptimes,path,sample_name,sample_basename)
    for key in count_store[sample_name].keys():
      if key not in anticount_store[sample_name].keys():
        anticount_store[sample_name][key] = [0 for i in range(0,window)]
    for key in rep_store[sample_basename].keys():
      if key not in antirep_store[sample_basename].keys():
        antirep_store[sample_basename][key] = [0 for i in range(0,window)]

# 2) Find middle 90% of genes in each sample, then find intersection
all_regs = []
all_rep_regs = []
for sample_name in sample_names:
  sample_regs = find_regions_to_keep(count_store,sample_name)
  all_regs.append(sample_regs)

for sample_basename in sample_basenames:
  sample_rep_regs = find_regions_to_keep(rep_store,sample_basename)
  all_rep_regs.append(sample_rep_regs)

final_reg_list = list(set.intersection(*map(set,all_regs)))
final_rep_reg_list = list(set.intersection(*map(set,all_rep_regs)))

# For each bin size:

for binsize in binsizes:
  suffixes = {
    'allregs': ('.' + str(binsize) + 'bpbin.' + endtype + '_metagene.allregions.txt'),
    'median': ('.' + str(binsize) + 'bpbin.' + endtype + '_metagene.median.txt'),
    'mean': ('.' + str(binsize) + 'bpbin.' + endtype + '_metagene.mean.txt'),
  }
  # 3) For each sample, extract genes, bin values, output counts and calcs
  all_sample_means = pd.DataFrame()
  all_rep_means = pd.DataFrame()
  for sample_name in sample_names:
    sample_means = sample_calc_output(
      count_store,
      sample_name,
      final_reg_list,
    )
    all_sample_means[sample_name] = sample_means
  colnames = list(all_sample_means.columns)
  colnames.sort()
  all_sample_means = all_sample_means[colnames]
  for sample_basename in sample_basenames:
    rep_means = sample_calc_output(
      rep_store,
      sample_basename,
      final_rep_reg_list,
    )
    all_rep_means[sample_basename] = rep_means
  colnames = list(all_rep_means.columns)
  colnames.sort()
  all_rep_means = all_rep_means[colnames]
  # 4) Output all sample data
  outfile = os.path.join(indir,'all_samples.' + str(binsize) + 'bpbin.' + endtype + '_metagene.mean.txt')
  all_sample_means.to_csv(
    outfile,
    sep='\t',
    header=True,
    index=True,
    quoting=csv.QUOTE_NONE,
  )
  outfile = os.path.join(indir,'all_samples_combined_replicates.' + str(binsize) + 'bpbin.' + endtype + '_metagene.mean.txt')
  all_rep_means.to_csv(
    outfile,
    sep='\t',
    header=True,
    index=True,
    quoting=csv.QUOTE_NONE,
  )

if endtype == '5prime':
  for binsize in binsizes:
    suffixes = {
      'allregs': ('.' + str(binsize) + 'bpbin.antisense.' + endtype + '_metagene.allregions.txt'),
      'median': ('.' + str(binsize) + 'bpbin.antisense.' + endtype + '_metagene.median.txt'),
      'mean': ('.' + str(binsize) + 'bpbin.antisense.' + endtype + '_metagene.mean.txt'),
    }
    # 3) For each sample, extract genes, bin values, output counts and calcs
    all_sample_means = pd.DataFrame()
    all_rep_means = pd.DataFrame()
    for sample_name in sample_names:
      sample_means = sample_calc_output(
        anticount_store,
        sample_name,
        final_reg_list,
      )
      all_sample_means[sample_name] = sample_means
    colnames = list(all_sample_means.columns)
    colnames.sort()
    all_sample_means = all_sample_means[colnames]
    for sample_basename in sample_basenames:
      rep_means = sample_calc_output(
        antirep_store,
        sample_basename,
        final_rep_reg_list,
      )
      all_rep_means[sample_basename] = rep_means
    colnames = list(all_rep_means.columns)
    colnames.sort()
    all_rep_means = all_rep_means[colnames]
    # 4) Output all sample data
    outfile = os.path.join(indir,'all_samples.' + str(binsize) + 'bpbin.antisense.' + endtype + '_metagene.mean.txt')
    all_sample_means.to_csv(
      outfile,
      sep='\t',
      header=True,
      index=True,
      quoting=csv.QUOTE_NONE,
    )
    outfile = os.path.join(indir,'all_samples_combined_replicates.' + str(binsize) + 'bpbin.antisense.' + endtype + '_metagene.mean.txt')
    all_rep_means.to_csv(
      outfile,
      sep='\t',
      header=True,
      index=True,
      quoting=csv.QUOTE_NONE,
    )

# metagene_pipeline_geneends_v2.py ends here
