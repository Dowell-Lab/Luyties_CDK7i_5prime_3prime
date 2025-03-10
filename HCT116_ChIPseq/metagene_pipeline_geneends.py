# metagene_pipeline_geneends.py
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
dataset = sys.argv[5]

### Initialize global variables
count_store,size_factors = {},{}
sample_names = []

if dataset == 'totalpolII_ser7_only':
  senselist = list(Path(indir).glob('**/5609_B2_N*.counts.bed'))
  senselist = senselist + list(Path(indir).glob('**/5609_B2_7*.counts.bed'))
elif dataset == 'totalpolII_ser7_IFN_SY_only':
  senselist = [
    os.path.join(indir,'5609_B2_N1.counts.bed'),
    os.path.join(indir,'5609_B2_N2.counts.bed'),
    os.path.join(indir,'5609_B2_71.counts.bed'),
    os.path.join(indir,'5609_B2_72.counts.bed'),
  ]
else:
  senselist = list(Path(indir).glob('**/*.counts.bed'))

if endtype == '5prime':
  binsizes = [20,50,100,150]
else:
  binsizes = [500,800,1000]

### Define functions
def readin_counts(sampdict,filename,sample):
  if sample not in sampdict.keys():
    sampdict[sample] = {}
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
      if line[7] != region['name']:
        if not firstline:
          sampdict = store_region(sampdict,counts,region,sample)
        firstline = False
        counts,region = set_region_info(region,line,firstline)
      region['coord'] = int(line[1])
      curr_coord = region['coord'] - region['start_coord']
      counts[curr_coord] = int(line[3])
  sampdict = store_region(sampdict,counts,region,sample)
  return sampdict

def set_region_info(reg,ln,first):
  cts = [0 for i in range(0,window)]
  if first:
    reg['name'] = 'A'
    reg['coord'],reg['start_coord'],reg['end_coord'] = 0,0,0
    reg['strand'] = '+'
  else:
    reg['name'] = ln[7]
    reg['start_coord'] = int(ln[5])
    reg['end_coord'] = int(ln[6])
    reg['strand'] = ln[9]
  return cts,reg

def store_region(sampcts,cts,reg,samp):
  if reg['strand'] == '-':
    cts.reverse()
  cts = list(np.array(cts)/size_factors[samp])
  cts = [int(i) if int(i) == i else i for i in cts]
  sampcts = add_counts(sampcts,cts,reg['name'],samp)
  return sampcts

def add_counts(store,ct,name,smp):
  if np.sum(ct) > 50: # Only include regions that have some counts
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
  bdf = pd.DataFrame(bct,columns=ctdf.columns)
  # Calculate values
  metacalcs = {
    'allregs': bdf,
    'mean': pd.DataFrame(bdf.mean(axis=1),columns=['mean_counts']),
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
    'mean': os.path.join(indir,smp + suffixes['mean']),
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

### Execute script

# 0) Read in size factors file
with open(size_factor_file) as f:
  for readline in f:
    line = tuple(readline.strip().split("\t"))
    size_factors[line[0]] = float(line[1])

# 1) Grab count numbers by looping through intersect files
for path in senselist:
  sample_name = str(path).split('.')[0].split('/')[-1]
  sample_names.append(sample_name)
  count_store = readin_counts(count_store,path,sample_name)

# 2) Find middle 90% of genes in each sample, then find intersection
all_regs = []
all_rep_regs = []
for sample_name in sample_names:
  sample_regs = find_regions_to_keep(count_store,sample_name)
  all_regs.append(sample_regs)

final_reg_list = list(set.intersection(*map(set,all_regs)))

# For each bin size:

for binsize in binsizes:
  suffixes = {
    'allregs': ('.' + str(binsize) + 'bpbin.' + endtype + '_metagene.allregions.txt'),
    'mean': ('.' + str(binsize) + 'bpbin.' + endtype + '_metagene.mean.txt'),
  }
  # 3) For each sample, extract genes, bin values, output counts and calcs
  all_sample_means = pd.DataFrame()
  for sample_name in sample_names:
    sample_means = sample_calc_output(
      count_store,
      sample_name,
      final_reg_list,
    )
    all_sample_means[sample_name] = pd.Series(sample_means['mean_counts'])
  colnames = list(all_sample_means.columns)
  colnames.sort()
  all_sample_means = all_sample_means[colnames]

  # 4) Output all sample data
  outfile = os.path.join(indir,'all_samples.' + str(binsize) + 'bpbin.' + endtype + '_metagene.mean.txt')
  all_sample_means.to_csv(
    outfile,
    sep='\t',
    header=True,
    index=True,
    quoting=csv.QUOTE_NONE,
  )

# metagene_pipeline_geneends_v2.py ends here
