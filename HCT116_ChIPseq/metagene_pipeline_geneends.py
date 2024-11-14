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

### Initialize global variables
count_store = {}
rep_store = {}
size_factors = {}
sample_names = []
sample_basenames = []
senselist = list(Path(indir).glob('**/*.counts.bed')) # For ChIP-seq
if endtype == '5prime':
  binsizes = [5,20,50,100]
else:
  binsizes = [500,800,1000]

### Define functions
def readin_counts(countdict,repdict,filename,sample,sample_base):
  if sample not in countdict.keys():
    countdict[sample] = {}
  if sample_base not in repdict.keys():
    repdict[sample_base] = {}
  with open(filename) as f:
    counts = [0 for i in range(0,window)]
    curr_gene = 'A'
    curr_coord = 0
    start_coord = 0
    strand = '+'
    firstline = True
    for line in f:
      countin = tuple(line.strip().split("\t"))
      if int(countin[1]) == curr_coord:
        continue
      if ((countin[7] == curr_gene) & (int(countin[5]) != start_coord)):
        continue
      if countin[7] != curr_gene:
        if not firstline:
          if strand == '-':
            counts.reverse()
          counts = list(np.array(counts)/size_factors[sample])
          counts = [int(i) if int(i) == i else i for i in counts]
          countdict = add_counts(countdict,sample,counts,curr_gene)
          repdict = add_replicate_counts(repdict,sample_base,counts,curr_gene)
        firstline = False
        counts = [0 for i in range(0,window)]
        curr_gene = countin[7]
        start_coord = int(countin[5])
        strand = countin[9]
      curr_coord = int(countin[1])
      coord = curr_coord - start_coord
      counts[coord] = int(countin[3])
  if strand == '-':
    counts.reverse()
  counts = list(np.array(counts)/size_factors[sample])
  counts = [int(i) if int(i) == i else i for i in counts]
  countdict = add_counts(countdict,sample,counts,curr_gene)
  repdict = add_replicate_counts(repdict,sample_base,counts,curr_gene)
  return countdict,repdict

def add_counts(countdict,sample,rawcounts,gene):
  if np.sum(rawcounts) > 50: # Only include genes that have some counts
    countdict[sample][gene] = rawcounts
  return countdict

def add_replicate_counts(repdict,sample_base,repcounts,gene):
  if gene in repdict[sample_base].keys():
    repdict[sample_base][gene] = [
      (repdict[sample_base][gene][i] + repcounts[i]) for i in range(0,len(repcounts))
    ]
    if np.sum(repcounts) < 51: # Only include genes that have some counts
      repdict[sample_base].pop(gene)
  else:
    repdict[sample_base][gene] = repcounts
  return repdict

def find_genes_to_keep(countdict,sample):
  current_counts = pd.DataFrame(countdict[sample])
  # Find total reads per gene
  colsums = current_counts.sum(axis=0)
  colsums_sorted = colsums.sort_values(ascending = False)
  shavepct = int(len(colsums_sorted)/20)
  # Find middle 95% of genes
  keep_genes = colsums_sorted[shavepct:(len(colsums_sorted) - shavepct + 1)]
  keep_genes = list(keep_genes.index)
  return keep_genes

def sample_calc_output(countdict,sample,genes):
  current_counts = pd.DataFrame(countdict[sample])
  # Extract common genes across samples
  current_counts = current_counts[genes]
  calculations = bin_and_calculate(current_counts)
  output_sample_data(calculations,sample)
  return calculations['mean']

def bin_and_calculate(countdf):
  binned = []
  for i in range(0,(window-binsize-1)):
    binsum = countdf.iloc[i:(i+binsize),:].sum(axis=0)
    binned.append(binsum/binsize)
  binneddf = pd.DataFrame(binned)
  # Calculate metagene values
  metacalcs = {
    'allgenes': binneddf,
    'med': binneddf.median(axis=1),
    'mean': binneddf.mean(axis=1),
  }
  return metacalcs

def output_sample_data(metacalcs,sample):
  # Make output filenames
  outfiles = {
    'allgenes': os.path.join(indir,sample + suffixes['allgenes']),
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

### Execute script

# 0) Read in size factors file
with open(size_factor_file) as f:
  for line in f:
    inline = tuple(line.strip().split("\t"))
    size_factors[inline[0]] = float(inline[1])

# 1) Grab count numbers by looping through intersect files
for path in senselist:
  sample_name = str(path).split('.')[0].split('/')[-1]
  sample_names.append(sample_name)
  sample_basename = sample_name.split('_')
  sample_basename = '_'.join(sample_basename[0:(len(sample_basename)-1)])
  if sample_basename not in sample_basenames:
    sample_basenames.append(sample_basename)
  count_store,rep_store = readin_counts(count_store,rep_store,path,sample_name,sample_basename)

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

# For each bin size:

for binsize in binsizes:
  suffixes = {
    'allgenes': ('.' + str(binsize) + 'bpbin.' + endtype + '_metagene.allgenes.txt'),
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
      final_gene_list,
    )
    all_sample_means[sample_name] = sample_means
  colnames = list(all_sample_means.columns)
  colnames.sort()
  all_sample_means = all_sample_means[colnames]
  for sample_basename in sample_basenames:
    rep_means = sample_calc_output(
      rep_store,
      sample_basename,
      final_rep_gene_list,
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

# metagene_pipeline_geneends.py ends here
