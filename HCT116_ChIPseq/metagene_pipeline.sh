#!/bin/bash

######## Define variables ########

### Main input variables

annotationbase=top_1000_filt

#dataset=totalpolII_ser7_only
dataset=totalpolII_ser7_IFN_SY_only
#dataset=all_samples

up3=1000 #upstream window for 3' metagene
down3=7500 #downstream window for 3' metagene
full3=8501 #full window size for 3' metagene (up + down + 1)
up5=1500 #upstream window for 5' metagene
down5=1500 #downstream window for 5' metagene
full5=3001 #full window size for 5' metagene (up + down + 1)

projectdir=/path/to/project/dir
processdir="$projectdir"/processing_files
scriptdir="$projectdir"/scripts/metagene

### I/O variables - uncomment correct set

indir=/path/to/bam/files
outdir="$projectdir"/metagenes/"$dataset"/"$annotationbase"
sizefactors="$processdir"/chipseq_norm_factors.txt
expt="$annotationbase"

scratch=/path/to/scratch/dir/"$dataset"_"$annotationbase"

### Variables set based on above variables

fullannotation="$processdir"/annotations/hg38_refseq_genenames_included.bed
annotation3="$scratch"/regions/"$annotationbase"_3prime_"$full3".bed
annotation5="$scratch"/regions/"$annotationbase"_5prime_"$full5".bed
annotationfull="$scratch"/regions/"$annotationbase"_fullgene_annotation.bed


######## Run scripts ########
### Run annotation creation script to generate correct windows for 5' and 3' metagenes
### (happens once for all samples)

# Set up initial files and directories
mkdir -p "$scratch"/regions "$outdir"

if [ -e "$annotation5" ] && [ -s "$annotation5" ]; then
  echo "annotation file already exists"

else
  rsync "$processdir"/annotations/"$annotationbase"_* "$scratch"/regions/
  rsync "$fullannotation" "$scratch"/regions/

  sbatch "$scriptdir"/metagene_pipeline_annotation_windows.sbatch \
    "$scratch"/regions \
    "$annotationbase" \
    "$annotation3" \
    "$annotation5" \
    "$up3" \
    "$down3" \
    "$up5" \
    "$down5"

  # Annotation creation is fast, so make sure it finishes
  sleep 30

  rsync "$annotation3" "$outdir"/
  rsync "$annotation5" "$outdir"/
  rsync "$annotationfull" "$outdir"/

fi

######## Run genomecov and python scripts for samples ########

### Run 5' metagene script
sbatch "$scriptdir"/metagene_pipeline_5prime.sbatch \
  "$scratch" \
  "$indir" \
  "$sizefactors" \
  "$annotation5" \
  "$full5" \
  "$projectdir" \
  "$outdir" \
  "$dataset"

# BAM creation doesn't take very long, so make sure no reading conflicts occur
sleep 210

### Run 3' metagene script
sbatch "$scriptdir"/metagene_pipeline_3prime.sbatch \
  "$scratch" \
  "$indir" \
  "$sizefactors" \
  "$annotation3" \
  "$full3" \
  "$projectdir" \
  "$outdir" \
  "$dataset"

# BAM creation doesn't take very long, so make sure no reading conflicts occur
sleep 180

### Run full gene processivity script
sbatch "$scriptdir"/metagene_pipeline_processivity.sbatch \
  "$scratch" \
  "$indir" \
  "$sizefactors" \
  "$annotationfull" \
  "$projectdir" \
  "$outdir" \
  "$dataset"
