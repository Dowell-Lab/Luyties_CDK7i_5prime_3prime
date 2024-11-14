#!/bin/bash
#SBATCH --job-name=metagene_annotation      # Job name
#SBATCH --mail-type=NONE         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1               # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=00:05:00         # Time limit hrs:min:sec
#SBATCH --partition=short
#SBATCH --mem=8gb               # Memory limit
#SBATCH --output=<path to outerr directory>/%x_%j.out
#SBATCH --error=<path to outerr directory>/%x_%j.err


##### OUTPUT RUN PARAMETERS #######################

# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors

##### LOAD MODULES AND RLIB LOCATION ##############

module load bedtools/2.28.0

##### SPECIFY VARIABLES ###########################

scratch=$1
annotationbase=$2
annotation3=$3
annotation5=$4
up3=$5
down3=$6
up5=$7
down5=$8
genome=<path to hg38 genome>/hg38.fa
chromsizes=<path to hg38 genome>/hg38.chrom.sizes
fullgene_annotation=hg38_refseq_genenames_included.bed

###################################################

######## 3' metagene annotation creation ########

# Pad window to each side and make sure no overlap with other genes or TES
bedtools slop -l "$up3" -r "$down3" -s \
    -i "$scratch"/"$annotationbase"_tes_annotation.bed \
    -g "$chromsizes" \
    > "$scratch"/tmp.bed

grep -v -- -1 "$scratch"/tmp.bed \
    > "$scratch"/tmp2.bed

bedtools intersect -wa \
    -a "$scratch"/"$fullgene_annotation" \
    -b "$scratch"/"$annotationbase"_tes_annotation.bed \
    > "$scratch"/"$fullgene_annotation"_subtract.bed

bedtools intersect -wa -v \
    -a "$scratch"/"$fullgene_annotation" \
    -b "$scratch"/"$fullgene_annotation"_subtract.bed \
    > "$scratch"/rest_of_genes_tocompare.bed

bedtools intersect -wa -v \
    -a "$scratch"/tmp2.bed \
    -b "$scratch"/rest_of_genes_tocompare.bed \
    > "$annotation3"

rm "$scratch"/tmp*
rm "$scratch"/"$fullgene_annotation"*
rm "$scratch"/rest_of_genes_tocompare.bed
rm "$scratch"/"$annotationbase"_tes_annotation.bed


######## 5' metagene annotation creation ########

# Pad window to each side and make sure no overlap with other genes or TSS
bedtools slop -l "$up5" -r "$down5" -s \
    -i "$scratch"/"$annotationbase"_tss_annotation.bed \
    -g "$chromsizes" \
    > "$scratch"/tmp.bed

grep -v -- -1 "$scratch"/tmp.bed \
    > "$scratch"/tmp2.bed

bedtools intersect -wa \
    -a "$scratch"/"$fullgene_annotation" \
    -b "$scratch"/"$annotationbase"_tss_annotation.bed \
    > "$scratch"/"$fullgene_annotation"_subtract.bed

bedtools intersect -wa -v \
    -a "$scratch"/"$fullgene_annotation" \
    -b "$scratch"/"$fullgene_annotation"_subtract.bed \
    > "$scratch"/rest_of_genes_tocompare.bed

bedtools intersect -wa -v \
    -a "$scratch"/tmp2.bed \
    -b "$scratch"/rest_of_genes_tocompare.bed \
    > "$annotation5"

rm "$scratch"/tmp*
rm "$scratch"/"$fullgene_annotation"*
rm "$scratch"/rest_of_genes_tocompare.bed
rm "$scratch"/"$annotationbase"_tss_annotation.bed


# Final time stamp
echo Time is `date`
