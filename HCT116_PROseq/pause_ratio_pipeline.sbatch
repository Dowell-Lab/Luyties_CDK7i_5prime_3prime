#!/bin/bash
#SBATCH --job-name=pause_ratio      # Job name
#SBATCH --mail-type=NONE         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=32              # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=1:00:00         # Time limit hrs:min:sec
#SBATCH --partition=short
#SBATCH --mem=30gb               # Memory limit
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

module load python/3.6.3
module load R/3.6.1
export R_LIBS_USER=~/.R/3.6.0/

##### SPECIFY VARIABLES ###########################

scratch=<path to scratch directory> #Don't put trailing slash on this
projectdir=<path to project directory>
scriptdir="$projectdir"/scripts/pause_ratio
processdir="$projectdir"/processing_files
outdir="$projectdir"/pause_ratio
filelist=featurecounts_input_combined.txt

###################################################

mkdir -p "$scratch"/regions
#mkdir -p "$scratch"/bams
mkdir -p "$scratch"/results
mkdir -p "$outdir"

# Transfer annotation files and file list (and bam files if not on scratch)
rsync "$processdir"/annotation/hg38_refseq_pauseratio_pauseregion.gtf "$scratch"/regions/
rsync "$processdir"/annotation/hg38_refseq_pauseratio_restofgene.gtf "$scratch"/regions/
rsync "$processdir"/featurecounts/"$filelist" "$scratch"/

# Run featurecounts
Rscript "$scriptdir"/featurecounts_pauseratio.r "$scratch" "$filelist"

# Run accumulation/calc script
python3 ./pauseratio.py "$scratch" "$filelist"

# Sync results back
rsync -ar "$scratch"/results "$outdir"

# Final time stamp
echo Time is `date`
