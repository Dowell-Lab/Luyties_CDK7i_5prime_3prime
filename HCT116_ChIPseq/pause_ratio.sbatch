#!/bin/bash
#SBATCH --job-name=pause_ratio      # Job name
#SBATCH --mail-type=NONE         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=32              # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=8:00:00         # Time limit hrs:min:sec
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

module load samtools/1.16.1
module load bedtools/2.28.0
module load python/3.6.3

##### SPECIFY VARIABLES ###########################

pause_window=200
projectdir=<path to project directory>
scratch=<path to scratch directory>/HCT116_ChIPseq/pause_ratio/pr"$pause_window" #Don't put trailing slash on this
scriptdir="$projectdir"/scripts/pause_ratio
pfilesdir="$projectdir"/processing_files
outdir="$projectdir"/pause_ratio/results_"$pause_window"
bamdir=<path to input directory>
normfactors="$pfilesdir"/chipseq_norm_factors.txt
genelist=top_1000_expr_filt.txt
annotation_pr=pauseratio_pauseregion_"$pause_window"bp.bed
annotation_rog=pauseratio_restofgene_"$pause_window"bp.bed

###################################################

mkdir -p "$scratch"/regions
#mkdir -p "$scratch"/bams
mkdir -p "$scratch"/results
mkdir -p "$outdir"

# Transfer annotation files and file list (and bam files if not on scratch)
rsync "$pfilesdir"/annotations/"$annotation_pr" "$scratch"/regions/
rsync "$pfilesdir"/annotations/"$annotation_rog" "$scratch"/regions/
rsync "$pfilesdir"/"$filelist" "$scratch"/
rsync "$pfilesdir"/"$genelist" "$scratch"/

for bamfile in "$bamdir"/5609_B2_*bam; do

  bn=$(basename "$bamfile" .sorted.bam_SBB.bam)
  basename="$scratch"/"$bn"
  locname=/localscratch/"$bn"

  if [ -e "$basename".counts.bed ] && [ -s "$basename".counts.bed ]; then
    continue

  else

    # Extract only reads in regions in 5' annotation file
    # Filtering steps as follows:
      # First read in pair only, only properly paired
      # No multimapped reads
      # Finds midpoints of fragments (awk statement)
        # Keep header intact (less than 9 fields)
        # Make odd fragment lengths even (adding 1 to positive, subtracting 1 from negative)
        # Adds fragment length/2 to leftmost coordinate
        # IMPORTANT! These are first reads, the coordinate is 5' on + strand and 3' on negative
          # Could correct for this on - strand by adding read length
          # But then to have the midpoint be counted in 5' genomecov, bedtools would add read length back
          # So here I keep it uncorrected
    samtools view -@ 16 -h -P "$bamfile" -L "$scratch"/regions/"$annotation_pr" \
      | samtools view -@ 16 -h -f 66 \
      | grep -P "(NH:i:1|^@)" \
      | awk 'BEGIN{FS=OFS="\t"};
             {if (NF < 9) {print $0} \
              else if (($9 % 2) != 0) \
                {if ($9 < 0) {newcol=($4+(($9-1)/2)); $4=newcol; print $0} \
                 else {newcol=($4+(($9+1)/2)); $4=newcol; print $0} \
                } \
              else {newcol=($4+($9/2)); $4=newcol; print $0} \
             }' \
      | samtools view -@ 16 -h -Sb - \
      > "$locname".bam

    samtools sort -@ 16 -o "$locname".sorted.bam "$locname".bam
    rm "$locname".bam

    bedtools coverage -counts -f 1 -s -sorted -a "$scratch"/regions/"$annotation_pr" -b "$locname".sorted.bam \
      > "$basename".cov.bed

  fi

done

bedtools coverage -counts -f 1 -s -sorted -a "$scratch"/regions/"$annotation_pr" \
  -b /localscratch/5609_B2_21.sorted.bam, /localscratch/5609_B2_22.sorted.bam, \
  /localscratch/5609_B2_23.sorted.bam, /localscratch/5609_B2_24.sorted.bam, \
  /localscratch/5609_B2_51.sorted.bam, /localscratch/5609_B2_52.sorted.bam, \
  /localscratch/5609_B2_53.sorted.bam, /localscratch/5609_B2_54.sorted.bam, \
  /localscratch/5609_B2_71.sorted.bam, /localscratch/5609_B2_72.sorted.bam, \
  /localscratch/5609_B2_73.sorted.bam, /localscratch/5609_B2_74.sorted.bam, \
  /localscratch/5609_B2_N1.sorted.bam, /localscratch/5609_B2_N2.sorted.bam, \
  /localscratch/5609_B2_N3.sorted.bam, /localscratch/5609_B2_N4.sorted.bam \
  > "$scratch"/all_samples.cov.bed

rm /localscratch/5609*bam

# Run accumulation/calc script
#python3 "$scriptdir"/pause_ratio.py "$scratch" "$filelist" "$normfactors" "$pause_window" # "$genelist"

# Sync results back
#rsync -ar "$scratch"/results/* "$outdir"/

# Final time stamp
echo Time is `date`
