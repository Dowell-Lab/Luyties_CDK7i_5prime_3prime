#!/bin/bash
#SBATCH --job-name=metagene_5prime      # Job name
#SBATCH --mail-type=NONE         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=16              # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=12:00:00         # Time limit hrs:min:sec
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

##### LOAD MODULES ##############

module load samtools/1.10
module load bedtools/2.28.0
module load python/3.6.3

##### SPECIFY VARIABLES ###########################

endtype="5prime"
scratch=$1
indir=$2
sizefactors=$3
annotation=$4
fullwindow=$5
projectdir=$6
outdir="$7"/metagene_"$endtype"
genome=<path to hg38 genome>/hg38.fa
chromsizes=<path to hg38 genome>/hg38.chrom.sizes
pyscript="$projectdir"/scripts/metagene/metagene_pipeline_geneends_v2.py

###################################################

mkdir -p "$scratch"/"$endtype" "$outdir"
mkdir -p /localscratch/"$endtype"

for bamfile in "$indir"/*bam; do

  # Find basename
  bn=$(basename "$bamfile" .sorted.bam)
  basename="$scratch"/"$endtype"/"$bn"
  locname=/localscratch/"$endtype"/"$bn"

  if [ -e "$basename".counts.bed ] && [ -s "$basename".counts.bed ]; then
    continue

  else
    # Extract only reads in regions in 5' annotation file
    samtools view -@ 16 -h "$bamfile" -L "$annotation" \
      | grep -P "(NH:i:1|^@)" \
      | samtools view -@ 16 -Sb -h - \
      > "$locname".bam
    samtools sort -@ 16 -o "$locname".sorted.bam "$locname".bam

    # Run genomecov to grab read coverage of 5' base of reads through those regions
    # Separately for + and - strands
    bedtools genomecov -dz -5 -ibam "$locname".sorted.bam -strand + \
      > "$locname".cov.pos_prelim.bed
    bedtools genomecov -dz -5 -ibam "$locname".sorted.bam -strand - \
      > "$locname".cov.neg_prelim.bed

    # Fix columns for intersecting 
    awk '{FS=OFS="\t"} {print $1,$2,$2+1,$3}' "$locname".cov.pos_prelim.bed \
      > "$locname".cov.pos.bed
    awk '{FS=OFS="\t"} {print $1,$2,$2+1,$3}' "$locname".cov.neg_prelim.bed \
      > "$locname".cov.neg.bed

    rm "$locname"*_prelim.bed

    # Run intersects to remove rest of genome and sort strands
    bedtools intersect -wa -wb \
      -a "$locname".cov.pos.bed \
      -b "$annotation" \
      > "$locname".int.bed

    awk '{if($10 == "+") print $0}' "$locname".int.bed > "$locname".cov.pos.bed
    awk '{if($10 == "-") print $0}' "$locname".int.bed > "$locname".cov.antipos.bed

    bedtools intersect -wa -wb \
      -a "$locname".cov.neg.bed \
      -b "$annotation" \
      > "$locname".int.bed

    awk '{if($10 == "-") print $0}' "$locname".int.bed > "$locname".cov.neg.bed
    awk '{if($10 == "+") print $0}' "$locname".int.bed > "$locname".cov.antineg.bed

    cat "$locname".cov.pos.bed "$locname".cov.neg.bed > "$basename".counts.bed
    cat "$locname".cov.antipos.bed "$locname".cov.antineg.bed > "$basename".antisense_counts.bed

  fi

  rm "$locname"*

done

rm -r /localscratch/"$endtype"

# Run python script to generate binned combined counts for metagenes
# Bin sizes are defined within python script
python3 "$pyscript" \
  "$scratch"/"$endtype" \
  "$sizefactors" \
  "$fullwindow" \
  "$endtype"

rsync "$scratch"/"$endtype"/*.txt "$outdir"/

# Final time stamp
echo Time is `date`
