#!/bin/bash
#SBATCH --job-name=metagene_processivity      # Job name
#SBATCH --mail-type=NONE         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=16              # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=06:00:00         # Time limit hrs:min:sec
#SBATCH --partition=short
#SBATCH --mem=8gb               # Memory limit
#SBATCH --output=/path/to/outerr/%x_%j.out
#SBATCH --error=/path/to/outerr/%x_%j.err


##### OUTPUT RUN PARAMETERS #######################

# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors

##### LOAD MODULES ##############

module load samtools/1.16.1
module load bedtools/2.28.0
module load python/3.6.3

##### SPECIFY VARIABLES ###########################

endtype="processivity"
scratch=$1
indir=$2
sizefactors=$3
annotation=$4
projectdir=$5
outdir="$6"/"$endtype"
genome=/path/to/genome/fasta
chromsizes=/path/to/genome/chromosome/sizes
pyscript="$projectdir"/scripts/metagene_pipeline_processivity.py

###################################################

mkdir -p "$scratch"/"$endtype" "$outdir"

for bamfile in "$indir"/*bam; do

  # Find basename
  bn=$(basename "$bamfile" .sorted.bam)
  basename="$scratch"/"$endtype"/"$bn"

  if [ -e "$basename".counts.bed ] && [ -s "$basename".counts.bed ]; then
    continue

  else
    # Extract only reads in regions in full gene annotation file
    samtools view -@ 16 -h "$bamfile" -L "$annotation" \
      | grep -P "(NH:i:1|^@)" \
      | samtools view -@ 16 -Sb -h - \
      > "$basename".bam
    samtools sort -@ 16 -o "$basename".sorted.bam "$basename".bam

    # Run genomecov to grab read coverage of 5' base of reads through those regions
    # Separately for + and - strands
    bedtools genomecov -dz -5 -ibam "$basename".sorted.bam -strand + \
      > "$basename".cov.pos_prelim.bed
    bedtools genomecov -dz -5 -ibam "$basename".sorted.bam -strand - \
      > "$basename".cov.neg_prelim.bed

    # Fix columns for intersecting 
    awk '{FS=OFS="\t"} {print $1,$2,$2+1,$3}' "$basename".cov.pos_prelim.bed \
      > "$basename".cov.pos.bed
    awk '{FS=OFS="\t"} {print $1,$2,$2+1,$3}' "$basename".cov.neg_prelim.bed \
      > "$basename".cov.neg.bed

    rm "$basename"*_prelim.bed

    # Run intersects to remove rest of genome and sort strands
    bedtools intersect -wa -wb \
      -a "$basename".cov.pos.bed \
      -b "$annotation" \
      > "$basename".int.bed

    awk '{if($10 == "+") print $0}' "$basename".int.bed > "$basename".cov.pos.bed

    bedtools intersect -wa -wb \
      -a "$basename".cov.neg.bed \
      -b "$annotation" \
      > "$basename".int.bed

    awk '{if($10 == "-") print $0}' "$basename".int.bed > "$basename".cov.neg.bed

    cat "$basename".cov.pos.bed "$basename".cov.neg.bed > "$basename".counts.bed

    rm "$basename".int.bed "$basename".cov*.bed
    rm "$basename"*bam

  fi

done

# Run python script to generate binned combined counts for metagenes
# Bin sizes are defined within python script
python3 "$pyscript" \
  "$scratch"/"$endtype" \
  "$sizefactors" \
  "$endtype"

rsync "$scratch"/"$endtype"/*.txt "$outdir"/

# Final time stamp
echo Time is `date`
