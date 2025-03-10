#!/bin/bash
#SBATCH --job-name=metagene_3prime      # Job name
#SBATCH --mail-type=NONE         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=16              # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=04:00:00         # Time limit hrs:min:sec
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

endtype="3prime"
scratch=$1
indir=$2
sizefactors=$3
annotation=$4
fullwindow=$5
projectdir=$6
outdir="$7"/metagene_"$endtype"
dataset=$8
genome=<path to genome fasta>
chromsizes=<path to genome chromosome size file>
pyscript=<path to scripts>/metagene_pipeline_geneends.py

###################################################

mkdir -p "$scratch"/"$endtype" "$outdir"

for bamfile in "$indir"/<bam_prefix>*bam; do

  # Find basename
  bn=$(basename "$bamfile" <bam suffix>)
  basename="$scratch"/"$endtype"/"$bn"

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
    samtools view -@ 16 -h -P "$bamfile" -L "$annotation" \
      | samtools view -@ 16 -h -f 66 \
      | grep -P "(NH:i:1|^@)" \
      | awk 'BEGIN{FS=OFS="\t"}; \
             {if (NF < 9) {print $0} \
              else if (($9 % 2) != 0) \
                {if ($9 < 0) {newcol=($4+(($9-1)/2)); $4=newcol; print $0} \
                 else {newcol=($4+(($9+1)/2)); $4=newcol; print $0} \
                } \
              else {newcol=($4+($9/2)); $4=newcol; print $0} \
             }' \
      | samtools view -@ 16 -h -Sb - \
      > "$basename".bam

    samtools sort -@ 16 -o "$basename".sorted.bam "$basename".bam

    # Run genomecov to grab read coverage of 5' base of reads through those regions
    bedtools genomecov -dz -5 -ibam "$basename".sorted.bam \
      > "$basename".cov.prelim.bed

    # Fix columns for intersecting
    awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$2+1,$3}' "$basename".cov.prelim.bed \
      > "$basename".cov.bed

    # Run intersects to remove rest of genome and sort strands
    bedtools intersect -wa -wb \
      -a "$basename".cov.bed \
      -b "$annotation" \
      > "$basename".counts.bed

    rm "$basename".cov* "$basename"*bam*

  fi

done

# Run python script to generate binned combined counts for metagenes
# Bin sizes are defined within python script
python3 "$pyscript" \
  "$scratch"/"$endtype" \
  "$sizefactors" \
  "$fullwindow" \
  "$endtype" \
  "$dataset"

rsync "$scratch"/"$endtype"/*.txt "$outdir"/

# Final time stamp
echo Time is `date`
