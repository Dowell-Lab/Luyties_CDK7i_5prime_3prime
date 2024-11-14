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

module load samtools/1.16.1
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
pyscript="$projectdir"/scripts/metagene/metagene_pipeline_geneends.py

###################################################

mkdir -p "$scratch"/"$endtype" "$outdir"

for bamfile in "$indir"/5609_B2_*bam; do

  # Find sample
  bn=$(basename "$bamfile" .sorted.bam_SBB.bam)
  basename="$scratch"/"$endtype"/"$bn"
  locname=/localscratch/"$endtype"/"$bn"
  mkdir /localscratch/"$endtype"

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

    # Run genomecov to grab read coverage of 5' base of reads through those regions
    # Separately for + and - strands
    bedtools genomecov -dz -5 -ibam "$locname".sorted.bam \
      > "$locname".cov.prelim.bed

    # Fix columns for intersecting
    awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$2+1,$3}' "$locname".cov.prelim.bed \
      > "$locname".cov.bed

    # Run intersects to remove rest of genome and sort strands
    bedtools intersect -wa -wb \
      -a "$locname".cov.bed \
      -b "$annotation" \
      > "$basename".counts.bed

    rm "$locname"*

  fi

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
