#!/bin/bash

# Call SVs with nanovar across whole genome, including unplaced contigs which will be removed after

# Nanovar indexes are picky and will crash if other indexes are present in the genome directory, so we need to provide a genome directory in which Nanovar can add its own indexes

# Genome indexing steps prior to SV calling are very LONG, so we run run these steps only once for the first sample, then we run others

# Run on MANITOU only in a conda env
# 
# parallel -a 02_infos/ind_ONT.txt -j 4 srun -c 10 -p medium --time=3-00:00:00 -J 03.1_nanovar_call_{} --mem=100G -o log/03.1_nanovar_call_{}_%j.log /bin/sh ./01_scripts/03.1_nanovar_call.sh {} &

# VARIABLES
SAMPLE=$1

GENOME_NV="03_genome/genome_NV/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=10

# 0. Create output dir
if [[ ! -d "$CALLS_DIR/nanovar/$SAMPLE" ]]
then
  mkdir "$CALLS_DIR/nanovar/$SAMPLE"
fi

# 1. Run NanoVar
nanovar $BAM_DIR/"$SAMPLE".bam $GENOME_NV $CALLS_DIR/nanovar/$SAMPLE -x ont -t $CPU 

# 2. Sort, remove SVs where END is < than POS, then compress and index
bcftools sort $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE".nanovar.pass.vcf | bcftools filter -e "POS > INFO/END" > $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_all_contigs.vcf
bgzip $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_all_contigs.vcf
tabix -p vcf $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_all_contigs.vcf.gz -f

# 3. Filter out unplaced contigs
bcftools view -R $CHR_BED $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_all_contigs.vcf.gz > $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE".vcf

# 3. Filter for PASS and PRECISE calls (should not do anything filter on PRECISE only because we use .pass vcf at previous step) and remove BNDs
bcftools filter -i 'FILTER="PASS" & SVTYPE!="BND"' $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE".vcf > $CALLS_DIR/nanovar/"$SAMPLE"_PASS.vcf