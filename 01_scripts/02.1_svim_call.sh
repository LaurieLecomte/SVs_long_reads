#!/bin/bash

# Call SVs with svim2 across whole genome, including unplaced contigs which will be removed after

# valeria
# parallel -a 02_infos/ind_ONT.txt -j 4 srun -c 4 -p ibis_small --time=1-00:00:00 -J 02.1_svim_call_{} --mem=15G -o log/02.1_svim_call_{}_%j.log /bin/sh ./01_scripts/02.1_svim_call.sh {} &

# VARIABLES
SAMPLE=$1

GENOME="03_genome/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=4

# 0. Create output dir
if [[ ! -d "$CALLS_DIR/svim/$SAMPLE" ]]
then
  mkdir "$CALLS_DIR/svim/$SAMPLE"
fi

# 1. Call SVs in whole genome
svim alignment $CALLS_DIR/svim/$SAMPLE $BAM_DIR/"$SAMPLE".bam $GENOME --insertion_sequences --read_names --sample $SAMPLE --max_consensus_length=50000

# 2. Sort, compress and index
bcftools sort $CALLS_DIR/svim/$SAMPLE/variants.vcf > $CALLS_DIR/svim/$SAMPLE/"$SAMPLE"_all_contigs.vcf
bgzip $CALLS_DIR/svim/$SAMPLE/"$SAMPLE"_all_contigs.vcf
tabix -p vcf $CALLS_DIR/svim/$SAMPLE/"$SAMPLE"_all_contigs.vcf.gz -f

# 3. Filter out unplaced contigs
bcftools filter -R $CHR_BED $CALLS_DIR/svim/$SAMPLE/"$SAMPLE"_all_contigs.vcf.gz > $CALLS_DIR/svim/"$SAMPLE".vcf