#!/bin/bash

# Call SVs with sniffles2 across whole genome, including unplaced contigs which will be removed after

# valeria
# parallel -a 02_infos/ind_ONT.txt -j 4 srun -c 4 -p ibis_small --time=1-00:00:00 -J 01.1_sniffles_call_{} --mem=15G -o log/01.1_sniffles_call_{}_%j.log /bin/sh ./01_scripts/01.1_sniffles_call.sh {} &

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
if [[ ! -d "$CALLS_DIR/sniffles/$SAMPLE" ]]
then
  mkdir "$CALLS_DIR/sniffles/$SAMPLE"
fi

# 1. Call SVs in whole genome
#sniffles --input $BAM_DIR/"$SAMPLE".bam --vcf $CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE"_all_contigs.vcf.gz --snf $CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE".snf --threads $CPU --reference $GENOME --sample-id $SAMPLE --output-rnames --combine-consensus --allow-overwrite

# 2. Sort, remove SVs where END is < than POS (usually happens if a SV is at POS 1 on an uplaced contig) and remove unplaced contigs
bcftools view -R $CHR_BED $CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE"_all_contigs.vcf.gz | bcftools filter -e "POS > INFO/END" > $CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE".vcf

# 3. Filter for PASS and PRECISE calls and remove BNDs and INVDUPs and SVs supported by less than 2 reads
bcftools filter -i 'FILTER="PASS" & PRECISE=1 & SVTYPE!="BND" & SVTYPE!="INVDUP" & SUPPORT > 1' $CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE".vcf > $CALLS_DIR/sniffles/"$SAMPLE"_PASS_PRECISE.vcf