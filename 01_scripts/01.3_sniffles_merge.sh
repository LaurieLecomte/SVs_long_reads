#!/bin/bash

# Merge refined sniffles SVs across samples with Jasmine

# manitou
# srun -c 10 -p medium --time=2-00:00:00 -J 01.3_sniffles_merge --mem=100G -o log/01.3_sniffles_merge_%j.log /bin/sh ./01_scripts/01.3_sniffles_merge.sh &

# valeria
# srun -c 10 -p ibis_medium --time=2-00:00:00 -J 01.3_sniffles_merge --mem=100G -o log/01.3_sniffles_merge_%j.log /bin/sh ./01_scripts/01.3_sniffles_merge.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=10

VCF_LIST_sniffles="$CALLS_DIR/sniffles/vcf_list.txt" # list of sniffles VCFs files
BAM_LIST="$CALLS_DIR/sniffles/bam_list.txt" # list of bam files for each sample


# 1. Make a list of sniffles VCF files to merge
ls -1 $CALLS_DIR/sniffles/*refined_dupToIns.vcf > $VCF_LIST_sniffles

# 2. Make a list of bam files
#ls -1 $BAM_DIR/*.bam > $BAM_LIST

# 3. Merge VCFs across samples 
jasmine file_list=$VCF_LIST_sniffles out_file="$MERGED_DIR/sniffles/sniffles_PASS_PRECISE_RSUPP2_refined.vcf" out_dir=$MERGED_DIR/sniffles genome_file=$GENOME --ignore_strand --ignore_merged_inputs --normalize_type --output_genotypes --allow_intrasample --mutual_distance --max_dist_linear=0.25 --threads=$CPU

#jasmine file_list=$VCF_LIST_sniffles out_file=$MERGED_DIR/sniffles/merged_refined_distlin0.1_dist50.vcf out_dir=$MERGED_DIR/sniffles genome_file=$GENOME bam_list=$BAM_LIST --mutual_distance --output_genotypes --normalize_type --allow_intrasample --ignore_strand threads=$CPU --max_dist_linear=0.1 --min_dist=50 


# Clean up
#rm 