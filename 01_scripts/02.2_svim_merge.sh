#!/bin/bash

# Merge SVs called by svim across samples with Jasmine and refine INS/DEL sequences with Iris


# valeria
# srun -c 10 -p ibis_medium --time=7-00:00:00 -J 02.2_svim_merge --mem=100G -o log/02.2_svim_merge_%j.log /bin/sh ./01_scripts/02.2_svim_merge.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=10

VCF_LIST_SVIM="$CALLS_DIR/svim/vcf_list.txt" # list of svim VCFs files
BAM_LIST="$CALLS_DIR/svim/bam_list.txt" # list of bam files for each sample

# 1. Make a list of svim VCF files to merge
ls -1 $CALLS_DIR/svim/*.vcf > $VCF_LIST_SVIM

# 2. Make a list of bam files
ls -1 $BAM_DIR/*.bam > $BAM_LIST

# 3. Merge samples VCFs and refine ALT sequences with Iris 
jasmine file_list=$VCF_LIST_svim out_file=$MERGED_DIR/svim/merged_refined_distlin0.1_dist50.vcf out_dir=$MERGED_DIR/svim genome_file=$GENOME bam_list=$BAM_LIST --mutual_distance --output_genotypes --normalize_type --allow_intrasample --ignore_strand threads=$CPU --run_iris iris_args=--keep_long_variants --max_dist_linear=0.1 --min_dist=50