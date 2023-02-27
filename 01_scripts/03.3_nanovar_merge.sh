#!/bin/bash

# Merge refined sniffles SVs across samples with Jasmine

# manitou
# srun -c 8 -p medium --time=2-00:00:00 -J 03.3_nanovar_merge --mem=100G -o log/03.3_nanovar_merge_%j.log /bin/sh ./01_scripts/03.3_nanovar_merge.sh &

# valeria
# srun -c 8 -p ibis_medium --time=2-00:00:00 -J 03.3_nanovar_merge --mem=100G -o log/03.3_nanovar_merge_%j.log /bin/sh ./01_scripts/03.3_nanovar_merge.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=8

VCF_LIST_NANOVAR="$CALLS_DIR/nanovar/vcf_list.txt" # list of nanovar VCFs files

# 1. Make a list of nanovar sample VCF files to merge
ls -1 $CALLS_DIR/nanovar/*refined_dupToIns.vcf > $VCF_LIST_NANOVAR

# 2. Merge VCFs across samples 
jasmine file_list=$VCF_LIST_NANOVAR out_file="$MERGED_DIR/nanovar/nanovar_PASS_RSUPP2_refined.vcf" out_dir=$MERGED_DIR/nanovar genome_file=$GENOME --ignore_strand --mutual_distance --allow_intrasample --output_genotypes --threads=$CPU

# 3. Convert INSs back to DUPs (out_file is the VCF to be postprocessed, will be modified in situ)
jasmine out_file="$MERGED_DIR/nanovar/nanovar_PASS_RSUPP2_refined.vcf" out_dir=$MERGED_DIR/svim out_dir=$MERGED_DIR/nanovar genome_file=$GENOME --threads=$CPU --dup_to_ins --postprocess_only


# Clean up
#rm 