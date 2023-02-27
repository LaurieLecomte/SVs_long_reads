#!/bin/bash

# Refine SVs called with svim2 in each sample
# Launch on valeria in module collection with jasmine and iris


#srun -c 4 -p ibis_small --time=1-00:00:00 -J 02.3.0_svim_merge --mem=100G -o log/02.3.0_svim_merge_%j.log /bin/sh ./01_scripts/02.3.0_svim_merge.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=4

VCF_LIST_SVIM="$CALLS_DIR/svim/vcf_list.txt" # list of sniffles VCFs files

# LOAD REQUIRED MODULES
#module load gcc python/3.10 bcftools/1.13 samtools/1.15 minimap2/2.24 blast+/2.13.0 bedtools/2.30.0 java/17.0.2 racon/1.4.13

# 1. Make a list of VCFs to merge 
ls -1 $CALLS_DIR/svim/*_10000_PASS_refined_sed.vcf > $VCF_LIST_SVIM

# after at the merging step, we'll convert back converted DUPs to DUPs :
#Merge SVs across samples: 
#jasmine file_list=$VCF_LIST_SVIM out_file=$MERGED_DIR/svim/merged_14062_13070_nopost.vcf genome_file=$GENOME --ignore_strand --mutual_distance --allow_intrasample --output_genotypes --threads=$CPU --dup_to_ins

#jasmine file_list=$VCF_LIST_SVIM out_file=$MERGED_DIR/svim/merged_14062_13070.vcf genome_file=$GENOME --ignore_strand --mutual_distance --allow_intrasample --output_genotypes 

# 2. First merge 
jasmine file_list=$VCF_LIST_SVIM out_file=$MERGED_DIR/svim/merged_14062_13070.vcf genome_file=$GENOME --ignore_strand --mutual_distance --allow_intrasample --output_genotypes --threads=$CPU

# 3. Then convert insertions back to duplications
jasmine out_file=$MERGED_DIR/svim/merged_14062_13070.vcf genome_file=$GENOME --threads=$CPU --dup_to_ins --postprocess_only

# merge without iris
ls -1 $CALLS_DIR/svim/*_PASS_10000.vcf > $CALLS_DIR/svim/vcf_list_noIris.txt

jasmine file_list=$CALLS_DIR/svim/vcf_list_noIris.txt out_file=$MERGED_DIR/svim/merged_14062_13070_noIris.vcf genome_file=$GENOME --ignore_strand --mutual_distance --allow_intrasample --output_genotypes --threads=$CPU