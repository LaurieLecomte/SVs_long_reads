#!/bin/bash

# Re-filter SVs called by svim before merging with other callers 

# valeria
# srun -c 1 -p ibis_small --time=1-00:00:00 -J 02.3_svim_filter --mem=20G -o log/02.3_svim_filter_%j.log /bin/sh ./01_scripts/02.3_svim_filter.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=1


# 1. Re-filter for PASS and wanted SVs types, then simplify VCF fields and sort
bcftools filter -i 'FILTER="PASS" & SVTYPE!="BND"' | bcftools annotate -x ^INFO/SVTYPE,INFO/SVLEN,INFO/END | bcftools sort > $FILT_DIR/svim/svim_PASS.vcf