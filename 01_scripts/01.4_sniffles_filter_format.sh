#!/bin/bash

# Re-filter SVs called by sniffles and format merged VCF before merging across callers 

# manitou
# srun -c 1 -p small --time=1-00:00:00 -J 01.4_sniffles_filter_format --mem=20G -o log/01.4_sniffles_filter_format_%j.log /bin/sh ./01_scripts/01.4_sniffles_filter_format.sh &

# valeria
# srun -c 1 -p ibis_small --time=1-00:00:00 -J 01.4_sniffles_filter_format --mem=20G -o log/01.4_sniffles_filter_format_%j.log /bin/sh ./01_scripts/01.4_sniffles_filter_format.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=1


# 1. Re-filter for PASS, PRECISE, number of supporting reads > 1 and wanted SVs types, then simplify VCF fields and sort
bcftools filter -i 'FILTER="PASS" & PRECISE=1 & SVTYPE!="BND" & SVTYPE!="INVDUP"' $MERGED_DIR/sniffles/sniffles_PASS_PRECISE_RSUPP2_refined.vcf | bcftools annotate -x ^INFO/SVTYPE,INFO/SVLEN,INFO/END,INFO/OLDTYPE | bcftools sort > $FILT_DIR/sniffles/sniffles_PASS_PRECISE_RSUPP2.vcf


# Clean up 
rm $MERGED_DIR/sniffles/sniffles_PASS_PRECISE_RSUPP2_refined_dupToIns.vcf
rm $MERGED_DIR/sniffles/sniffles_PASS_PRECISE_RSUPP2_refined_noGenotypes.vcf