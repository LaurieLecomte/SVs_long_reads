#!/bin/bash

# Re-filter SVs called by svim and format merged VCF before merging across callers 

# manitou
# srun -c 1 -p small --time=1-00:00:00 -J 02.4_svim_filter_format --mem=20G -o log/02.4_svim_filter_format_%j.log /bin/sh ./01_scripts/02.4_svim_filter_format.sh &

# valeria
# srun -c 1 -p ibis_small --time=1-00:00:00 -J 02.4_svim_filter_format --mem=20G -o log/02.4_svim_filter_format_%j.log /bin/sh ./01_scripts/02.4_svim_filter_format.sh &


# VARIABLES
GENOME="03_genome/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=1


# 1. Prepare file with new read names
## Extract old names
bcftools query -l $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.vcf > $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.old_names

## Extract new names
sed -E 's/[0-9]+\_([A-Za-z0-9]+)/\1/' $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.old_names > $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.new_names

## Paste old and new names in a single file to supply to bcftools reheader
paste -d "\t" $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.old_names $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.new_names > $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.old_new_names

# 2. Re-filter for PASS, number of supporting reads > 1 and wanted SVs types, then simplify VCF fields and sort
bcftools filter -i 'FILTER="PASS" & SVTYPE!="BND"' $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.vcf | bcftools annotate -x ^INFO/SVTYPE,INFO/SVLEN,INFO/END,INFO/OLDTYPE,INFO/IRIS_REFINED | bcftools reheader -s $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.old_new_names | bcftools sort > $FILT_DIR/svim/svim_PASS_RSUPP2.vcf

# Clean up 
rm $MERGED_DIR/svim/svim_PASS_RSUPP2_refined_dupToIns.vcf
rm $MERGED_DIR/svim/svim_PASS_RSUPP2_refined_noGenotypes.vcf
rm $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.old_names
rm $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.new_names
rm $MERGED_DIR/svim/svim_PASS_RSUPP2_refined.old_new_names