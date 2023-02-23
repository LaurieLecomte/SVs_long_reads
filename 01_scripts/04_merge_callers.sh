#!/bin/bash

# Merge SV calls across the 3 callers used
# Jasmine must be installed in current env or session

# valeria
# srun -c 4 -p ibis_small --mem=100G -J 04_merge_callers -o log/04_merge_callers_%j.log /bin/sh 01_scripts/04_merge_callers.sh &

# manitou
# srun -c 4 -p small --mem=100G -J 04_merge_callers -o log/04_merge_callers_%j.log /bin/sh 01_scripts/04_merge_callers.sh &
 
# VARIABLES
GENOME="03_genome/genome.fasta"
CALLS_DIR="05_calls"
FILT_DIR="07_filtered"

MERGED_UNION_DIR="08_merged_union"
FILT_UNION_DIR="09_filtered_union"

SNIFFLES_VCF="$FILT_DIR/sniffles/sniffles_PASS_PRECISE.vcf"
SVIM_VCF="$FILT_DIR/svim/svim_PASS.vcf"
NANOVAR_VCF="$FILT_DIR/nanovar/nanovar_PASS.vcf"

VCF_LIST="02_infos/callers_VCFs.txt"
MERGED_VCF="$MERGED_UNION_DIR/merged_sniffles_svim_nanovar.vcf"

REGIONS_EX="02_infos/excl_chrs.txt"

CPU=4


# 0. Remove VCF list from previous trials if any
if [[ -f $VCF_LIST ]]
then
  rm $VCF_LIST
fi

# 1. List all merged VCFs
touch $VCF_LIST
 
echo $SNIFFLES_VCF > $VCF_LIST
echo $SVIM_VCF >> $VCF_LIST
echo $NANOVAR_VCF >> $VCF_LIST

# 2. Merge SV calls accross samples, using predefined parameters (same as merging for the 3 short reads callers)
jasmine file_list=$VCF_LIST out_file=$MERGED_VCF out_dir=$MERGED_UNION_DIR genome_file=$GENOME --ignore_strand --ignore_merged_inputs --normalize_type --output_genotypes --allow_intrasample --mutual_distance --max_dist_linear=0.25 --threads=$CPU