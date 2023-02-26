#!/bin/bash

# Refine SVs called with sniffles2 in each sample
# Launch on valeria in module collection with jasmine and iris

# manitou
# parallel -a 02_infos/ind_ONT.txt -j 4 srun -c 10 -p medium --time=3-00:00:00 -J 01.2_sniffles_refine_{} --mem=80G -o log/01.2_sniffles_refine_{}_%j.log /bin/sh ./01_scripts/01.2_sniffles_refine.sh {} &

# valeria
# parallel -a 02_infos/ind_ONT.txt -j 4 srun -c 10 -p ibis_medium --time=7-00:00:00 -J 01.2_sniffles_refine_{} --mem=100G -o log/01.2_sniffles_refine_{}_%j.log /bin/sh ./01_scripts/01.2_sniffles_refine.sh {} &


# VARIABLES
SAMPLE=$1

GENOME="03_genome/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=10

# LOAD REQUIRED MODULES
#module load gcc python/3.10 bcftools/1.13 samtools/1.15 minimap2/2.24 blast+/2.13.0 bedtools/2.30.0 java/17.0.2 racon/1.4.13

# 1. Replace DUP:TANDEM by DUP (not mandatory since sniffles calls DUPs only - we keep this step so that VCFs are processed in the same way across callers)
sed -E 's/SVTYPE=DUP\:TANDEM/SVTYPE=DUP/' "$CALLS_DIR/sniffles/"$SAMPLE"_PASS_PRECISE.vcf" | sed -E 's/\<DUP\:TANDEM\>/DUP/' > "$CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE"_PASS_PRECISE_correctedDUPs.vcf"

# 2. Convert duplications to insertions temporarily
echo "$CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE"_PASS_PRECISE_correctedDUPs.vcf" > $CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE".txt

jasmine file_list=$CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE".txt out_dir=$CALLS_DIR/sniffles/$SAMPLE genome_file=$GENOME out_file=$CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE"_PASS_PRECISE_correctedDUPs_dupToIns.vcf --dup_to_ins --preprocess_only

# 3. Refine with iris
iris genome_in=$GENOME vcf_in=$CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE"_PASS_PRECISE_correctedDUPs_dupToIns.vcf reads_in=$BAM_DIR/"$SAMPLE".bam vcf_out=$CALLS_DIR/sniffles/"$SAMPLE"_PASS_PRECISE_refined_dupToIns.vcf --out_dir=$CALLS_DIR/sniffles/$SAMPLE --keep_long_variants --also_deletions --threads=$CPU


# Clean up 
rm $CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE"_PASS_PRECISE_correctedDUPs.vcf
rm $CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE".txt
rm $CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE"_PASS_PRECISE_correctedDUPs_dupToIns.vcf
rm $CALLS_DIR/sniffles/$SAMPLE/resultsstore.txt
#rm $CALLS_DIR/sniffles/$SAMPLE/"$SAMPLE"_dupToIns.txt
rm $CALLS_DIR/sniffles/$SAMPLE/list_dupToIns.txt