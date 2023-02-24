#!/bin/bash

# Call SVs with nanovar across whole genome, including unplaced contigs which will be removed after

# Nanovar indexes are picky and will crash if other indexes are present in the genome directory, so we need to provide a genome directory in which Nanovar can add its own indexes

# Genome indexing steps prior to SV calling are very LONG, so we run run these steps only once for the first sample, then we run others

# Run on MANITOU only in a conda env
# 
# parallel -a 02_infos/ind_ONT.txt -j 4 srun -c 10 -p medium --time=3-00:00:00 -J 03.1_nanovar_call_{} --mem=100G -o log/03.1_nanovar_call_{}_%j.log /bin/sh ./01_scripts/03.1_nanovar_call.sh {} &


# Possible adjustment : first run nanovar on a single sample to produce required indexes, then re-run for other (or all) samples

# VARIABLES
SAMPLE=$1

GENOME_NV="03_genome/genome_NV/genome.fasta"
CHR_BED="02_infos/chrs.bed.gz"
BAM_DIR="04_bam"
CALLS_DIR="05_calls"
MERGED_DIR="06_merged"
FILT_DIR="07_filtered"

CPU=10

# 0. Create output dir
#if [[ ! -d "$CALLS_DIR/nanovar/$SAMPLE" ]]
#then
#  mkdir "$CALLS_DIR/nanovar/$SAMPLE"
#fi

# 1. Run NanoVar
#nanovar $BAM_DIR/"$SAMPLE".bam $GENOME_NV $CALLS_DIR/nanovar/$SAMPLE -x ont -t $CPU 

# 2. Sort, remove SVs where END is < than POS (usually happens if a SV is at POS 1 on an uplaced contig), then compress and index
bcftools sort $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE".nanovar.pass.vcf | bcftools filter -e "POS > INFO/END" > $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_all_contigs.vcf
bgzip $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_all_contigs.vcf -f
tabix -p vcf $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_all_contigs.vcf.gz -f

# 3. Filter out unplaced contigs
bcftools view -R $CHR_BED $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_all_contigs.vcf.gz > $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE".vcf

# 4. Filter for PASS calls, SVs other than BNDs and supported by more than 2 reads
bcftools filter -i 'FILTER="PASS" & SVTYPE!="BND" & SR > 1' $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE".vcf > $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_PASS.vcf

# 5. Add read names
## Extract required info from VCF
bcftools query -f '%CHROM\t%POS\t%ID\t%END\n' $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_PASS.vcf > $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_PASS.table

## Extract and format reads from sv_support_reads.tsv
Rscript 01_scripts/utils/nanovar_add_rnames.R $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_PASS.table $CALLS_DIR/nanovar/$SAMPLE/sv_support_reads.tsv $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_PASS.annot

## Annotate VCF with read names information
bgzip $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_PASS.annot -f
tabix -s1 -b2 -e4 $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_PASS.annot.gz -f

bcftools annotate -a $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_PASS.annot.gz -h 02_infos/annot.hdr -c CHROM,POS,ID,END,RNAMES $CALLS_DIR/nanovar/$SAMPLE/"$SAMPLE"_PASS.vcf > $CALLS_DIR/nanovar/"$SAMPLE"_PASS.vcf
