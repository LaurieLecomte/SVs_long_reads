#!/bin/bash

# Call SVs with nanovar across whole genome, including unplaced contigs which will be removed after

# Nanovar indexes are picky and will crash if other indexes are present in the genome directory, so we need to provide a genome directory in which Nanovar can add its own indexes

# Genome indexing steps prior to SV calling are very LONG, so we run run these steps only once for the fi
# valeria
#
# parallel -a 02_infos/ind_ONT.txt -j 4 srun -c 4 -p ibis_small --time=1-00:00:00 -J 02.1_svim_call_{} --mem=15G -o log/02.1_svim_call_{}_%j.log /bin/sh ./01_scripts/02.1_svim_call.sh {} &

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
if [[ ! -d "$CALLS_DIR/nanovar/$SAMPLE" ]]
then
  mkdir "$CALLS_DIR/nanovar/$SAMPLE"
fi

# 1. Run NanoVar
nanovar $BAM_DIR/"$SAMPLE".bam $GENOME_NV $CALLS_DIR/nanovar/$SAMPLE -x ont -t $CPU 