#!/bin/bash

# Prerequisites: 
# Bowtie2, seqkit, pigz, samtools, trimgalore, and sratoolkit

SAMPLE=$1 # SRA run ID (e.g., SRR######)
NAME="5BR.genes.ffn"
IND=../Indices # path to bowtie2 index
THREADS=10 # Number of threads
OUT=../Output # path to output directory
TMP=../tmp # Path to temporary directory
HUMAN=/gpfs/gsfs12/users/Irp-jiang/share/DB_Share/human/chm13v2.0/chm13v2.0 #bowtie human reference path
#BOWTIE2_INDEXES=/data/dufaultthompskm/uroR/08.GTDB-marker-genes/bac120_marker_genes_reps_r214/fna

mkdir -p $OUT
mkdir -p $TMP
# Download reads and trim them using trimgalore
echo "Downloading reads"
echo $SAMPLE
fasterq-dump -v -e $THREADS -t $TMP -O ${TMP}/ $SAMPLE

echo "Downloaded files:"
ls $TMP

pigz -p $THREADS ${TMP}/${SAMPLE}_1.fastq
pigz -p $THREADS ${TMP}/${SAMPLE}_2.fastq

echo "Trimming reads"
trim_galore -o ${TMP}/ -j $THREADS --paired ${TMP}/${SAMPLE}_1.fastq.gz ${TMP}/${SAMPLE}_2.fastq.gz

# Map reads to the human genome reference and remove likely contaminants
bowtie2 -p $THREADS -x $HUMAN -1 ${TMP}/${SAMPLE}_1_val_1.fq.gz \
    -2 ${TMP}/${SAMPLE}_2_val_2.fq.gz \
     | samtools view -bS | samtools fastq -@ $THREADS -f 12 -F 256 \
     -1 ${TMP}/${SAMPLE}_1.fq -2 ${TMP}/${SAMPLE}_2.fq -

#mv ${TMP}/${SAMPLE}_1_val_1.fq.gz ${TMP}/${SAMPLE}_1.fq.gz
#mv ${TMP}/${SAMPLE}_2_val_2.fq.gz ${TMP}/${SAMPLE}_2.fq.gz

pigz -p $THREADS ${TMP}/${SAMPLE}_1.fq
pigz -p $THREADS ${TMP}/${SAMPLE}_2.fq

# Summarize total number of reads in the QC'd reads
cat ${TMP}/${SAMPLE}_1.fq.gz | seqkit -j $THREADS stats -T > ${TMP}/${SAMPLE}.stats
cat ${TMP}/${SAMPLE}_2.fq.gz | seqkit -j $THREADS stats -T >> ${TMP}/${SAMPLE}.stats
TOTAL=$(head -n 2 ${TMP}/${SAMPLE}.stats | tail -n 1| cut -f 4)

# ../02.indices/3BHSDHI.ffn  ../02.indices/5BR.ffn  ../02.indices/BAIH.ffn
for NAME in "3-beta_HSDH-I" "5-beta_reductase" "delta6-reductase" ; do
    # Map reads to the reductase index
    echo "Mapping Reads"
    bowtie2 --no-unal -p $THREADS -x ${IND}/${NAME} -1 ${TMP}/${SAMPLE}_1.fq.gz -2 ${TMP}/${SAMPLE}_2.fq.gz \
        | samtools view -bS | samtools sort > ${TMP}/${SAMPLE}.${NAME}.bam
    COUNT=$(samtools idxstats ${TMP}/${SAMPLE}.${NAME}.bam | awk '{ SUM += $3} END {print SUM}')
    samtools idxstats ${TMP}/${SAMPLE}.${NAME}.bam > ${TMP}/${SAMPLE}.${NAME}.genemap

    mv ${TMP}/${SAMPLE}.${NAME}.mapping.summary ${OUT}/
    mv ${TMP}/${SAMPLE}.${NAME}.genemap ${OUT}/

done