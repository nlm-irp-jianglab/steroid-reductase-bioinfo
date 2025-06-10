#!/bin/bash

module load bowtie
module load seqkit
module load pigz
module load samtools
module load trimgalore
module load sratoolkit

SAMPLE=$1 # SRA run ID (e.g., SRR######)
NAME="5BR.genes.ffn"
IND=/gpfs/gsfs12/users/Irp-jiang/share/keith/mgx_prog/indices/ # path to bowtie2 index for gene (change gene to your file name)
THREADS=10 # Number of threads
OUT=/gpfs/gsfs12/users/Irp-jiang/share/keith/mgx_prog/out_2/ # path to output directory
TMP=/lscratch/$SLURM_JOB_ID # Path to temporary directory
HUMAN=/gpfs/gsfs12/users/Irp-jiang/share/DB_Share/human/chm13v2.0/chm13v2.0 #bowtie human reference path
#BOWTIE2_INDEXES=/data/dufaultthompskm/uroR/08.GTDB-marker-genes/bac120_marker_genes_reps_r214/fna

# Download reads and trim them using trimgalore
echo "Downloading reads"
echo $SAMPLE
fasterq-dump -v -e $THREADS -t $TMP -O /lscratch/$SLURM_JOB_ID/ $SAMPLE

echo "Downloaded files:"
ls /lscratch/$SLURM_JOB_ID/

pigz -p $THREADS /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fastq
pigz -p $THREADS /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fastq

echo "Trimming reads"
trim_galore -o /lscratch/$SLURM_JOB_ID/ -j $THREADS --paired /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fastq.gz /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fastq.gz

# Map reads to the human genome reference and remove likely contaminants
bowtie2 -p $THREADS -x $HUMAN -1 /lscratch/$SLURM_JOB_ID/${SAMPLE}_1_val_1.fq.gz \
    -2 /lscratch/$SLURM_JOB_ID/${SAMPLE}_2_val_2.fq.gz \
     | samtools view -bS | samtools fastq -@ $THREADS -f 12 -F 256 \
     -1 /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fq -2 /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fq -

#mv /lscratch/$SLURM_JOB_ID/${SAMPLE}_1_val_1.fq.gz /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fq.gz
#mv /lscratch/$SLURM_JOB_ID/${SAMPLE}_2_val_2.fq.gz /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fq.gz

pigz -p $THREADS /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fq
pigz -p $THREADS /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fq

# Summarize total number of reads in the QC'd reads
cat /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fq.gz | seqkit -j $THREADS stats -T > /lscratch/$SLURM_JOB_ID/${SAMPLE}.stats
cat /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fq.gz | seqkit -j $THREADS stats -T >> /lscratch/$SLURM_JOB_ID/${SAMPLE}.stats
TOTAL=$(head -n 2 /lscratch/$SLURM_JOB_ID/${SAMPLE}.stats | tail -n 1| cut -f 4)

# ../02.indices/3BHSDHI.ffn  ../02.indices/5BR.ffn  ../02.indices/BAIH.ffn
for NAME in "3BHSDHI.ffn" "5BR.ffn" "BAIH.ffn" ; do
    # Map reads to the bilirubin reductase index
    echo "Mapping Reads"
    bowtie2 --no-unal -p $THREADS -x ${IND}/${NAME} -1 /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fq.gz -2 /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fq.gz \
        | samtools view -bS | samtools sort > /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.bam
    COUNT=$(samtools idxstats /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.bam | awk '{ SUM += $3} END {print SUM}')
    samtools idxstats /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.bam > /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.genemap

    mv /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.mapping.summary ${OUT}/
    mv /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.genemap ${OUT}/

done