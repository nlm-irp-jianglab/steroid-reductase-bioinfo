# Bioinformatics Analysis for "Gut Bacteria Encode Reductases that Biotransform Steroid Hormones"
Code relevant to the Arp et al. 2025 _Nature Communications_ paper, "Gut Bacteria Encode Reductases that Biotransform Steroid Hormones."

HMMs for steroid reductases are available in the ProkFunFind repository: ```https://github.com/nlm-irp-jianglab/ProkFunFind/tree/master/data/Steroid_reductase```

## Installation
### Operating System
This read mapping was tested on Rocky Linux 8.7 (Green Obsidian)

### Installation Requirements
The following tools must be installed to replicate the read mapping step. (Versions tested on)
+ bowtie2 (v2.5.3)
+ seqkit (v2.7.0)
+ pigz (v2.8)
+ samtools (v1.21)
+ trimgalore (v0.6.7)
+ sratoolkit (v3.2.1)

## Steps to map reads to the steroid reductase index

1. Clone this repository ```git clone https://github.com/nlm-irp-jianglab/steroid-reductase-bioinfo```
2. Change directory ```cd steroid-reductase-bioinfo/Scripts/```
3. Make script an executable ```chmod +x Map_Reads.sh```
4. Modify ```Map_Reads.sh``` to change the human index directory
5. Run ```./Map_Reads.sh {YOUR SRR NUMBER HERE}``` (e.g. ERR636349). Your output table will be stored in ```steroid-reductase-bioinfo/Output```

The script will map the reads of a metagenomics sample (SRR) to our 3 bowtie2 indices: Δ4-3-ketosteroid 5β-reductase, 3β-hydroxysteroid dehydrogenase/Δ5-4 isomerase, and Δ6-3-ketosteroid reductase.

### Output
There will be 3 outputs files for each ERR number, showing the number of reads mapped to each gene in the bowtie2 index. An example for ERR636349 is in the Output Folder. There is a file for each gene identified in the study (3-beta HSDH-I, 5-beta reductase, and delta 6-reductase)

+ ```steroid-reductase-bioinfo/Output/ERR636349.3-beta_HSDH-I.genemap```
+ ```steroid-reductase-bioinfo/Output/ERR636349.5-beta_reductase.genemap```
+ ```steroid-reductase-bioinfo/Output/ERR636349.delta6-reductase.genemap```

The columns of the output tables are:
1. Reference sequence name
2. Reference sequence length (in bp)
3. Number of reads mapped to this reference
4. Number of reads *unmapped* but whose mate maps to this reference

## Reproducibility
All SRR Accession numbers for the 1549 metagenomic samples used in the read mapping are listed in Supplementary Table 3.