# steroid-reductase-bioinfo
Code relevant to the Arp et al. 2025 _Nature Communications_ paper, "Gut Bacteria Encode Reductases that Biotransform Steroid Hormones."

HMMs for steroid reductases are available in the ProkFunFind repository: ```https://github.com/nlm-irp-jianglab/ProkFunFind/tree/master/data/Steroid_reductase```

## Read Mapping
### Installation Requirements
The following tools must be installed to replicate the read mapping step.
+ bowtie2
+ seqkit
+ pigz
+ samtools
+ trimgalore
+ sratoolkit

### Steps to map reads to the steroid reductase index

1. Clone this repository ```git clone https://github.com/nlm-irp-jianglab/steroid-reductase-bioinfo```
2. Change directory ```cd steroid-reductase-bioinfo/Scripts/```
3. Make script an executable ```chmod +x Map_Reads.sh```
4. Modify ```Map_Reads.sh``` to change the human index directory
5. Run ```./Map_Reads.sh {YOUR SRR NUMBER HERE}```. Your output table will be stored in ```steroid-reductase-bioinfo/Output```

The script will map the reads of a metagenomics sample (SRR) to our 3 bowtie2 indices: Δ4-3-ketosteroid 5β-reductase, 3β-hydroxysteroid dehydrogenase/Δ5-4 isomerase, and Δ6-3-ketosteroid reductase.