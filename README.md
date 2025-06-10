# steroid-reductase-bioinfo
Code relevant to the Arp et al. 2025 _Nature Communications_ paper, "Gut Bacteria Encode Reductases that Biotransform Steroid Hormones."

HMMs for steroid reductases are available in the ProkFunFind repository: ```https://github.com/nlm-irp-jianglab/ProkFunFind/tree/master/data/Steroid_reductase```

Steps to run the read mapping:

1. Clone this repository ```git clone https://github.com/nlm-irp-jianglab/steroid-reductase-bioinfo```
2. Change directory ```cd steroid-reductase-bioinfo/Scripts/```
3. Make script an executable ```chmod +x Map_Reads.sh```
4. Modify ```Map_Reads.sh``` to change the human index directory
5. Run ```./Map_Reads.sh {YOUR SRR NUMBER HERE}```. Your output table will be stored in ```steroid-reductase-bioinfo/Output```