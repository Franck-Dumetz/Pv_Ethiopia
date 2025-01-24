#Differential gene expression

Software requirements:
-  hisat2-2.2.1
-  picard
-  R
-  samtools-1.20
 
##Read mapping to the P. vivax reference genome
Build the hisat index 
```
hisat2-2.2.1/hisat2-buid PlasmoDB-68_PvivaxP01_Genome.fasta PvivaxP01_hisat_index
```
then use [Map2Pv.sh](https://github.com/Franck-Dumetz/Pv_Ethiopia/blob/main/DEG_analysis/Map2Pv.sh) to map using Hisat2 and remove duplicate using Picard
for many samples, use a slurm array and use [Hisat_summary_table.py](https://github.com/Franck-Dumetz/Pv_Ethiopia/blob/main/DEG_analysis/Hisat_summary_table.py) to generate a table summarasing all hisat log output in a table.

[FeatureCount.sh](https://github.com/Franck-Dumetz/Pv_Ethiopia/blob/main/DEG_analysis/FeatureCount.sh) with generate a count table to use in R for the DEG analysis


