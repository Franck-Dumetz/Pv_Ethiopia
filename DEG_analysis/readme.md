## Differential gene expression analysis

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

We used cibersortX and used gene deconvolution to extrapolate the stage composition of each sample as in Tebben et al mSystem 2022(PMID: 35862820 PMCID: PMC9426464 DOI: 10.1128/msystems.00258-22).
Cybersort run: 
- batch correction mode: B-mode
- Run mode: relative
- 
