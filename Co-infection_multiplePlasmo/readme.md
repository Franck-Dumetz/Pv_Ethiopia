## Checking for co-infection with different species of Plasmodium

The miochondrial genome is too close to between Pv and Pk to distinguish them. The distinction can only be done using the whole genome. <br />
Remove the mitochondrial genome from the fasta file of the 5 species and create a concatenated fasta file. <br /> 
Then use [Map_5sp_SlurmArray.sh](https://github.com/Franck-Dumetz/Pv_Ethiopia/blob/main/Co-infection_multiplePlasmo/Map_5sp_SlurmArray.sh) to map the reads to the 5sp genome fasta. Samtools excludes the combined flags Combined: 2048 + 4 + 256 = 2308:
  - 2048 (0x800): Suppress supplementary alignments.
  - 4 (0x004): Exclude unmapped reads.
  - 256 (0x100): Exclude secondary alignments. <br />

We need to count how manyreads maps to each portion go the different genomes
Firstly, use [count_specie.py](https://github.com/Franck-Dumetz/Pv_Ethiopia/blob/main/Co-infection_multiplePlasmo/count_species.py) to read through all the bam files and export csv file for each samples.
Secondly, extract all headers from initial fasta file to have a reference => what header for what species
```
grep '>' Pk_GCF_000006355_no_MIT.fasta | awk '{print $1}' | sed 's/^>//' > Pk_names.txt  
grep '>' PvivaxP01_genome_no_MIT.fasta | awk '{print $1}' | sed 's/^>//' > Pv_names.txt
grep '>' PmUG01_genome_no_MIT.fasta | awk '{print $1}' | sed 's/^>//' > Pm_names.txt
grep '>' Pfalciparum3D7_genome_no_MIT.fasta | awk '{print $1}' | sed 's/^>//' > Pf_names.txt
grep '>' PocGH01_genome_no_MIT.fasta | awk '{print $1}' | sed 's/^>//' > Po_names.txt
```
Finally, use [Species_RC_summary.py](https://github.com/Franck-Dumetz/Pv_Ethiopia/blob/main/Co-infection_multiplePlasmo/Species_RC_summary.py) to generate a summary table with all the samples with the number of reads per species and the percentage that value represents in the total sample. <br />

