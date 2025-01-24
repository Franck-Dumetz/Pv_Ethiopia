## Checking for co-infection with different species of Plasmodium

The miochondrial genome is too close to between Pv and Pk to distinguish them. The distinction can only be done using the whole genome. <br />
Remove the mitochondrial genome from the fasta file of the 5 species and create a concatenated fasta file. <br /> 
Then use [Map_5sp_SlurmArray.sh](https://github.com/Franck-Dumetz/Pv_Ethiopia/blob/main/Co-infection_multiplePlasmo/Map_5sp_SlurmArray.sh) to map the reads to the 5sp genome fasta. Samtools excludes the combined flags Combined: 2048 + 4 + 256 = 2308:
  - 2048 (0x800): Suppress supplementary alignments.
  - 4 (0x004): Exclude unmapped reads.
  - 256 (0x100): Exclude secondary alignments. <br />
