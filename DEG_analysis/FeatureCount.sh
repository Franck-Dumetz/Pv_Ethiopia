#!/bin/bash
#SBATCH --job-name=featureCounts                      # Job name
#SBATCH --output=/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/featureCounts.out   # Standard output log
#SBATCH --error=/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/featureCounts.err    # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=fdumetz@som.umaryland.edu
#SBATCH --cpus-per-task=32				# Number of CPUs per task
#SBATCH --mem=30G                                        # Memory per node
#SBATCH --partition=all
# Using featureCount to generate a count table

# Paths to resources
ANNOTATION="/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/PlasmoDB-68_PvivaxP01.gff"
OUTPUT="Pv_Ethiopia_PvCount.txt"
THREADS=32

# Collect all BAM files from multiple directories
BAM_FILES=$(find /local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/bam_Pv2P01_dedup -name "*.bam")

# Run featureCounts
featureCounts -a $ANNOTATION \
              -o $OUTPUT \
              -t exon \
              -g gene_id \
              -T $THREADS \
              -p \
              -s 2 \
              -C --largestOverlap \
              -M --primary \
              $BAM_FILES
