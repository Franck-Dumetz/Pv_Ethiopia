#!/bin/bash

# Using featureCount to generate a count table

# Paths to resources
ANNOTATION="/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/PlasmoDB-68_PvivaxP01.gff"
OUTPUT="Pv_Ethiopia_PvAligned.txt"
THREADS=8

# Collect all BAM files from multiple directories
BAM_FILES=$(find /local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/bam_Pv_aligned/bam_dupRM -name "*.bam")

# Run featureCounts
featureCounts -a $ANNOTATION \
              -o $OUTPUT \
              -t exon \
              -g gene_id \
              -T $THREADS \
              $BAM_FILES
