#!/bin/bash

# Path to input BAM files and output directory
INPUT_DIR="/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/bam_files"      # Replace with the directory containing BAM files
OUTPUT_DIR="/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/Dup_rm"              # Replace with the directory for duplicate-removed BAM files
METRICS_DIR="/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/Dup_rm/metrics"            # Directory to store metrics files for each BAM
PICARD_JAR="/usr/local/packages/picard-2.25.3/picard.jar"         # Path to Picard JAR file
THREADS=4                                # Number of threads for Picard (if multithreading is enabled)

# Create output and metrics directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$METRICS_DIR"

# Iterate over all BAM files in the input directory
for BAM_FILE in "$INPUT_DIR"/*.bam; do
    # Extract the base filename (without directory and extension)
    BASENAME=$(basename "$BAM_FILE" .bam)
    
    echo "Processing $BAM_FILE..."
    
    # Run Picard MarkDuplicates to remove duplicates
    java -jar "$PICARD_JAR" MarkDuplicates \
        INPUT="$BAM_FILE" \
        OUTPUT="$OUTPUT_DIR/${BASENAME}_dedup.bam" \
        METRICS_FILE="$METRICS_DIR/${BASENAME}_metrics.txt" \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=/tmp

    # Index the deduplicated BAM file
    samtools index "$OUTPUT_DIR/${BASENAME}_dedup.bam"
done

echo "Duplicate removal completed. Results saved in $OUTPUT_DIR"
