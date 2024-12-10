#!/bin/bash

# Paths and settings
HISAT2_INDEX="/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/P01_hisat_index/PvivaxP01_hisat_index"  # Replace with your HISAT2 index base name
PARENT_DIR="/local/projects-t4/aberdeen2ro/SerreDLab-4/raw_reads/2024-12-10_Eugenia_Lo"  # Parent directory containing subfolders with FASTQ files
ALIGN_OUTPUT_DIR="bam_Pv_aligned"  # Directory for aligned BAM files
DEDUP_OUTPUT_DIR="bam_Pv_dedup"  # Directory for duplicate-removed BAM files
METRICS_DIR="bam_Pv_dedup/metrics"  # Directory for Picard metrics
PICARD_JAR="/usr/local/packages/picard-2.25.3/picard.jar"  # Path to Picard JAR file
THREADS=16  # Number of threads for HISAT2 and Picard

# Create necessary directories
mkdir -p "$ALIGN_OUTPUT_DIR"
mkdir -p "$DEDUP_OUTPUT_DIR"
mkdir -p "$METRICS_DIR"

# Loop through each subfolder in the parent directory
for SUBFOLDER in "$PARENT_DIR"/*; do
    if [ -d "$SUBFOLDER" ]; then
        echo "Processing folder: $SUBFOLDER"

        # Find paired FASTQ files
        FILE1=$(find "$SUBFOLDER" -type f -name "*_R1_trimmed.fastq.gz")
        FILE2=$(find "$SUBFOLDER" -type f -name "*_R2_trimmed.fastq.gz")

        # Check if both files exist
        if [[ -z "$FILE1" || -z "$FILE2" ]]; then
            echo "Skipping $SUBFOLDER: Missing paired files"
            continue
        fi

        # Extract the base folder name to use for output naming
        BASENAME=$(basename "$SUBFOLDER")

        echo "Processing $FILE1 and $FILE2..."

        # Step 1: Map the paired-end sequences using HISAT2
        /usr/local/packages/hisat2-2.2.1/hisat2 -x "$HISAT2_INDEX" \
                                                -1 "$FILE1" \
                                                -2 "$FILE2" \
                                                -S "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.sam" \
                                                --max-intronlen 5000 \
                                                -p "$THREADS" \
                                                2> "$ALIGN_OUTPUT_DIR/${BASENAME}_log.txt"

        echo "Mapping completed for $BASENAME. Converting SAM to BAM..."

        # Convert SAM to BAM and sort
        samtools view -bS "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.sam" > "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.bam"
        samtools sort "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.bam" -o "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv_sorted.bam"
        samtools index "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv_sorted.bam"

        # Remove intermediate SAM and unsorted BAM files
        rm "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.sam" "$ALIGN_OUTPUT_DIR/${BASENAME}_Pv.bam"

        # Step 2: Remove duplicates using Picard
        echo "Removing duplicates for $BASENAME..."

        java -jar "$PICARD_JAR" MarkDuplicates \
            INPUT="$ALIGN_OUTPUT_DIR/${BASENAME}_Pv_sorted.bam" \
            OUTPUT="$DEDUP_OUTPUT_DIR/${BASENAME}_Pv_dedup.bam" \
            METRICS_FILE="$METRICS_DIR/${BASENAME}_metrics.txt" \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=SILENT \
            TMP_DIR=/tmp

        # Index the deduplicated BAM file
        samtools index "$DEDUP_OUTPUT_DIR/${BASENAME}_Pv_dedup.bam"

        echo "Duplicate removal completed for $BASENAME."
    fi

done

echo "All processes completed. Aligned files saved in $ALIGN_OUTPUT_DIR and deduplicated files saved in $DEDUP_OUTPUT_DIR."
