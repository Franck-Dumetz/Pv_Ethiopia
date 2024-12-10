#!/bin/bash

# Set paths to your HISAT2 index, input parent directory, and output directory
HISAT2_INDEX="/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/P01_hisat_index/PvivaxP01_hisat_index"  # Replace with the base name of your HISAT2 index
PARENT_DIR="/local/projects-t4/aberdeen2ro/SerreDLab-4/raw_reads/2024-12-10_Eugenia_Lo"      # Parent directory containing subfolders with FASTQ files
OUTPUT_DIR="bam_Pv_aligned"          # Directory to store output SAM/BAM files
THREADS=16                            # Number of threads for HISAT2

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each subfolder in the parent directory
for SUBFOLDER in "$PARENT_DIR"/*; do
    if [ -d "$SUBFOLDER" ]; then
        echo "Processing folder: $SUBFOLDER"
        
        # Find paired FASTQ files in the subfolder
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
        
        # Map the paired-end sequences using HISAT2 and output a SAM file
        /usr/local/packages/hisat2-2.2.1/hisat2 -x "$HISAT2_INDEX" \
                                                -1 "$FILE1" \
                                                -2 "$FILE2" \
                                                -S "$OUTPUT_DIR/${BASENAME}_Pv.sam" \
                                                --max-intronlen 5000 \
                                                -p "$THREADS"\
                                                2> "$OUTPUT_DIR/${BASENAME}_log.txt"                             
          fi

echo "Mapping complete. Check logs in $OUTPUT_DIR for alignment details."
        # Convert SAM to BAM for smaller file size and indexing
        samtools view -bS "$OUTPUT_DIR/${BASENAME}_Pv.sam" > "$OUTPUT_DIR/${BASENAME}_Pv.bam"
        samtools sort "$OUTPUT_DIR/${BASENAME}_Pv.bam" -o "$OUTPUT_DIR/${BASENAME}_Pv_sorted.bam"
        samtools index "$OUTPUT_DIR/${BASENAME}_Pv_sorted.bam"
        
        # Remove intermediate SAM and unsorted BAM files to save space
        rm "$OUTPUT_DIR/${BASENAME}_Pv.sam" "$OUTPUT_DIR/${BASENAME}_Pv.bam"
done

echo "Mapping complete. Results saved in $OUTPUT_DIR"

