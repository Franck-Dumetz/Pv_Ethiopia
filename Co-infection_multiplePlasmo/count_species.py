#!/bin/env python3

import pysam
import sys
from collections import Counter

# Load your BAM file and FASTA file
bam_file = "/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/5sp_aligned/KMHC-384T_5spsorted.bam"  # Replace with your BAM file path
output_file = "/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/ReadCount_384T_5sp.txt"  # Output file for the table

# Count occurrences of values in the second column
def count_column2_values(bam_path):
    column2_counts = Counter()

    with pysam.AlignmentFile(bam_path, "rb", threads=4) as bam:  # Use multiple threads for efficiency
        for read in bam.fetch(until_eof=True):
            if not read.is_unmapped:  # Ignore unmapped reads
                reference_name = bam.get_reference_name(read.reference_id)
                column2_counts[reference_name] += 1

    return column2_counts

# Write results to a file
def write_to_file(output_path, column2_counts):
    with open(output_path, "w") as out_file:
        out_file.write("Value\tCount\n")
        for value, count in column2_counts.items():
            out_file.write(f"{value}\t{count}\n")

# Main execution
def main():
    try:
        column2_counts = count_column2_values(bam_file)
        write_to_file(output_file, column2_counts)
        print(f"Counts have been written to {output_file}")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
