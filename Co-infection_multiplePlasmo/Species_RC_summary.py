#!/bin/env python3

import os
import pandas as pd

# Directory containing RC files
input_folder = "/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/read_count_5sp"
output_file = "/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/5sp_RC_summary.csv"
pf_file = "/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/Pf_names.txt"
pv_file = "/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/Pv_names.txt"
pk_file = "/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/Pk_names.txt"
pm_file = "/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/Pm_names.txt"
po_file = "/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/Po_names.txt"

# Load species-specific chromosome names into sets
def load_species_names(file_path):
    with open(file_path, "r") as f:
        return set(line.strip() for line in f)

pf_names = load_species_names(pf_file)
pv_names = load_species_names(pv_file)
pk_names = load_species_names(pk_file)
pm_names = load_species_names(pm_file)
po_names = load_species_names(po_file)

# Initialize a dictionary to store results
summary_data = {
    "Pf": [],
    "Pv": [],
    "Pk": [],
    "Pm": [],
    "Po": []
}
files = []

# Process all RC files in the input folder
for rc_file in os.listdir(input_folder):
    if rc_file.endswith(".txt"):
        file_path = os.path.join(input_folder, rc_file)
        files.append(rc_file)

        # Read the RC file into a DataFrame
        df = pd.read_csv(file_path, sep="\t")

        # Initialize counters for each species
        species_counts = {
            "Pf": 0,
            "Pv": 0,
            "Pk": 0,
            "Pm": 0,
            "Po": 0
        }

        # Calculate total reads for each species
        for _, row in df.iterrows():
            species = row["Species"]
            count = row["Count"]

            if species in pf_names:
                species_counts["Pf"] += count
            elif species in pv_names:
                species_counts["Pv"] += count
            elif species in pk_names:
                species_counts["Pk"] += count
            elif species in pm_names:
                species_counts["Pm"] += count
            elif species in po_names:
                species_counts["Po"] += count

        # Append counts for each species
        for species, count in species_counts.items():
            summary_data[species].append(count)

# Create a DataFrame with species as rows and files as columns
summary_df = pd.DataFrame(summary_data, index=files).transpose()

# Calculate percentages for each sample
percentage_df = summary_df.div(summary_df.sum(axis=0), axis=1) * 100

# Interleave counts and percentages in the output
combined_data = []
for col in summary_df.columns:
    combined_data.append(summary_df[col])
    combined_data.append(percentage_df[col].rename(f"{col} (%)"))

combined_df = pd.concat(combined_data, axis=1)

# Write the combined DataFrame to a CSV file
combined_df.to_csv(output_file, index_label="Species")

print(f"Species read counts summary with percentages for all files has been written to {output_file}")
