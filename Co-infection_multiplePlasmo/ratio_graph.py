#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

# Load the second dataset from the new file
file_path_2 = '/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/Book4.txt'
data_2 = pd.read_csv(file_path_2, sep='\t')

# Save the plot as a PDF
output_pdf_path = '/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/MultiSpecies/Ratio_species.Pv.pdf'

# Set the BAM File column as the index for better labeling in the plot
data_2.set_index('BAM File', inplace=True)

# Plotting the stacked bar graph for the new dataset
data_2.plot(kind='bar', stacked=True, figsize=(15, 8))

# Adding labels and title
#plt.title('Stacked Bar Graph of Pf/Pv, Pk/Pv, Pm/Pv, Po/Pv Ratios', fontsize=16)
plt.xlabel('BAM File', fontsize=14)
plt.ylabel('Ratio', fontsize=14)
plt.legend(title="Ratios", fontsize=10)
plt.xticks(rotation=45, ha='right')

# Save the figure
plt.savefig(output_pdf_path)

print(f"Graph saved to: {output_pdf_path}")