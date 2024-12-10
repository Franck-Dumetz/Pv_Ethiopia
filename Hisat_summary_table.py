import pandas as pd
import re
import os

def parse_alignment_log(log_file):
    """Parse a single alignment log file to extract metrics."""
    metrics = {}

    with open(log_file, 'r') as file:
        for line in file:
            if "reads; of these:" in line:
                total_reads = re.search(r"^(\d+)", line)
                metrics["Total Reads"] = int(total_reads.group(1)) if total_reads else "N/A"
            if "were paired; of these:" in line:
                paired_reads = re.search(r"\(([\d\.]+%)\)", line)
                metrics["Paired Reads (%)"] = paired_reads.group(1) if paired_reads else "N/A"
            if "aligned concordantly 0 times" in line:
                concordant_0 = re.search(r"\(([\d\.]+%)\)", line)
                metrics["Aligned Concordantly 0 Times (%)"] = concordant_0.group(1) if concordant_0 else "N/A"
            if "aligned concordantly exactly 1 time" in line:
                concordant_1 = re.search(r"\(([\d\.]+%)\)", line)
                metrics["Aligned Concordantly Exactly 1 Time (%)"] = concordant_1.group(1) if concordant_1 else "N/A"
            if "aligned concordantly >1 times" in line:
                concordant_gt1 = re.search(r"\(([\d\.]+%)\)", line)
                metrics["Aligned Concordantly >1 Times (%)"] = concordant_gt1.group(1) if concordant_gt1 else "N/A"
            if "aligned discordantly 1 time" in line:
                discordant = re.search(r"\(([\d\.]+%)\)", line)
                metrics["Aligned Discordantly 1 Time (%)"] = discordant.group(1) if discordant else "N/A"
            if "aligned 0 times" in line and "mates" in line:
                mates_0 = re.search(r"\(([\d\.]+%)\)", line)
                metrics["Aligned 0 Times (Mates) (%)"] = mates_0.group(1) if mates_0 else "N/A"
            if "aligned exactly 1 time" in line and "mates" in line:
                mates_1 = re.search(r"\(([\d\.]+%)\)", line)
                metrics["Aligned Exactly 1 Time (Mates) (%)"] = mates_1.group(1) if mates_1 else "N/A"
            if "aligned >1 times" in line and "mates" in line:
                mates_gt1 = re.search(r"\(([\d\.]+%)\)", line)
                metrics["Aligned >1 Times (Mates) (%)"] = mates_gt1.group(1) if mates_gt1 else "N/A"
            if "overall alignment rate" in line:
                overall_rate = re.search(r"([\d\.]+%)", line)
                metrics["Overall Alignment Rate (%)"] = overall_rate.group(1) if overall_rate else "N/A"

    return metrics

def summarize_logs(log_dir):
    """Summarize metrics from multiple log files into a single table."""
    summary_data = {}
    for log_file in os.listdir(log_dir):
        if log_file.endswith("_log.txt"):
            # Get file name without '_log.txt'
            sample_name = log_file.replace("_log.txt", "")
            log_path = os.path.join(log_dir, log_file)
            
            # Parse the log file
            metrics = parse_alignment_log(log_path)
            
            # Add to summary data
            summary_data[sample_name] = metrics

    # Convert summary data to a DataFrame
    summary_df = pd.DataFrame.from_dict(summary_data, orient="index")
    return summary_df

# Path to the directory containing log files
log_dir = "/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/bam_Pv_aligned/log_files"  # Replace with the directory containing log files

# Generate the summary table
summary_table = summarize_logs(log_dir)

# Save the summary table to a CSV file
summary_table.to_csv("/local/projects-t3/SerreDLab-3/fdumetz/Pv_ethiopia/bam_Pv_aligned/log_files/alignment_summary_table.csv")

# Display the summary table
print(summary_table)
