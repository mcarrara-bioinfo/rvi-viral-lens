import os
import csv
import argparse

# Create an argument parser
parser = argparse.ArgumentParser(description="Generate a CSV file from a directory of fastq files.")
parser.add_argument("directory", help="Directory path containing fastq files")

# Parse the command-line arguments
args = parser.parse_args()

# Get the absolute path of the provided directory
directory = os.path.abspath(args.directory)

# Check if the directory exists
if not os.path.isdir(directory):
    print("Error: The specified directory does not exist.")
else:
    # Initialize an empty list to store the data
    data = []

    # Loop through files in the directory
    idx = 0
    for filename in os.listdir(directory):
        if filename.endswith("_R1.fq.gz") or filename.endswith("_R2.fq.gz"):
            # Extract sample_id and virus_id from the file name
            parts = filename.split('_')
            sample_id = parts[0]
            virus_id = parts[1]

            # Determine if it's R1 or R2
            if filename.endswith("_R1.fq.gz"):
                r1_file = os.path.join(directory, filename)
                r2_file = os.path.join(directory, filename.replace("_R1.fq.gz", "_R2.fq.gz"))
            else:
                r2_file = os.path.join(directory, filename)
                r1_file = os.path.join(directory, filename.replace("_R2.fq.gz", "_R1.fq.gz"))

            # Append the data to the list
            data.append([idx, sample_id, r1_file, r2_file, virus_id])
            idx +=1 
    # Define the CSV output file
    output_file = os.path.join(directory, 'manifest.csv')

    # Write the data to a CSV file
    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["index", "sample_id", "fastq_file_r1", "fastq_file_r2", "virus_id"])
        csv_writer.writerows(data)

    print(f"CSV file '{output_file}' has been generated.")