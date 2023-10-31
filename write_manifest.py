import os
import csv
import argparse

# Create an argument parser
parser = argparse.ArgumentParser(description="Generate a CSV file from a directory of fastq files.")
parser.add_argument("directory", help="Directory path containing fastq files")
parser.add_argument("-fq1_ext", help=" forward fq file naming pattern", default="_R1.fq")
parser.add_argument("-fq2_ext", help=" reverse fq file naming pattern", default="_R2.fq")
parser.add_argument("-get_taxid", help="don't look for taxid on file names", action="store_true")

# Parse the command-line arguments
args = parser.parse_args()

# Get the absolute path of the provided directory
directory = os.path.abspath(args.directory)
ext_1 = args.fq1_ext
ext_2 = args.fq2_ext
get_taxid = args.get_taxid

# Check if the directory exists
if not os.path.isdir(directory):
    print("Error: The specified directory does not exist.")
else:
    # Initialize an empty list to store the data
    data = []

    # Loop through files in the directory
    idx = 0
    for filename in os.listdir(directory):
        if filename.endswith(ext_1) or filename.endswith(ext_2):
            # Extract sample_id and virus_id from the file name
            parts = filename.split('_')
            sample_id = parts[0]
            if get_taxid == False:
                virus_id = parts[1]

            # Determine if it's R1 or R2
            if filename.endswith(ext_1):
                r1_file = os.path.join(directory, filename)
                r2_file = os.path.join(directory, filename.replace(ext_1, ext_2))
            else:
                r2_file = os.path.join(directory, filename)
                r1_file = os.path.join(directory, filename.replace(ext_2, ext_1))

            # Append the data to the list

            if get_taxid == False:
                rows_values = [idx, sample_id, r1_file, r2_file]
            if get_taxid == True:
                rows_values = [idx, sample_id, r1_file, r2_file, virus_id]

            data.append(rows_values)
            idx +=1

    # Define the CSV output file
    output_file = os.path.join(directory, 'manifest.csv')

    if len(data) == 0:
        print("ERROR: No data was found")
        exit(1)

    # Write the data to a CSV file
    col_names = ["index","sample_id","reads_1","reads_2","taxid"]

    if get_taxid == False:
        col_names.remove("taxid")

    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(col_names)
        csv_writer.writerows(data)

    print(f"CSV file '{output_file}' has been generated.")