import os
import csv
import argparse
import glob
from pathlib import Path

# Create an argument parser
parser = argparse.ArgumentParser(description="Generate a CSV file from a directory of fastq files.")
parser.add_argument("glob_str", help="A glob expresion which should return all fastq files of interest (ex: output/*/reads_by_taxon/*.extracted_{1,2}.fq")
parser.add_argument("-fq1_ext", help=" forward fq file naming pattern", default="_R1.fq")
parser.add_argument("-fq2_ext", help=" reverse fq file naming pattern", default="_R2.fq")
parser.add_argument("-get_taxid", help="don't look for taxid on file names", action="store_true")
parser.add_argument("-filename_sep", default='.',
    help="separator char used to grab metadata from filenames (default='.')")
parser.add_argument("-manifest_out", help="path for the manifest csv", default=None)

# Parse the command-line arguments
args = parser.parse_args()
files_found = [os.path.abspath(fl) for fl in glob.glob(args.glob_str)]
ext_1 = args.fq1_ext
ext_2 = args.fq2_ext
get_taxid = args.get_taxid
mtdt_sep = args.filename_sep
mnf_out_path = args.manifest_out

# SANITY CHECK
if files_found == None:
    print("ERROR: no items provided for `--files`")
    print("       try to use a glob expression (ex: `--files mydir/*_{1,2}.fq`)")
    exit(1)

else:
    try:
        assert(len(files_found) > 0)
        print(f" > {len(files_found)} total files found")
    except(AssertionError):
        print("ERROR: no files provided")

'''
# FOR SOME REASON, DOES NOT WORK ON THE FARM =O. 
# exists() just return false even if the file exist on lustre
# check if all files exists
for fl in files_found:
    try:
        assert(os.path.exists(fl))
        #assert(fl.is_file())
    except(AssertionError):
        print(f"Error: '{fl}' does not exist.")
        print(fl.is_file())
        print(stop)

else:
'''

# Initialize an empty list to store the data
data = []

# Loop through files in the directory
idx = 0
for fl in files_found:
        # get filename
        split_path_lst = fl.split("/")
        source_dir = "/".join(split_path_lst[0:-1])
        filename = split_path_lst[-1]

        if filename.endswith(ext_1): #if ext1 exist, ext2 must exist
            fq1_flpath = fl
            # Check if pair exist
            fq2_filename = filename.replace(ext_1, ext_2)
            fq2_flpath = source_dir+"/"+fq2_filename
            
            try:
                assert(fq2_flpath in files_found)
            except(AssertionError):
                print(f"WARN: Skipping {fl}. No pair available")
                print(f"fq2 expected = {fq2_flpath}")
                continue
            
            # Extract sample_id and virus_id from the file name
            parts = filename.split(mtdt_sep)
            try:
                if get_taxid == False:
                    assert(len(parts) > 1)
                if get_taxid == True:
                    assert(len(parts) > 2)
            except(AssertionError):
                print(f"WARN: separator provided {mtdt_sep} doesn't return expected partition for {f}")
                print(f"      partition = {parts}")
                continue

            sample_id = f"{parts[0]}"
            # get meta data
            if "#" in sample_id:
                sample_id = sample_id.replace("#","-")

            if get_taxid == True:
                taxid = parts[1]

            # Append the data to the list
            if get_taxid == False:
                rows_values = [idx, sample_id, fq1_flpath, fq2_flpath]
            if get_taxid == True:
                rows_values = [idx, sample_id, fq1_flpath, fq2_flpath, taxid]

            data.append(rows_values)
            idx +=1

print(f" > {idx+1} samples on the manifest")
# Define the CSV output file
if mnf_out_path == None:
        output_file = os.path.join(os.getcwd(), 'manifest.csv')

if mnf_out_path != None:
        output_file = mnf_out_path
    
# SANITY CHECK
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