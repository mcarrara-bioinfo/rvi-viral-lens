#!/usr/bin/env python3
import argparse

def read_trf(file_path):
    """Read a TRF file and return a dictionary of read identifiers with their corresponding lines."""
    reads = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('@'):  # Read identifier lines start with '@'
                read_id = line.strip()
                # Remove the trailing number (1 or 2) from the read identifier
                read_id_base = read_id.rsplit(' ', 1)[0]  # Get the part before the last space
                reads[read_id_base] = read_id  # Store the original read ID
    return reads

def generate_combined_trf(trf1_path, trf2_path, output_path, unpaired_path):
    """Generate a new TRF file that includes equivalent reads from both TRF files."""
    reads_trf1 = read_trf(trf1_path)
    reads_trf2 = read_trf(trf2_path)
    # store reads which lose their pairs
    unpaired_reads_fl = open(unpaired_path,"w")

    with open(output_path, 'w') as output_file:
        for read_id_base in set(reads_trf1.keys()).union(reads_trf2.keys()):
            if read_id_base in reads_trf1:
                output_file.write(reads_trf1[read_id_base] + '\n')
                try:
                    r2_sfx = reads_trf1[read_id_base].split(" ")[1].replace("1:","2:")
                    r2_read_equiv = f"{reads_trf1[read_id_base].split(' ')[0]} {r2_sfx}"
                    if r2_read_equiv not in reads_trf2:
                        #print(r2_read_equiv +" | "+reads_trf1[read_id_base])
                        output_file.write(r2_read_equiv + '\n')
                        unpaired_reads_fl.write(r2_read_equiv + '\n')
                except(IndexError):
                    print(f"WARN: {reads_trf1[read_id_base]} is not a read id. skipping it")

            if read_id_base in reads_trf2:
                output_file.write(reads_trf2[read_id_base] + '\n')
                try:
                    r1_sfx = reads_trf2[read_id_base].split(" ")[1].replace("2:","1:")
                    r1_read_equiv = f"{reads_trf2[read_id_base].split(' ')[0]} {r1_sfx}"
                    if r1_read_equiv not in reads_trf1:
                        #print(r1_read_equiv + " | " + reads_trf2[read_id_base])
                        output_file.write(r1_read_equiv + '\n')
                        unpaired_reads_fl.write(r1_read_equiv + '\n')
                except(IndexError):
                    print(f"WARN: {reads_trf2[read_id_base]} is not a read id. skipping it")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a combined TRF file from two TRF files.")
    parser.add_argument("-f1", "--trf1", required=True, help="Path to the first TRF file.")
    parser.add_argument("-f2", "--trf2", required=True, help="Path to the second TRF file.")
    parser.add_argument("-o", "--output", required=True, help="Output path for the combined TRF file.")
    parser.add_argument("-u", "--unpaired", required=True, help="Output path for the unpaired reads ids list.")

    args = parser.parse_args()

    generate_combined_trf(args.trf1, args.trf2, args.output, args.unpaired)

    print(f"Combined TRF file created at: {args.output}")