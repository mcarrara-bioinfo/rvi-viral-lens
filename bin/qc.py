#!/usr/bin/env python3

from Bio import SeqIO
import csv
import argparse

def get_covered_positions_from_file(depth_file_path: str, min_depth: int) -> int:
    """
    Parses a tab-delimited file of per-position depths and returns the total number 
    positions where the depth (3rd column) exceeds or is equal to a minimum depth 
    threshold value.
    """
    counter = 0

    with open(depth_file_path, "r") as depth_file_reader:
        for line in depth_file_reader:
            depth = line.split("\t")[2]

            if int(depth) >= min_depth:
                counter = counter + 1

    return counter

def get_N_positions(fasta: SeqIO) -> list:
    """
    Returns a list of all 'n' positions in a FASTA SeqIO object
    """
    n_pos = [pos for pos, letter in enumerate(fasta.seq.lower()) if letter == 'n']

    return n_pos

def get_pct_N_bases(fasta: SeqIO) -> float:
    """
    Determines the % of bases of a FASTA that are 'n'.
    """
    count_N = len(get_N_positions(fasta))

    pct_N_bases = count_N / len(fasta.seq) * 100

    return pct_N_bases

def get_largest_N_gap(fasta: SeqIO) -> list:
    """
    Returns the largest a list of all 'n' positions in a FASTA SeqIO object
    """
    # List of all 'n' positions in the FASTA
    n_pos = get_N_positions(fasta)
    n_pos = [0] + n_pos + [len(fasta.seq)]

    # From both start and end, simultaneously iterate through all positions and get all of the gaps between these values
    n_gaps = [j-i for i, j in zip(n_pos[:-1], n_pos[1:])]

    # Return largest 'n' gap
    return sorted(n_gaps)[-1]

def assess_depths_and_consensus(min_depth: int, fasta_file: str, ref_file: str, depths_file: str) -> dict:
    """
    Generates a dictionary of QC values
        - pct_N_bases, largest_N_gap, qc_pass are calculated from the consensus FASTA
        - pct_covered_bases is calculated from the depths file and the reference genome FASTA
    """
    # Default values for several consensus FASTA metrics
    pct_N_bases   = 0
    largest_N_gap = 0
    qc_pass       = "FALSE"

    # Read in FASTA and if not empty then perform QC calculations
    fasta = SeqIO.read(fasta_file, "fasta")
    if len(fasta.seq) != 0:

        pct_N_bases = get_pct_N_bases(fasta)
        largest_N_gap = get_largest_N_gap(fasta)

        # QC PASS / FAIL
        if largest_N_gap >= 10000 or pct_N_bases < 50.0:
            qc_pass = "TRUE"

    # Used by depth calculations
    depth_covered_bases = get_covered_positions_from_file(depths_file, min_depth)

    # Depth calculations
    ref_length = len(SeqIO.read(ref_file, "fasta").seq)
    pct_covered_bases = depth_covered_bases / ref_length * 100

    # Assemble dictionary
    # The order of keys is important
    pairs = [('pct_N_bases', "{:.2f}".format(pct_N_bases)),
             ('pct_covered_bases', "{:.2f}".format(pct_covered_bases)),
             ('longest_no_N_run', largest_N_gap),
             ('fasta', fasta_file),
             ('qc_pass', qc_pass)]
        
    return dict(pairs)

def generate_qc_file(args: argparse.ArgumentParser.parse_args):
    # Get QC values for a pair of bam-fasta files
    qc_values = assess_depths_and_consensus(args.minimum_depth, args.fasta, args.ref, args.depths_file)

    # Get the original dictionary keys in the order they were inserted
    column_names = list(qc_values)

    # Parse input flagstat file
    with open(args.flagstat_file) as flagstat_reader:
        flagstat_values = [line.split(" ")[0] for line in flagstat_reader.readlines()]

    # Add columns from flagstat file to QC dict
    qc_values['num_aligned_reads'] = flagstat_values[5]
    qc_values['total_mapped_reads'] = flagstat_values[4]
    qc_values['total_unmapped_reads'] = flagstat_values[10]

    # Add non-QC columns to QC dict
    qc_values['sample_name'] = args.sample
    qc_values['bam'] = args.bam

    # Add columns to header at specific indices
    column_names.insert(0, 'sample_name')
    column_names.insert(4, 'num_aligned_reads')
    column_names.insert(6, 'bam')
    column_names.extend(['total_mapped_reads', 'total_unmapped_reads'])

    # Add ivar_md column if correct flag supplied
    if args.ivar_md != None:
        qc_values['ivar_md'] = args.ivar_md
        column_names.insert(-1, 'ivar_md')

    # Write output header & QC columns to a CSV file
    with open(args.outfile, 'w') as csvfile:
        header = column_names
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        writer.writerow(qc_values)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outfile', required=True, type=str,
        help='''The path of the output QC summary file''')
    parser.add_argument('--sample', required=True, help='Sample name.')
    parser.add_argument('--ref', required=True, type=str,
        help='''The path of the reference FASTA file.''')
    parser.add_argument('--bam', required=True, type=str,
        help='''The path of the aligned & filtered BAM file.''')
    parser.add_argument('--fasta', required=True, type=str,
        help='''The path of a consensus fasta file produced by ivar using the
                minimum depth given by the --ivar_amd argument, required.''')
    parser.add_argument('--depths_file', required=True, type=str,
        help='''Tab-delimited file of lines each with 3 columns: contig name, 
                position, & number of reads aligned at this position.''')
    parser.add_argument('--flagstat_file', required=True, type=str,
        help='''Output file generated by the samtools flagstat command''')
    parser.add_argument('--minimum_depth', required=False, type=int, default=10,
        help='''Minimum depth value to be used for filtering when counting the 
                number of covered positions, optional.''')
    parser.add_argument('--ivar_md', required=False, type=int, default=None,
        help='''Minimum depth value used for ivar when generating the consensus
                file given by the --fasta argument, optional.''')

    args = parser.parse_args()
    generate_qc_file(args)

if __name__ == "__main__":
    main()
