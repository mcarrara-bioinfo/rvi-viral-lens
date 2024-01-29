process compute_percentage_covered_bases {
    tag "${meta.id}"
    input:
        tuple val(meta), path(consensus_fa), path(ref_fasta)
    output:
        tuple val(meta), stdout

    script:
$/
#!/usr/bin/python3
from Bio import SeqIO

def read_fasta_file(file_path):
    #read fasta file account for split lines
    with open(file_path, "r") as file:
        sequence = "".join([str(record.seq) for record in SeqIO.parse(file, "fasta")])
    return sequence

# read the FASTA sequence
consensus_sequence = read_fasta_file("${consensus_fa}")
ref_sequence = read_fasta_file("${ref_fasta}")

# number of ns in the consensus sequence
num_n = consensus_sequence.count('N')

# length of the sequence
consensus_length_sequence = len(consensus_sequence)
reference_length_sequence = len(ref_sequence)

# number of nucleotides that are not n
num_non_n = consensus_length_sequence - num_n

# percentage of covered bases
pct_covered_bases = (num_non_n/ reference_length_sequence)*100

# stdout
print(f"{pct_covered_bases}|{consensus_length_sequence}|{num_n}")
/$
}

/*
---
def load_fasta_seq(fasta_file):
    # Read the FASTA file and extract the sequence
    sequence = ''
    with open(fasta_file, 'r') as file:
        for line in file:
            # NOTE: Assumes theres is only one sequence
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

# load consensus FASTA file
fasta_file_path = "${consensus_fa}"
consensus_seq = load_fasta_seq(fasta_file_path)
consensus_seq_size = len(consensus_seq)

# load reference fasta
ref_fasta_file_path = "${ref_fasta)}"
---

ref_seq_size = len(load_fasta_seq(ref_fasta_file_path))


# get number of valid calls for consensus 
#1 - Count the number of 'N' characters in the sequence
n_count = sum(1 for base in consensus_seq if base.upper() == 'N')

# 2 - Compute number of valid calls (everything non-N is considered a valid call)
non_n_count = consensus_seq_size - n_count

# 3 - Compute number of difference seq sizes between consensus and ref
delta_n_bases = ref_seq_size - consensus_seq_size

# 4 - Compute ref coverage
coverage = non_n_count / ref_seq_size

print(f"{consensus_seq_size}|{non_n_couts}|{delta_n_bases}|{coverage}")

# get number of valid calls for consensus 
#1 - Count the number of 'N' characters in the sequence
n_count = sum(1 for base in consensus_seq if base.upper() == 'N')

# 2 - Compute number of valid calls (everything non-N is considered a valid call)
non_n_count = consensus_seq_size - n_count

# 3 - Compute number of difference seq sizes between consensus and ref
delta_n_bases = ref_seq_size - consensus_seq_size

# 4 - Compute ref coverage
coverage = non_n_count / ref_seq_size

print(f"{consensus_seq_size}|{non_n_couts}|{delta_n_bases}|{coverage}")


    # Specify the path to your FASTA file
    fasta_file=${consensus_fa}
    ref_fasta=${ref_fasta}
    # Extract the sequence (excluding the header) from the FASTA file
    sequence=\$(tail -n +2 "\$fasta_file" | tr -d '\n')
    ref_sequence=\$(tail -n +2 "\$ref_fasta" | tr -d '\n')

    # Count the number of non-N characters in the sequence
    non_n_count=\$(echo "\$sequence" | tr -cd 'ACGTacgt' | wc -c)

    # Count the total number of characters in the sequence
    total_count=\$(echo "\$sequence" | wc -c)

    # Calculate the fraction of valid calls (non-N)
    fraction_valid_calls=\$(awk 'BEGIN {printf \$non_n_count / \$total_count}')

    # Print the result
    echo "Fraction of valid calls (non-N) in \$fasta_file: \$fraction_valid_calls"
*/