# QC Script for BAM and FASTA Files

## Overview

This Python script generates a quality control (QC) summary report for a sample by analyzing a BAM file, a consensus FASTA file, a reference FASTA, and a per-position depth file. The QC metrics include the percentage of N bases in the consensus, the largest contiguous gap of N bases, and the percentage of reference bases covered at a minimum depth. The script also integrates alignment statistics from a SAMtools `flagstat` output to provide additional insights on the quality of the aligned reads.

## Inputs

1. **Consensus FASTA File** (`--fasta`): The consensus sequence file.

2. **Reference FASTA File** (`--ref`): The reference genome against which the reads were aligned and used for depth calculations.

3. **BAM File** (`--bam`): The aligned and filtered BAM file, which contains the reads aligned to the reference.

4. **Per-position Depths File** (`--depths_file`): A tab-delimited file listing read depth per position in the alignment.

5. **SAMtools Flagstat File** (`--flagstat_file`): The output from the `samtools flagstat` command, providing alignment statistics.

6. **Sample Name** (`--sample`): The name of the sample being processed.

7. **Output File** (`--outfile`): The path where the QC summary report (in CSV format) will be written.

8. **Minimum Depth** (`--minimum_depth`, optional): The minimum depth threshold used when calculating covered bases. Default is `10`.

9. **Ivar Minimum Depth** (`--ivar_md`, optional): The minimum depth used by ivar when generating the consensus FASTA.

## Outputs

The script generates a CSV file containing various QC metrics for the sample. The output columns include:

- `sample_name`: Name of the sample.
- `pct_N_bases`: Percentage of N bases in the consensus FASTA.
- `pct_covered_bases`: Percentage of the reference genome covered at or above the minimum depth threshold.
- `longest_no_N_run`: Length of the largest contiguous region without N bases in the consensus.
- `num_aligned_reads`: Number of aligned reads (from `SAMtools flagstat`).
- `bam`: Path to the input BAM file.
- `total_mapped_reads`: Total number of mapped reads (from `SAMtools flagstat`).
- `total_unmapped_reads`: Total number of unmapped reads (from `SAMtools flagstat`).
- `qc_pass`: Indicates whether the sample passed QC based on N content and gap size criteria.
- `ivar_md`: If applicable, the minimum depth used by ivar.

## Key Features

- **N Content Analysis**: The script calculates the percentage of bases that are 'N' in the consensus FASTA and finds the longest contiguous stretch of N bases.

- **Depth-Based Coverage**: It calculates the percentage of reference genome positions covered by reads at or above a specified depth threshold, using the per-position depths file.

- **Alignment Statistics**: The script reads the `SAMtools flagstat` file to report the number of aligned, mapped, and unmapped reads.

- **QC Pass/Fail**: The QC status is determined by evaluating N content (e.g., `TRUE` if the largest N gap exceeds `10,000` or if `N` bases constitute less than 50% of the sequence).

## Workflow

1. **Input Parsing**: The script takes input arguments using argparse for the BAM, FASTA, reference, and depth files, as well as additional parameters like the minimum depth threshold and output file path.

2. **QC Metric Calculation**:
   - It calculates the percentage of N bases and identifies the longest N gap in the consensus sequence.
   - It counts the positions in the depth file where the depth exceeds or equals the minimum threshold and calculates the coverage percentage.

3. **Flagstat Integration**: The script reads alignment statistics from the SAMtools flagstat file.

4. **QC Report Generation**: It assembles the calculated QC metrics into a dictionary, formats them, and writes them as a CSV report.

## Example Command

```bash
./generate_qc.py --outfile sample123.qc.csv --sample sample123 --ref ref.fasta --bam sample123.bam --fasta sample123.consensus.fasta --depths_file sample123.depths.txt --flagstat_file sample123.flagstat.txt --minimum_depth 10 --ivar_md 5
```

This command processes the provided files and generates a QC summary report in `sample123.qc.csv`. The minimum depth for coverage calculations is set to 10, and ivar was run with a minimum depth of 5.
