# Workflow Documentation: `GENERATE_CLASSIFICATION_REPORT.nf`

## Overview

The `GENERATE_CLASSIFICATION_REPORT` workflow generates a classification report based on metadata associated with sequencing samples. This workflow collects metadata from each sample, formats the data into a report line, and aggregates these lines into a final classification report file.

## Key Processes:

- **Report Line Generation**: Each sample's metadata is processed to generate a line of text summarizing the classification data.
- **Report Aggregation**: All report lines are aggregated into a single report file.

## Workflow Inputs and Outputs

### Inputs

- **Metadata Channel (`meta_ch`)**: A channel containing metadata for each sample. The metadata must include the following keys:
  - `sample_id`: Unique identifier for the sample.
  - `taxid`: Taxonomic ID of the sample.
  - `ref_selected`: Reference sequence used for the sample, formatted with commas but replaced by pipes (`|`).
  - `virus_name`: Name of the virus detected.
  - `virus_subtype`: Subtype of the virus, if applicable (defaults to `'None'` if `null`).
  - `flu_segment`: Influenza segment information, if applicable (defaults to `'None'` if `null`).
  - `percentage_genome_coverage`: Percentage of the genome covered by aligned reads.
  - `total_mapped_reads`: Number of reads that mapped successfully.
  - `longest_no_N_segment`: Length of the longest segment without ambiguous bases ('N').
  - `percentage_of_N_bases`: Proportion of ambiguous bases in the consensus sequence.

### Outputs

- **Classification Report (`write_classification_report.out`)**: A text file containing aggregated classification data for all samples, formatted as a CSV.

## Workflow Steps

1. **Generate Report Lines**
The workflow begins by generating a report line for each sample:

- Each sample's metadata is processed to create a formatted string, representing a line in the report.
- `Null` values for `virus_subtype` and `flu_segment` are replaced with `'None'` to ensure consistency in the output.
- Each line is formatted as a CSV string with fields separated by commas:

```
sample_id,virus,report_name,virus_name,taxid,ref_selected,flu_segment,virus_subtype,sample_subtype,percentage_genome_coverage,total_mapped_reads,longest_no_N_segment,percentage_of_N_bases
```

The generated lines are collected into a single channel (`report_lines_ch`).

2. **Write Classification Report**

The `write_classification_report` process takes the aggregated report lines and writes them to a classification report file:

- All per-sample report lines are concatenated into a single text file.
- The final report is outputted through the `write_classification_report.out` channel.

3. **Emit Output**
The workflow emits the final classification report file through the output channel (`write_classification_report.out`).

## Function Definitions

### Metadata Structure Requirements

The metadata (`meta`) passed into the workflow must include the following keys all keys listed Inputs section of this documentation

### Example Manifest File

The manifest file should be in CSV format and include headers matching the required metadata keys. Example:

```csv
sample_id,taxid,ref_selected,virus_name,virus_subtype,flu_segment,percentage_genome_coverage,total_mapped_reads,longest_no_N_segment,percentage_of_N_bases
sample_001,12345,ref1.fasta,Influenza A,H1N1,Segment 1,95.5,10000,500,2.5
sample_002,67890,ref2.fasta,Influenza B,,Segment 2,90.0,9500,450,5.0
```

Example Usage

To run the workflow, provide the required manifest file via Nextflow run command:

```bash
nextflow run generate_classification_report.nf --manifest_file <path_to_manifest.csv>
Replace <path_to_manifest.csv> with the path to your manifest file.
```

## Notes

- The running the workflow is useful for debuging it in isolation, but has no usage outside of this pipeline context.
- Ensure that all required metadata fields are correctly formatted and available in the input channel.
- The workflow currently handles null values for virus_subtype and flu_segment by substituting 'None', ensuring that these fields are not left empty in the final report.
- Review the output of the QC script to verify that the values are parsed correctly, especially if the script's output format changes.