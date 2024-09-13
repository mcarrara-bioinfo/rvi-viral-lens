# Workflow Documentation: `COMPUTE_QC_METRICS`

## Overview

The `COMPUTE_QC_METRICS` workflow is designed to compute quality control (QC) metrics for consensus sequences generated from sequencing data. The workflow processes each sample's data to evaluate the quality and coverage of the generated consensus sequences. The QC metrics include the percentage of bases covered, the percentage of N bases, the longest segment without N bases, and read alignment statistics (total reads aligned, unmapped and mapped).

## Key Processes

- **QC Metrics Calculation**: Uses a custom QC script (check QC script documentation) to compute various quality metrics based on the aligned BAM file and the consensus FASTA sequence.

## Workflow Inputs and Outputs

### Inputs

- **QC Metrics Input Channel (`qc_metrics_In_ch`)**: A channel containing tuples of metadata and the consensus FASTA file.
  - Metadata (`meta`) must include the following keys:
    - `id`: Unique identifier combining sample ID and taxonomic ID.
    - `taxid`: Taxonomic ID of the sample.
    - `sample_id`: Sample identifier.
    - `ref_files`: Paths to reference genome files.

### Outputs

- `qc_Out_ch`: A channel containing tuples with updated metadata that includes computed QC metrics.

## Workflow Steps

1. **Prepare Input for QC Script**
The workflow starts by preparing the input channel for the QC script:

- It maps metadata and the consensus FASTA file from the `qc_metrics_In_ch` channel.
- The BAM file and reference genome path are extracted from the metadata and used as inputs to the QC script.
- The output format is a tuple containing metadata, BAM file, consensus FASTA file, and reference genome path (`qc_script_In_ch`).

2. **Compute QC Metrics**
The run_qc_script process computes QC metrics using the prepared input channel (`qc_script_In_ch`):
- QC metrics calculated include:
  - **Percentage of N bases**: Proportion of ambiguous bases in the consensus sequence.
  - **Percentage Genome Coverage**: Proportion of the genome covered by aligned reads.
  - **Longest No-N Segment**: The longest continuous segment without N bases.
  - **Total Aligned Reads**: Number of reads aligned to the reference.
  - **Total Unmapped Reads**: Number of reads that did not align.
  - **Total Mapped Reads**: Number of reads that mapped successfully.

3. **Populate Metadata with QC Values**
The output of the QC script is used to populate the metadata with the computed QC metrics:
- The script output is tokenized and values are extracted to update the metadata.
- The channel is updated to include only metadata and the BAM file (`qc_Out_ch`).

4. **Emit Output**
The workflow emits the final output channel (`qc_Out_ch`) containing updated metadata with the QC metrics.

## Function Definitions

### Metadata Structure Requirements

The metadata (`meta`) passed through the workflow must include:

- `id`: A unique identifier that combines sample ID and taxonomic ID.
- `taxid`: The taxonomic ID of the sample.
- `sample_id`: The identifier for the sample.
- `ref_files`: Paths to the reference genome files associated with the sample.

### Example QC Metrics Script Output Parsing

The `run_qc_script` output is parsed as follows:

- The script output is a comma-separated string.
- Values are extracted based on their positions in the string:
  - `percentage_of_N_bases`: 2nd token.
  - `percentage_genome_coverage`: 3rd token.
  - `longest_no_N_segment`: 4th token.
  - `total_aligned_reads`: 5th token.
  - `total_unmapped_reads`: 11th token.
  - `total_mapped_reads`: 9th token.

## Notes

The QC script's output is expected to be in a specific format; any changes in this format may require adjustments in the parsing logic within the workflow.