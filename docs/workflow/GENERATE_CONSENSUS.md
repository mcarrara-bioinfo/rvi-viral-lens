# Workflow Documentation: `GENERATE_CONSENSUS`

## Overview

The `GENERATE_CONSENSUS` workflow performs read alignment and consensus sequence generation for sequencing data. It processes paired-end reads by aligning them to reference genomes using BWA, followed by consensus calling with iVar. This workflow is designed to take in sequencing data for different samples and taxonomic IDs, process them, and produce consensus sequences.

## Key Processes

- **BWA Alignment**: Aligns sequencing reads to the provided reference genome.
- **Consensus Calling with iVar**: Generates consensus sequences from the aligned reads.

## Workflow Inputs and Outputs

### Inputs

- **Sample Taxid Channel (`sample_taxid_ch`)**: A channel containing tuples of metadata and paired-end FASTQ files. Metadata (`meta`) must include the following keys:
  - `id`: Unique identifier combining sample ID and taxonomic ID.
  - `taxid`: Taxonomic ID of the sample.
  - `sample_id`: Sample identifier.
  - `ref_files`: Paths to reference genome files.

### Outputs

- `run_ivar.out`: A channel containing tuples of metadata and the generated consensus FASTA file.

## Workflow Steps

1. **Prepare Input for BWA Alignment**
The workflow starts by preparing the input channel for BWA alignment:

- It maps metadata and reads from the `sample_taxid_ch` channel.
- The channel output format is a tuple containing metadata, reads, and reference genome paths (`bwa_input_ch`).

> **DEV NOTE**: we get the reference files from `meta`, we should change that to have the files explicitly stated on the input channels. Make sense to get it from meta, if we start from a consensus manifest, but we should do that only under that condition.

1. **BWA Alignment and Post-Processing**

The `bwa_alignment_and_post_processing` process aligns reads to the reference genomes using BWA. It outputs:

- A sorted BAM file and its corresponding index file (`.bai`).

3. **Prepare Input for iVar Consensus Generation**
The output from BWA alignment is further processed:

- The BAM file is added to the metadata.
- The channel is prepared for input into the iVar process (`ivar_in_ch`).

> **DEV NOTE**: maybe we should consider to put `.bam` file explicitly in the emitted channel.

4. Generate Consensus with iVar
The `run_ivar` process takes the input channel (`ivar_in_ch`) containing metadata and aligned BAM files:

- It generates consensus sequences in FASTA format based on the BAM alignments.

5. **Emit Output**
The consensus FASTA files are emitted from the workflow via the output channel (`run_ivar.out`).

## Function Definitions

### `parse_consensus_mnf_meta(consensus_mnf)`

- **Purpose**: Parses a manifest file to create a channel of metadata and paired-end FASTQ files for each sample.
- **Input**: `consensus_mnf` (path to the manifest file)
- **Output**: A channel containing tuples of metadata and FASTQ file pairs.

**Example Manifest File Format:**

The manifest file should be a CSV with the following columns:

- `sample_id`: Sample identifier.
- `taxid`: Taxonomic ID of the sample.
- `ref_files`: Semicolon-separated list of reference genome paths.
- `reads_1`: Path to the first FASTQ file (paired-end).
- `reads_2`: Path to the second FASTQ file (paired-end).

### `check_generate_consensus_params()`

- **Purpose**: Placeholder function for parameter checking. Currently, it does not perform any checks.
- **Output**: Returns the number of errors found (currently always 0).

> **DEV NOTE**: Maybe is time to add that functionality

### Example Usage

To run the workflow, provide the required parameters in a Nextflow run command:

```bash
nextflow run generate_consensus.nf \
    --manifest <path_to_manifest.csv> \
    --db_path <path_to_kraken_db> \
    --db_library_fa_path <path_to_library_fasta>
```

Replace `<path_to_manifest.csv>`, `<path_to_kraken_db>`, and `<path_to_library_fasta>` with the appropriate paths for your analysis.

### Notes

- Ensure that reference files are correctly specified and accessible in the manifest file.
- The BWA process requires indexed reference genomes; ensure that reference files are pre-processed as needed for BWA.
- The iVar process will output consensus sequences in FASTA format; ensure that the output path is properly set if required.
