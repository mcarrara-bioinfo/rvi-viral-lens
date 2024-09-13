# Workflow Documentation: `SORT_READS_BY_REF.nf`

## Overview

The `SORT_READS_BY_REF` workflow processes paired-end sequencing reads by sorting them according to taxonomic classifications obtained from Kraken2. This workflow uses a manifest file to process multiple samples and produces sorted by taxid FASTQ files for each sample and classification reports.

## Key Processes

1. **Run Kraken**: Classifies reads using Kraken2 against a specified database.
2. **Sort Reads**: Uses Kraken2Ref to sort reads by taxonomic ID and extracts them into separate FASTQ files.
3. **Merge FASTQ Parts**: Merges split FASTQ parts if necessary.
4. **Collect Reference Files**: Retrieves reference sequences based on taxonomic IDs for downstream analysis.

## Workflow Inputs and Outputs

### Inputs

- **Manifest File (`mnf_path`)**: A CSV file containing sample metadata and paths to paired-end FASTQ files.
- **Kraken Database Path (`params.db_path`)**: Path to the Kraken database.
- **Kraken2Ref Library Fasta (`params.db_library_fa_path`)**: Optional. Path to the Kraken2Ref library fasta file. If none is provided, it assumes there is a `${params.db_path}/library/library.fna`.

### Outputs

- `sample_taxid_ch`: Channel containing tuples of metadata and sorted reads per taxonomic ID.
- `sample_pre_report_ch`: Channel containing pre-reports with sample-level summaries.

## Workflow Steps

1. **Parse Manifest File**
The `parse_mnf` function reads the manifest file and extracts sample metadata and paths to paired-end FASTQ files. It outputs a channel with tuples of metadata and FASTQ pairs.

2. **Filter Out Samples with Empty FASTQ Files**
Samples with FASTQ files smaller than 5000 bytes are considered empty and are logged as warnings. The remaining samples are processed further.

3. **Run Kraken**
The `run_kraken` process classifies reads using the Kraken database specified in `params.db_path`. The output includes:

- Kraken output files
- Classified and unclassified FASTQ file pairs
- Kraken reports

4. **Sort Reads with Kraken2Ref**
The workflow sorts reads according to taxonomic classification:

- **Sort Reads**: The `run_k2r_sort_reads` process sorts reads by taxonomic IDs based on Kraken outputs.
- **Prepare Channels for Dumping FASTQs**: Determines which samples need splitting based on the number of reads.
- **Dump FASTQs and Generate Pre-Reports**: The `run_k2r_dump_fastqs_and_pre_report` process dumps sorted reads into FASTQ files and generates preliminary reports.
- **Merge FASTQ Parts**: The `concatenate_fqs_parts` process merges split FASTQ parts where applicable.

5. **Collect Reference Files**
The workflow retrieves reference files for each unique taxonomic ID using the `get_taxid_reference_files` process. These reference files are then combined with the corresponding sample data.

## Function Definitions

### `parse_mnf(consensus_mnf)`

- **Purpose**: Parses the manifest file to create a channel of metadata and FASTQ file pairs.
- **Input**: consensus_mnf (path to the manifest file)
- **Output**: Channel with tuples of metadata and FASTQ file pairs.

### `check_sort_reads_params()`

- **Purpose**: Checks for necessary parameters and validates paths to ensure they exist. Logs errors if any required parameters are missing.
- **Output**: Number of errors encountered during the checks.

### Parameter Checking

The function `check_sort_reads_params` ensures that required parameters (`db_path` and `manifest`) are provided and valid. Errors are logged for missing or invalid inputs.

### Example Usage

To run the workflow, provide the required parameters in a Nextflow run command:

```bash
nextflow run sort_reads_by_ref.nf \
    --manifest <path_to_manifest.csv> \
    --db_path <path_to_kraken_db> \
    --db_library_fa_path <path_to_library_fasta> \
    --k2r_max_total_reads_per_fq <max_reads_per_fastq>
```

Replace `<path_to_manifest.csv>`, `<path_to_kraken_db>`, `<path_to_library_fasta>`, and `<max_reads_per_fastq>` with appropriate paths and values for your analysis.

### Notes

- Ensure that the Kraken database is correctly formatted and accessible.
- Check the manifest file for correct formatting and paths to input files.
- The workflow can handle large datasets by splitting and merging FASTQ files based on the maximum reads threshold (`params.k2r_max_total_reads_per_fq`).