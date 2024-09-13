# Nextflow Process Documentation: `run_kraken2ref_pre_report.nf`

## Process: `run_k2r_sort_reads`

### Overview

This process runs Kraken2Ref to parse Kraken reports and sort reads by taxonomic classification. It generates JSON files that map taxonomic IDs to read IDs and performs sorting based on the decomposed taxonomy tree if available.

### Tags and Labels

- Tag: `${meta.id} - ${task.attempt} - ${task.memory}` - Identifies the process run by sample ID, number of attempts, and memory requested.

- Labels: `kraken2ref`, `mem_k2r_escalate` - Defines resource requirements for the process, namely which container to use and escalating memory needs for Kraken2Ref operations.

check [ADD CONF DOCS]() for more details on the memory escalation strategy

### Input

- Input Tuple: `tuple val(meta), path(kraken_output), path(kraken_report)`
  - `meta`: Metadata including identifiers like sample ID.
  - `kraken_output`: Path to the Kraken output file containing classification results.
  - `kraken_report`: Path to the Kraken report file summarizing classifications.

### Output

- Output Tuple: `tuple val(meta), path("${meta.id}_tax_to_reads.json"), path("${meta.id}_decomposed.json"), optional: true, emit: json_files`
  - JSON files mapping taxonomic IDs to read IDs (`${meta.id}_tax_to_reads.json`)
  - decomposed JSON file (`${meta.id}_decomposed.json`). The decomposed JSON file is optional and only generated if applicable.

### Publish Directory

- Publish Directory: `${params.results_dir}/${meta.id}/reads_by_taxon/`

Outputs are published in the results directory organized by sample ID, ensuring files are accessible and systematically stored.

### Script

The script section runs Kraken2Ref commands to parse and sort reads based on Kraken reports:

```bash
kraken2ref -s ${meta.id} parse_report -i ${kraken_report} -o ./ -t ${params.min_reads_for_taxid}

# If the decomposed JSON file is generated, proceed to sort reads
if [ -e "${meta.id}_decomposed.json" ]; then
    kraken2ref -s ${meta.id} sort_reads -k ${kraken_output} -r ./${meta.id}_decomposed.json -m tree -u
else
    echo "Warning: JSON file does not exist. Skipping downstream commands."
fi
```

### Script Breakdown

1. **Parsing Kraken Report**: `kraken2ref -s ${meta.id} parse_report -i ${kraken_report} -o ./ -t ${params.min_reads_for_taxid}`
Parses the Kraken report to generate JSON mapping of taxonomic IDs to read IDs. The `-t` option sets the minimum number of reads for a taxonomic ID to be included.

2. **Sorting Reads**: `kraken2ref -s ${meta.id} sort_reads -k ${kraken_output} -r ./${meta.id}_decomposed.json -m tree -u`
Sorts reads using the decomposed JSON if available, based on a taxonomy tree structure.
If the decomposed JSON file is missing, sorting is skipped with a warning message.

> **DEV NOTE**: we check for decomposed json file, at this point there should not be empty fastq files, so I think we should let the pipeline brake if this is the case.

---

## Process: `run_k2r_dump_fastqs_and_pre_report`

### Overview

The `run_k2r_dump_fastqs_and_pre_report` process uses Kraken2Ref to extract classified reads into FASTQ files and generate a preliminary report based on the taxonomic classification data. It processes classified reads and produces a detailed report for further analysis. The command kraken2ref `-s ${prefix}` dump_fastqs runs the `dump_fastqs` function of Kraken2Ref, which extracts reads that match specific taxonomic IDs into FASTQ files.

### Tags and Labels

- Tag: `${meta.id} - ${task.attempt} - ${task.memory}`
- Labels: `kraken2ref`, `mem_k2r_escalate`
- Memory: Dynamically set based on the `params.k2r_dump_fq_mem` [default 6GB] value to handle potentially large datasets.

### Input

- Input Tuple: `tuple val(meta), path(classified_fqs), path(json_tax_to_readsid), path(decomposed_json), path(kraken_report)`
  - `classified_fqs`: Paths to paired FASTQ files containing classified reads.
  - `json_tax_to_readsid`: JSON file mapping taxonomic IDs to read IDs.
  - `decomposed_json`: Decomposed JSON file representing detailed taxonomic hierarchy.
  - `kraken_report`: Kraken report summarizing classification data.

### Output

- Output Tuples:
  - `tuple val(meta), path("*_R{1,2}.fq"), optional: true, emit: fq_files`
    - FASTQ files with extracted reads by taxonomic classification.
  - `path("${meta.id}_pre_report.tsv"), optional: true, emit: report_file`
    - Preliminary TSV report file containing classification summary data.

### Script

The script runs Kraken2Ref to extract classified reads into FASTQ files and generate a preliminary report:

```bash
if [ "!{meta.splitted}" == "true" ]; then
    part=$(echo !{fq_1}| awk -F'[.]' '{print $(NF-1)}')
    prefix="${part}-!{meta.id}"
else
    prefix="!{meta.id}"
fi

kraken2ref -s ${prefix} dump_fastqs \
        -fq1 !{fq_1} -fq2 !{fq_2} \
        --tax_to_readsid_path !{json_tax_to_readsid} \
        -o ./ --fq_load_mode !{params.k2r_fq_load_mode} \
        -r !{decomposed_json}
```

### Script Breakdown

This script is designed to handle the processing of classified FASTQ files using Kraken2Ref. It dynamically sets a prefix for the output files based on whether the input FASTQ files are split parts and then runs Kraken2Ref to extract reads into separate FASTQ files for further analysis.

1. **Determine Prefix Based on File Splitting**
Purpose: Sets the prefix for output files based on whether the FASTQ splitted files are splitted or not.
   1. if `meta.splitted` variable is set to "true", indicating that the FASTQ files are split into parts. then, it extracts the part identifier from the filename of `fq_1` and the prefix is set as `${part}-!{meta.id}`.
   2. If not split, the prefix is set as the ID: `!{meta.id}`.

2. **Run Kraken2Ref to Extract FASTQs**
Purpose: Runs Kraken2Ref to extract reads based on taxonomic classification into separate FASTQ files.
Parameters:
    - `-s ${prefix}`: Sets the prefix for the output files, determined in the previous step.
    - `-fq1 !{fq_1} -fq2 !{fq_2}`: Specifies the input paired FASTQ files.
    - `--tax_to_readsid_path !{json_tax_to_readsid}`: Provides the path to the JSON file that maps taxonomic IDs to read IDs, guiding which reads to extract.
    - `-o ./`: Sets the output directory to the current directory (`./`).
    - `--fq_load_mode !{params.k2r_fq_load_mode}`: Specifies the mode for loading FASTQ files, defined by the parameter `params.k2r_fq_load_mode`.
    - `-r !{decomposed_json}`: Uses the decomposed JSON file to provide detailed taxonomic structure for sorting reads.

3. **Generate preliminary report**
Purpose: Produce a preliminary report in TSV format summarizing the classification data.

---

## Process: `concatenate_fqs_parts`

### Overview

This process concatenates FASTQ files from multiple parts into final combined FASTQ files for each taxonomic classification. This process ensures that all parts corresponding to the same taxonomic ID are merged into single files.

### Tags and Labels

- Tag: Not explicitly defined.
- Labels: Not explicitly defined.

### Input

- Input Tuple: `tuple val(id), path(fq_parts)`
  - `id`: Identifier for the sample or taxonomic classification.
  - `fq_parts`: Paths to FASTQ file parts that need to be concatenated.

### Output

- Output Tuple: `tuple val(id), path("*_R{1,2}.fq")`
    -Concatenated FASTQ files for each taxonomic classification.

### Script

The script concatenates all matching FASTQ parts into final output files:

```bash
#!/bin/bash

# Loop through all FASTQ files in the current directory
for file in *_R[12].fq; do
    # Extract taxid from the file name
    taxid=$(echo "$file" | awk -F'[_]' '{print $(NF-1)}')

    # Define the output filenames
    fq_1_out="!{id}_${taxid}_R1.fq"
    fq_2_out="!{id}_${taxid}_R2.fq"

    # Check if the output files already exist
    if [ -f "$fq_1_out" ]; then
        echo "Concatenation for !{id}-${taxid} already done"
    else
        # Concatenate all matching files into the output files
        echo "Assembling $fq_1_out and $fq_2_out"
        cat *-!{id}_${taxid}_R1.fq >> "$fq_1_out"
        cat *-!{id}_${taxid}_R2.fq >> "$fq_2_out"
    fi
done
```

### Script Breakdown

1. **File Loop**:
Iterates over FASTQ files with suffix `_R1.fq` and `_R2.fq`, assuming paired-end sequencing data.

2. **Output Filename Construction**:
Constructs output filenames based on sample or taxonomic ID and checks if concatenation has already been completed to avoid redundancy.

3. **Concatenation**:
Concatenates parts into the designated output files using cat, appending parts that match the specified naming convention for each taxonomic classification.
