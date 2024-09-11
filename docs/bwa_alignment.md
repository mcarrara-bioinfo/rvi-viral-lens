# Nextflow Process Documentation: `bwa_alignment_and_post_processing`

## Overview

The process is responsible for mapping sequencing reads to a reference genome using BWA (Burrows-Wheeler Aligner), followed by post-processing steps including conversion to BAM format, sorting, and indexing. At current version of this pipeline, this process generates high-quality, sorted BAM files that are essential for consensus sequence generation.

## Process

This process aligns reads to a reference genome and outputs indexed BAM files. Below is a detailed breakdown of the process components:

### Tags and Labels

- Tag: Not explicitly tagged in this process, but meta information is utilized for sample and taxonomic identifiers.

- Label: No label is applied to this process, therefore, it is running on base container [Labels documentation[TODO]]().

### Input

- Input Tuple: `tuple val(meta), path(fastq), path(ref_files)`
  - `meta`: Metadata associated with the sample. The process assumes there is:
    - `taxid` : taxonomic ID
    - `sample_id`: sample ID
    - 
  - `fastq`: Path to the FASTQ file(s) containing sequencing reads to be aligned.
  - `ref_files`: Path to the reference files required for alignment, including the reference FASTA (`${meta.taxid}.fa`).

### Output

- Output Tuple: `tuple val(meta), path("*.sorted.bam*")`
  - Sorted BAM file (`*.sorted.bam`).
  - Index file for the sorted BAM (`*.sorted.bam.bai`).

### Publish Directory

- Publish Directory: `${params.results_dir}/${meta.sample_id}/${meta.taxid}/`

The results are published to a directory organized by sample ID and taxonomic ID, facilitating the organized storage and easy retrieval of alignment results.

### Script

The script section includes the shell commands executed within the process:

```bash
ref_fa="${meta.taxid}.fa"

set -e
set -o pipefail

# Run BWA to align reads to the reference genome
bwa mem ${ref_fa} ${fastq} > ${meta.id}.sam

# Convert SAM to BAM format
samtools view -S -b ${meta.id}.sam -o ${meta.id}.bam

# Sort the BAM file by leftmost coordinates
samtools sort ${meta.id}.bam -o ${meta.id}.sorted.bam

# Generate BAM index
samtools index ${meta.id}.sorted.bam
```

### Script Breakdown

1. **Error Handling**:

- `set -e`: Exits the script if any command returns a non-zero status, ensuring that errors are caught early.
- `set -o pipefail`: Ensures that errors in piped commands are correctly propagated.

2. **BWA Alignment**:

- Command: `bwa mem ${ref_fa} ${fastq} > ${meta.id}.sam`
  - Aligns the input sequencing reads (`${fastq}`) to the specified reference genome (`${ref_fa}`), producing a SAM file (`${meta.id}.sam`) containing the raw alignments.

3. **Convert SAM to BAM**:

- Command: `samtools view -S -b ${meta.id}.sam -o ${meta.id}.bam`
  - Converts the SAM file to BAM format, which is a compressed binary version of the alignment file, making it more efficient for storage and processing.

4. **Sort BAM File**:

- Command: `samtools sort ${meta.id}.bam -o ${meta.id}.sorted.bam`
  - Sorts the BAM file by the leftmost coordinates of the alignments, a necessary step for many downstream analyses.

5. **Generate BAM Index**:

- Command: `samtools index ${meta.id}.sorted.bam`
  - Indexes the sorted BAM file, creating an index file (`*.sorted.bam.bai`) that allows for fast access to specific regions of the alignments.
