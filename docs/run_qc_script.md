# Nextflow Process Documentation: `run_qc_script`

## Overview

This process iperforms quality control (QC) on sequence data using a custom Python QC script (available at `bin/` dir). This process generates a CSV file with QC metrics and a PNG image visualizing sequencing depth, which are essential for assessing the quality and reliability of the sequencing data. The script is a modified version of [the `qc.py` of the ARTIC pipeline approach](https://gitlab.internal.sanger.ac.uk/malariagen1/ncov2019-artic-nf/-/blob/main/bin/qc.py?ref_type=heads).

## Process

This process runs a QC analysis on the input BAM, FASTA, and reference files, outputting QC metrics and a sequencing depth plot. Below is a detailed breakdown of the process components:

### Tags and Labels

- Tag: `${meta.id}` – The tag identifies the process run by the sample ID, facilitating tracking and logging within the pipeline.

### Input

- Input Tuple: `tuple val(meta), path(bam), path(fasta), path(ref)`
  - `meta`: Metadata associated with the sample, assumes the following keys  are available: 
    - `id`: internal pipeline ID
    - `sample_id`: sample ID
    - `taxid`: taxonomic ID
  - `bam`: Path to the BAM file containing aligned sequencing reads.
  - `fasta`: Path to the consensus FASTA file for the sample.
  - `ref`: Path to the reference file against which the reads are aligned.

### Output

- Output Tuple: `tuple val(meta), path("${meta.id}.qc.csv"), path("${meta.id}.depth.png"), stdout`
  - A CSV file containing QC metrics (`${meta.id}.qc.csv`).
  - A PNG file visualizing the depth of sequencing (`${meta.id}.depth.png`).
  - Standard output (`stdout`) prints the first row of the QC CSV file for quick verification.

### Publish Directory

- Publish Directory: `${params.results_dir}/${meta.sample_id}/${meta.taxid}/`
  - The results are published in a structured directory organized by sample ID and taxonomic ID, making it easy to locate and review QC results.

### Script

```bash
qcSetting="--illumina"

qc.py ${qcSetting} \
    --outfile ${meta.id}.qc.csv \
    --sample ${meta.id} \
    --ref ${ref} \
    --bam ${bam} \
    --fasta ${fasta} \
    --ivar_md ${params.ivar_min_depth}

# Print the first row of the QC CSV for quick verification
sed -n "2p" ${meta.id}.qc.csv
```

### Script Breakdown

1. **QC Script Execution**:

- Command: `qc.py ${qcSetting} --outfile ${meta.id}.qc.csv --sample ${meta.id} --ref ${ref} --bam ${bam} --fasta ${fasta} --ivar_md ${params.ivar_min_depth}`

  - The command runs a custom Python script (`qc.py`) for QC analysis, with the following parameters:
    - `--illumina`: A setting flag for the QC script, suggesting it is configured for Illumina sequencing data.
    - `--outfile`: Specifies the output file name for QC metrics in CSV format.
    - `--sample:` Specifies the sample ID for the analysis.
    - `--ref`: Provides the path to the reference file.
    - `--bam`: Supplies the BAM file with aligned reads.
    - `--fasta`: Provides the path to the consensus FASTA file for the sample.
    - `--ivar_md`: Sets the minimum depth threshold for ivar, passed as a parameter (`${params.ivar_min_depth}`).

For more details, check [bins DOCs [TODO])](TODO)

1. **Extracting QC Results**:

- Command: `sed -n "2p" ${meta.id}.qc.csv`
  - This command prints the second line of the QC CSV file, the first row of data (excluding headers). It provides a quick check of the QC output for monitoring and verification purposes.
