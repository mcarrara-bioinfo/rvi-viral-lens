# Nextflow Process Documentation: `run_ivar`

## Overview

The `run_ivar` process in this Nextflow pipeline is designed to reproduce the ARTIC pipeline approach for generating consensus sequences from BAM files using the samtools and ivar tools. Its main objective is create consensus sequences based on specified depth and frequency thresholds. The resulting files are organized by sample and taxonomic ID in the specified results directory.

## Parameters

- `params.ivar_min_depth`: The minimum read depth required to make a base call in the consensus sequence. If the depth at a given position is below this threshold, an `'N'` will be used instead. The default value is `10`.

- `params.ivar_freq_threshold`: Minimum frequency threshold (0 - 1) to call consensus. Bases with a frequency below this threshold are not called in the consensus. The default value is 0.75.

## Process Description

The `run_ivar` process generates a consensus sequence from a sorted BAM file using samtools mpileup and ivar consensus. Below is a breakdown of the key components of this process:

## Tags and Labels

- Tag: `meta.id` – This tag is used for tracking and logging purposes within the Nextflow pipeline, allowing identification of process runs by sample ID.

- Label: `"ivar"` – This label is used for resource configuration in Nextflow, typically to specify compute resources like memory or CPU for the ivar process.

## Input

- Input Tuple: tuple val(meta), path(bams)
  - `meta`: Metadata associated with the sample, containing identifiers such as sample ID and taxonomic ID.
  - `bams`: Path to the BAM file(s) to be processed.

## Output

- Output Tuple: `tuple val(meta), path("${meta.id}.consensus.fa")`
  - The output consists of the metadata and the path to the generated consensus FASTA file (.fa) for the given sample.

## Publish Directory

- Publish Directory: `${params.results_dir}/\${meta.sample_id}/\${meta.taxid}/`
  - The results are published to a directory structure based on sample ID and taxonomic ID, copying files with extensions .fa, .tsv, and .txt.

## Script

The script section includes the commands executed within the process:

```bash
set -e
set -o pipefail

samtools mpileup -aa -A -B -d 0 -Q0 ${sorted_bam} | \
  ivar consensus -t ${params.ivar_freq_threshold} -m ${params.ivar_min_depth} -n N -p ${meta.id}.consensus
```

## Command Breakdown

- `set -e`: This command ensures that the script exits immediately if a command exits with a non-zero status, providing error handling during execution.

- `set -o pipefail`: This option causes the pipeline to return the exit status of the last command in the pipeline that failed, which helps in debugging errors.
- `samtools mpileup`: This command generates a pileup of reads from the BAM file. The options used are:
  - `-aa`: Output absolutely all positions, including unused reference sequences and zero depth.
  - `-A`: Do not skip anomalous read pairs in variant calling. 
    - equiavelent to `--count-orphans`
    - Anomalous read pairs are those marked in the FLAG field as paired in sequencing but without the properly-paired flag set.
  - `-B`: Disable BAQ ([Base Alignment Quality](https://academic.oup.com/bioinformatics/article/27/8/1157/227268)).
  - `-d 0`: Set the maximum depth per file to unlimited. 
    - At a position, read maximally INT reads per input file. Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. Passing zero for this option sets it to the highest possible value, effectively removing the depth limit. [8000]
  - `-Q0`: Set the minimum base quality to 0.
    - equivalent to `--min-BQ`.
    - Minimum base quality for a base to be considered. Note base-quality 0 is used as a filtering mechanism for overlap removal which marks bases as having quality zero and lets the base quality filter remove them. Hence using `--min-BQ 0` will make the overlapping bases reappear, albeit with quality zero.

- `ivar consensus`: This command creates a consensus sequence from the pileup data:
  - `-t ${params.ivar_freq_threshold}`: Sets the threshold for base frequency.
  - `-m ${params.ivar_min_depth}`: Minimum depth to call consensus.
  - `-n N`: Set `N` as the character to print in regions with less than minimum coverage.
  - `-p ${meta.id}.consensus`: Specifies the prefix for the output consensus sequence file.

Information about parameters were obtained from [samtools mpileup](http://www.htslib.org/doc/samtools-mpileup.html) and [ivar documentation](https://andersen-lab.github.io/ivar/html/manualpage.html#**autotoc_md19)
