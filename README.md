[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.1-23aa62.svg)](https://www.nextflow.io/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) ![Nf-test](https://img.shields.io/badge/NFtest-%E2%89%A50.8.4-23aa62.svg?labelColor=0000)

![](./docs/assets/vira_pipeline_logo_placeholder.png)

The **VIRAL_PIPELINE** is a Nextflow pipeline developed under the context of the [RVI project]() by [GSU]() and its main goal is to identify the presence of Flu, SARS-CoV-2 and RSV and obtain, if possible, high quality consensus sequences for those virus. For more details, check the [[ADD REFERENCE PAPER]]()

> [THE CURRENT LOGO IS A **PLACEHOLDER** AND MUST BE UPDATED TO THE FINAL ONE BEFORE OPEN SOURCE THIS]

---
## Contents
- [Contents](#contents)
- [Pipeline Summary](#pipeline-summary)
- [How to Cite](#how-to-cite)
- [Quick Start](#quick-start)
- [Installation](#installation)
  - [dependencies](#dependencies)
  - [build containers](#build-containers)
- [Usage](#usage)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Configuration Sections](#configuration-sections)
  - [Parameters (`params`)](#parameters-params)
    - [input](#input)
    - [Pipeline flow controls](#pipeline-flow-controls)
    - [Kraken Database Parameters\*\*](#kraken-database-parameters)
    - [Kraken2Ref Handling\*\*](#kraken2ref-handling)
    - [Kraken2ref Report Filter](#kraken2ref-report-filter)
    - [iVar Parameters](#ivar-parameters)
    - [Virus Subtyping](#virus-subtyping)
    - [Consensus Generation](#consensus-generation)
    - [Output Directory](#output-directory)
    - [General](#general)
  - [Containerization](#containerization)
  - [Process Settings (`process`)](#process-settings-process)
  - [Profiles](#profiles)
    - [`sanger_standard` Profile](#sanger_standard-profile)
- [Unit Tests](#unit-tests)
- [Pipeline components documentation](#pipeline-components-documentation)
  - [Processes](#processes)
      - [`run_kraken`](#run_kraken)
      - [`run_k2r_sort_reads`](#run_k2r_sort_reads)
      - [`run_k2r_dump_fastqs_and_pre_report`](#run_k2r_dump_fastqs_and_pre_report)
      - [`concatenate_fqs_parts`](#concatenate_fqs_parts)
      - [`get_taxid_references`](#get_taxid_references)
      - [`bwa_alignment_and_post_processing`](#bwa_alignment_and_post_processing)
      - [`run_ivar`](#run_ivar)
      - [`run_pagolin`](#run_pagolin)
      - [`run_qc_script`](#run_qc_script)
      - [`write_classification_report`](#write_classification_report)
- [Licence](#licence)

---

## Pipeline Summary

The pipeline takes a manifest containing  **fastq pairs file** paths and a **kraken detabase** as inputs (check [Inputs section](#inputs) for more details) and outputs a **classification report**, **consensus sequences** and a **collection of intermediate files** (check [[ADD OUTPUT DOCUMENTATION LINK]]() for more details). Here is an broad overview of the pipeline logic

1. **Sort Reads**: The initial step is sort reads using `kraken2` for each fastq pairs according to the database provided. The classified reads is used as input to `kraken2ref` which will generate one pair of fastq files per taxid found.
   - An option to split big files is provided (check [Parameter section](#parameters)).
   - For a more in depth technical description check [`SORT_READS_BY_REF` workflow documentation](./docs/workflow/SORT_READS_BY_REF.md)

2. **Generate Consensus**: After all samples been classified, all references observed for that samples batch are fetch from the `kraken database` (or an arbitrary fasta file provided by the user). The classfied reads are aligned to their respective references (via `bwa`). The alignment is used as input for `ivar` to obtain a consensus sequence.
   - For a more in depth technical description check [`GENERATE_CONSENSUS` workflow documentation](./docs/workflow/GENERATE_CONSENSUS.md)
3. **Compute QC**: QC metrics are computed via `samtools` and a custom script (`qc.py`)
   - For a more in depth technical description check [`COMPUTE_QC_METRICS` workflow documentation](./docs/workflow/COMPUTE_QC_METRICS.md)

4. **SARS-CoV-2 Subtyping**: SARS-CoV-2 subtyping can be done if present on the sample
   - For a more in depth technical description check [`SCOV2_SUNTYPING` workflow documentation](./docs/workflow/SCOV2_SUBTYPING.md)

![](./docs/assets/metro_map.png)

---

[**(&uarr;)**](#contents)


## How to Cite

[ADD CITATION]

[**(&uarr;)**](#contents)

---
## Quick Start

Assuming [dependencies](#dependencies) are installed on the system:

1. Setup the pipeline:

```bash
# clone the repo
git clone https://gitlab.internal.sanger.ac.uk/malariagen1/viral/viral_pipeline.git
cd viral_pipeline/
# build containers
cd containers/
sudo singularity build base_container.sif baseContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
sudo singularity build pangolin.sif pangolinContainer.sing
sudo singularity build kraken.sif krakenContainer.sing
sudo singularity build kraken2ref.sif kraken2ref.sing
```

2. Run pipeline

You will need a manifest and a kraken database (check [Inputs section](#inputs) and [Usage section](#usage) for more details)

```bash
PIPELINE_CODES=<path to viral pipeline>
MANIFEST=<path to my manifest>
kraken_db_path=<path to my kraken DB>
PIPELINE_CONTAINERS=<path to my containers dir>

nextflow run ${PIPELINE_CODES}/main.nf --manifest ${MANIFEST} \
    --db_path ${kraken_db_path} \
    --outdir ./output/ \
    --containers_dir ${PIPELINE_CONTAINERS} \
    -with-trace -with-report -with-timeline \
    -profile sanger_local \
    -resume
```

[**(&uarr;)**](#contents)

---
## Installation

### dependencies

- [Nextflow](https://www.nextflow.io) (tested on `23.10.1`)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) (required to use Singularity Containers, tested on ``ce version 3.11.4``)

> We strongly recommend to run the pipeline using the containers recipe provided at `containers/` subdir.

If not using containers, all the software needs to be available at run time. Here is a list of all the softwares versions currently in use in the pipeline as set on each container.

- **Base container**

  - Samtools = `1.8`
  - BWA = `0.7.17`
  - biopython = `1.79`
  - pysam = `0.22.0`
  - pandas = `1.1.5`
  - matplotlib = `3.3.4`

- **Ivar container**
  - samtools = `1.11`
  - BWA = `0.7.17`
  - iVar = `1.4.3`

- **Kraken2Ref cotainer**:
  - pytest = `6.2.2`
  - importlib-resources = `5.1.0`
  - flake8 = `7.0.0`
  - pandas = `2.1.4`
  - cached-property = `1.5.2`
  - scipy = `1.12.0`
  - kraken2ref = `v2.0.0`

- **Kraken containers**
  - kraken2 = `v2.1.3`
  - kraken_tools = `v1.2`
  - biopython = `v1.81`

> NOTE: If not running under the provided containers, is strongly recommended to use the same versions described here.

### build containers

Recipes for the Singularity container used on this pipeline are available on this repository at `containers/` dir. To build the containers, run the commands bellow.

```{bash}
cd containers/
sudo singularity build base_container.sif baseContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
sudo singularity build pangolin.sif pangolinContainer.sing
sudo singularity build kraken.sif krakenContainer.sing
sudo singularity build kraken2ref.sif kraken2ref.sing
```

> NOTE: Currently we only support running on Singularity local containers. We should add Docker container registries.

[**(&uarr;)**](#contents)

## Usage

1. **Generate manifest**
For convenience, a script to generate the manifest is provided on this repo:

```{bash}
python write_manifest.py ./path/to/my/fastqs_dir/ -fq1_ext my_r1_ext -fq2_ext my_r2_ext
```

2. **Run pipeline**

```{bash}
nextflow run /path/to/rvi_consensus_gen/main.nf --manifest /path/to/my/manifest.csv \
        --db_path /path/to/my/kraken_db \
        --results_dir outputs/ \
        --containers_dir /path/to/my/containers_dir/ \
        -profile sanger_stantard -resume -with-trace -with-report -with-timeline
```

Optionally is possible to start the pipeline from **GENERATE_CONSENSUS**

```{bash}
nextflow run ${CHECKOUT}/main.nf --entry_point consensus_gen \
        --consensus_mnf sorted_manifest.csv
        --results_dir $LAUNCHDIR/outputs/ \
        --containers_dir /path/to/containers_dir/ \
        -profile sanger_standard-resume -with-trace -with-report -with-timeline
```

[**(&uarr;)**](#contents)

## Inputs

- `manifest` : CSV Manifest of input fastq file pairs.
  - Must have `sample_id`,`reads_1` and `reads_2` collumns
  - If you have your set of fastq pairs in a single dir, a script (`write_manifest.py`) is provided to facilitate this process.

- `db_path` : Path of a valid [kraken2 database](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#kraken-2-databases)

check [Input documentation](./docs/input.md) for more details.

> **NOTE**: if using the `consensus_gen` entry point, the manifest must containing pair ended reads and genome reference files for each pair. must have `sample_id`, `taxid`, `ref_files`, `reads_1`, `reads_2`

[**(&uarr;)**](#contents)

## Outputs

The **main output** files of the pipeline are:

- **Report CSV**: A csv file summarizing consensus sequences information and qc metrics obtained from each input samples.
  - location: `<output_dir>/classification_report.csv`
- **Consensus sequences per taxid**: the consensus sequence obtained for a given `taxid` and `sample_id` combination.
  -  location: `<output_dir/<sample_id>/<sample_id>.<taxid>.consensus.fa`

The pipeline outputs a set of **intermediate files** which are key to generate the information of the main output and can be helpful for futher investigation, for a more detailed description check [intermediate files documentation](./docs/output.docs).

The output file tree should look like the tree bellow: 
```bash
<output_dir>/
├── <sample_id>
│   ├── <taxid>
│   │   ├── <sample_id>.<taxid>.consensus.fa
│   │   ├── <sample_id>.<taxid>.qc.csv
│   │   ├── <sample_id>.<taxid>.sorted.bam
│   │   └── <sample_id>.<taxid>.sorted.bam.bai
│   ├── [...]
│   ├── <sample_id>.class_seqs_1.fq
│   ├── <sample_id>.class_seqs_2.fq
│   ├── <sample_id>.kraken.output
│   ├── <sample_id>.report.txt
│   ├── <sample_id>.unclass_seqs_1.fq
│   ├── <sample_id>.unclass_seqs_2.fq
├── [...]
├── classification_report.csv
```
[**(&uarr;)**](#contents)


## Configuration Sections

### Parameters (`params`)

The `params` block defines key user-modifiable settings for the workflow.

#### input
- `manifest`: Path to the manifest file (default: `null`).

#### Pipeline flow controls

- `entry_point` [Default = `"sort_reads"`]: Entry point for the workflow. It can be set to `"sort_reads"` or `"consensus_gen"`.
- `do_scov2_subtyping` [Default = `true`] : Switch SARS-CoV-2 on and off.
- `scv2_keyword` [Default = `"Severe acute respiratory syndrome coronavirus 2"`] : keyword to identified samples with SARS-CoV-2. This string is obtained from the kraken2 database, therefore, it should be in line with the database in use.

#### Kraken Database Parameters**

- `db_path`: Path to the Kraken database.
- `db_library_fa_path` (OPTIONAL): Path to the Kraken database library FASTA file.
  - By default, it assumes there is a `${params.db_path}/library/library.fna`.

#### Kraken2Ref Handling**

- `k2r_fq_load_mode`: Loading mode for Kraken2 fastq files (either `full` or `chunks`).
  - Default: `"full"`.
- `k2r_max_total_reads_per_fq`: Maximum number of reads to process per fastq file.
  - Default: `10,000,000`.
- `k2r_dump_fq_mem`: Memory allocated for dumping fastq files.
  - Default: `"6 GB"`.

#### Kraken2ref Report Filter

- `min_reads_for_taxid`: Minimum number of reads required to assign a taxonomic ID.
  - Default: `100`.

#### iVar Parameters

- `ivar_min_depth`: Minimum depth for consensus calling.
  - Default: `10`.
- `ivar_freq_threshold`: Sets base frequency threshold for consensus calling.
  - Default: `0.75`.

#### Virus Subtyping

- `scv2_keyword`: Keyword to identify SARS-CoV-2 sequences. Any taxid name equal to the string set by this parameter will be considered as SCOV2 and subjected to specific SARS-CoV-2 subtyping.
  - Default: `"Severe acute respiratory syndrome coronavirus 2"`.

- `do_scov2_subtyping`: Boolean flag to enable or disable SARS-CoV-2 subtyping via Pangolin. Check [SARS0COV2 documentation](./workflow/SCOV2_SUBTYPING.md) for more details
  - Default: `true`.

#### Consensus Generation

- `consensus_mnf` [OPTIONAL]: Path to the consensus manifest file (default: `null`).

#### Output Directory

- `results_dir`: Directory where output files will be published.
  - Default: `"$launchDir/output/"`.

#### General

- `containers_dir` [DEFAULT =  `containers/` dir of this repository] : By default, the pipeline relies on Singularity containers and __assumes__ all containers are present on this directory and were named on a specific manner
- `results_dir` [DEFAULT = `$launchDir/output/`] : set where output files should be published. By default, it will write files to an `output/` dir (if not existent, it will be created) at pipeline launch directory.

### Containerization

Currently, the pipeline only provide Singularity containers.

**Docker**

- `enabled`: Flag to enable or disable Docker.
  - Default: `false`.

**Singularity**

- `enabled`: Flag to enable or disable Singularity.
  - Default: `true`.

By default, the current version of the pipeline assumes all singularity containers are available on specific paths with specific names defined at `./conf/containers.config`.

- Currently, `"$projectDir/containers/"` is the default location for the containers, it can be changed by the user using `containers_dir` parameter.

All containers used by this pipeline recipes can be found at `./containers/` dir of this repository. The current containers in use are:

- `base_container.sif`: the default container for all processes unless overridden.
- `ivar.sif`: container to be used for processes using the `ivar` label
- `kraken.sif`: container to be used for processes using the `kraken` label
- `pangoling.sif`: container to be used for processes using the `pangolin` label
- `kraken2ref.sif`: container to be used for processes using the `kraken2ref` label.

### Process Settings (`process`)

- `cache='lenient'`: Defines the cache behavior for processes, allowing cached results to be reused even when minor changes occur.
- `executor='local'`: The default executor is set to local, meaning processes will run on the local machine unless otherwise specified.

### Profiles

Profiles define environment-specific configurations. Currently, there is a single predefined profile: `sanger_standard`. User should consider to write their own, check [profiles Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) for more details

#### `sanger_standard` Profile

This profile is optimized for running jobs on the Sanger Institute's infrastructure and includes settings for using Singularity, job execution on the LSF scheduler, and handling of specific processes.

**Singularity Settings**

- `autoMounts`: Automatically mounts paths required for job execution.
- `cacheDir`: Cache directory set to "`$PWD`" (the current working directory).
- `runOptions`: Singularity run options to bind necessary paths (`/lustre`, `/nfs`, `/software`, `/data/`).

**Process-Specific Settings**
Inherited from the global process configuration
- `cache='lenient'`: this makes `resume` more tolerant to changes in files attributes, such as timestamps. This is handy when using distributed file system which uses [Network File System protocols](https://en.wikipedia.org/wiki/Network_File_System), check [Nextflow documentation](https://www.nextflow.io/docs/latest/cache-and-resume.html#inconsistent-file-attributes) for more details.
- `executor='local'`: Default executor is local unless overridden.

**LSF Job Execution**
Certain processes are configured to run on the LSF scheduler with specific resource allocations:

- `run_k2r_sort_reads`:
  - Executor: `lsf`
  - CPU cores: `1`

- `run_k2r_dump_fastqs_and_pre_report`:
  - Executor: `lsf`
  - CPU cores: `1`

- `run_kraken`:
  - Executor: `lsf`
  - CPU cores: `16`
  - Memory: `1 GB`

**Job Naming and Memory Management**

- `jobName`: Custom job name format for LSF jobs, based on the task name and tag (`"RVI-preprocess - $task.name - $task.tag"`).

- `perJobMemLimit=true`: Ensures that memory limits are set per job.

[**(&uarr;)**](#contents)

---

## Unit Tests

The workflow & process unit tests for this pipeline are written in the [nf-test](https://www.nf-test.com/) (`v0.8.4`) Nextflow testing framework.

**Running Tests**

The following command if entered from the repository top-level directory can be used to execute all of the per-process & per-workflow unit tests:

**Run all tests**

```{bash}
nf-test test ./
```

**Run individual module/workflow test**

```{bash}
nf-test test tests/<modules or workflows>/<module_to_test>.nf.test
```

[**(&uarr;)**](#contents)

---
## Pipeline components documentation

For a more in depth technical documentation of all the processes and workflows can be found at `docs/` dir, here you will find documentation for:

### Processes


##### `run_kraken`

Executes Kraken2 on paired-end FASTQ files, producing outputs that include the classification results, classified and unclassified reads, and a summary report.

##### `run_k2r_sort_reads`

This process runs Kraken2Ref to parse Kraken reports and sort reads by taxonomic classification. It generates JSON files that map taxonomic IDs to read IDs and performs sorting based on the decomposed taxonomy tree if available.

##### `run_k2r_dump_fastqs_and_pre_report`

Extract classified reads into FASTQ files and generate a preliminary report based on the taxonomic classification data. It processes classified reads and produces a detailed report for further analysis.

##### `concatenate_fqs_parts`
This process concatenates FASTQ files from multiple parts into final combined FASTQ files for each taxonomic classification. This process ensures that all parts corresponding to the same taxonomic ID are merged into single files.

##### `get_taxid_references`

Retrieves sequences for a given taxid from a source FASTA file and indexes them for further analysis.

##### `bwa_alignment_and_post_processing`

Mapping sequencing reads to a reference genome using BWA (Burrows-Wheeler Aligner), followed by post-processing steps including conversion to BAM format, sorting, and indexing. At current version of this pipeline, this process generates high-quality, sorted BAM files that are essential for consensus sequence generation.

##### `run_ivar`

Generates a consensus sequence from a sorted BAM file using samtools mpileup and ivar consensus.


##### `run_pagolin`

The process runs the Pangolin tool on a consensus FASTA file to determine the SARS-CoV-2 lineage and extracts relevant metadata from the output.

##### `run_qc_script`

This process runs a QC analysis on the input BAM, FASTA, and reference files, outputting QC metrics.

##### `write_classification_report`

This process creates a classification report from provided report lines, ensuring the output is correctly formatted and free from common formatting issues such as unexpected character encodings.

[**(&uarr;)**](#contents)
**Workflow**

- [SORT_READS_BY_REF.nf](./docs/workflow/SORT_READS_BY_REF.md)
- [GENERATE_CONSENSUS.nf](./docs/workflow/GENERATE_CONSENSUS.md)
- [COMPUTE_QC_METRICS.nf](./docs/workflow/COMPUTE_QC_METRICS.md)
- [SCOV2_SUBTYPING.nf](./docs/workflow/SCOV2_SUBTYPING.md)

---

## Licence

[ADD LICENSE]
