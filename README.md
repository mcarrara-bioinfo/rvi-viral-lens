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
- [Parameters](#parameters)
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

## Parameters

- **General**

  - `containers_dir` [DEFAULT =  `containers/` dir of this repository] : By default, the pipeline relies on Singularity containers and __assumes__ all containers are present on this directory and were named on a specific manner
  - `results_dir` [DEFAULT = `$launchDir/output/`] : set where output files should be published. By default, it will write files to an `output/` dir (if not existent, it will be created) at pipeline launch directory.

- **Pipeline flow controls**

  - `entry_point` [Default = `"sort_reads"`] : Defines the entry point of the pipeline. The `sort_reads`
  - `scv2_keyword` [Default = `"Severe acute respiratory syndrome coronavirus 2"`] : keyword to identified samples with SARS-CoV-2. This string is obtained from the kraken2 database, therefore, it should be in line with the database in use.
  - `do_scov2_subtyping` [Default = `true`] : Switch SARS-CoV-2 on and off.

- **Specific processes**

  - `get_taxid_references`:
    - `db_library_fa_path` [Default = `null`]: A fasta file containing reference sequences for the taxids present on the kraken database. It assumes there is a `$db_`
  - `raken2ref`:
    - `k2r_fq_load_mode` [Default = `"full"`] : kraken2ref load fastq into memory mode ["full"/"chunk"]
    - `k2r_max_total_reads_per_fq` [Default = `10000000`] : set maximum number of reads to be accepted any classified fastq file pair. Files excedding that limit will be splitted before `run_k2r_dump_fastq`.
    - `k2r_dump_fq_mem` [Default = `"6 GB"`] : Memory to be requested by `run_k2r_dump_fq`, adjust according to `k2r_max_total_reads_per_fq`.
    - `min_reads_for_taxid` [Default = `100`] : min reads number to be considered.
  - `ivar`
    - `ivar_min_depth` [Default = `10`] : Minimum depth to call consensus
    - `ivar_freq_threshold` [Default = `0.75`] : Minimum frequency threshold(0 - 1) to call consensus.
  - `qc`
    - `qc_minimum_depth` [Default = `10`] : Minimum depth value to be used for filtering when counting the number of covered positions, optional.

> NOTE: `Kraken2ref` have a escalation memory strategy based on linear regression, check [k2r_memory_escalation documentation](./docs/k2r_memory_escalation.md) for more details.

[**(&uarr;)**](#contents)

## Unit Tests
---
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
