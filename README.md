[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
---
# Viral Pipeline

This pipeline accepts `fastq` files as inputs and leverages `Kraken2` to perform read sorting, `ivar` for consensus sequence generation, and computes quality control (QC) metrics for each reference sequence identified in a given sample. Additionally, it facilitates SARS-CoV-2 subtyping. The primary output is a comprehensive classification report detailing all findings.

## Pipeline introduction

---

## Pipeline Summary

---

## Quitck Start

---

## Installation

---
### dependencies


- [Nextflow](https://www.nextflow.io) (tested on v23.10.0)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) (required to use Singularity Containers)

> We strongly recommend to run the pipeline using the containers recipe provided at `containers/` subdir.

If not using containers, all the software needs to be available at run time. A list of this softwares can be found [here](./doc/requirements.md)

### build containers

#### Singularity

Recipes for the container used on this pipeline are available on this repository at `containers/` dir. To build the containers, run the command bellow.


```{bash}
cd containers/
sudo singularity build base_container.sif baseContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
sudo singularity build pangolin.sif pangolinContainer.sing
sudo singularity build kraken.sif krakenContainer.sing
sudo singularity build kraken2ref.sif kraken2ref.sing
```

## Usage

### Generate manifest

For convenience, a script to generate the manifest is provided on this repo

```{bash}
python write_manifest.py ./path/to/my/fastqs_dir/ -fq1_ext my_r1_ext -fq2_ext my_r2_ext
```

### Run pipeline

```{bash}
nextflow run /path/to/rvi_consensus_gen/main.nf --manifest /path/to/my/manifest.csv \
        --db_path /path/to/my/kraken_db \
        --results_dir outputs/ \
        --containers_dir /path/to/my/containers_dir/ \
        -profile sanger_stantard -resume -with-trace -with-report -with-timeline
```

### Running from **GENERATE_CONSENSUS**

```{bash}
nextflow run ${CHECKOUT}/main.nf --entry_point consensus_gen \
        --consensus_mnf sorted_manifest.csv
        --results_dir $LAUNCHDIR/outputs/ \
        --containers_dir /path/to/containers_dir/ \
        -profile sanger_standard-resume -with-trace -with-report -with-timeline
```


## Required Input

- `manifest` : Manifest if input fastq files
- `db_path` : Path of kraken2 database


### sort_reads

- kraken database path
- manifest for fastqs

### consensus_gen

- Manifest containing pair ended Reads and genome reference files for each pair
- manifest (required if starting from this entry point)

## Parameters

### Input data

- `manifest` [REQUIRED] : Manifest containing sample id and fastq pairs paths.

### General

- `containers_dir` [DEFAULT =  `containers/` dir of this repository] : By default, the pipeline relies on Singularity containers and __assumes__ all containers are present on this directory and were named on a specific manner
- `results_dir` [DEFAULT = `$launchDir/output/`] : set where output files should be published. By default, it will write files to an `output/` dir (if not existent, it will be created) at pipeline launch directory.

### Pipeline flow controls

- `entry_point` [Default = "sort_reads"] : Defines the entry point of the pipeline. The `sort_reads`
- `scv2_keyword` [Default = "Severe acute respiratory syndrome coronavirus 2"] : keyword to identified samples with SARS-CoV-2. This string is obtained from the kraken2 database, therefore, it should be in line with the database in use.
- `do_scov2_subtyping` [Default = true] : Switch SARS-CoV-2 on and off.

### Specific processes

- `db_library_fa_path` [Default = null]:

#### Kraken2ref

- `k2r_fq_load_mode` [Default = "full"] : kraken2ref load fastq into memory mode ["full"/"chunk"]
- `k2r_max_total_reads_per_fq` [Default = 10000000] : set maximum number of reads to be accepted any classified fastq file pair. Files excedding that limit will be splitted before `run_k2r_dump_fastq`.
- `k2r_dump_fq_mem` [Default = "6 GB"] : Memory to be requested by `run_k2r_dump_fq`, adjust according to `k2r_max_total_reads_per_fq`.
- `min_reads_for_taxid` [Default = 100] : min reads number to be considered.

Kraken2ref have a escalation memory strategy based on linear regression.

### ivar params

- `ivar_min_depth` [Default = 10] : Minimum depth to call consensus
- `ivar_freq_threshold` [Default = 0.75] : Minimum frequency threshold(0 - 1) to call consensus.

### 
> TODO: add documentation for all pipeline parameters

## Unit Tests
---
The workflow & process unit tests for this pipeline are written in the [nf-test](https://www.nf-test.com/) (v0.8.4) Nextflow testing framework. Nf-test will need to be installed to run the tests.

### Running Tests

The following command if entered from the repository top-level directory can be used to execute all of the per-process & per-workflow unit tests:

#### Run all tests

```{bash}
nf-test test ./
```

#### Run individual module/workflow test

```{bash}
nf-test test tests/<modules or workflows>/<module_to_test>.nf.test
```
