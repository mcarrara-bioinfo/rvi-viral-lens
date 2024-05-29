# Viral Pipeline

This pipeline accepts `fastq` files as inputs and leverages `Kraken2` to perform read sorting, `ivar` for consensus sequence generation, and computes quality control (QC) metrics for each reference sequence identified in a given sample. Additionally, it facilitates SARS-CoV-2 subtyping. The primary output is a comprehensive classification report detailing all findings.

## Installation

### dependencies

- [Nextflow](https://www.nextflow.io) (tested on v23.10.0)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/)

### build containers

Recipes for the container used on this pipeline are available on this repository at `containers/` dir

```{bash}
cd containers/
sudo singularity build base_container.sif baseContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
sudo singularity build pangolin.sif pangolinContainer.sing
sudo singularity build kraken.sif krakenContainer.sing
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
        -profile sanger_local -resume -with-trace -with-report
```

### Running from **GENERATE_CONSENSUS**

```{bash}
nextflow run ${CHECKOUT}/main.nf --entry_point consensus_gen \
        --consensus_mnf sorted_manifest.csv     
        --results_dir $LAUNCHDIR/outputs/ \
        --containers_dir /path/to/containers_dir/ \
        -profile sanger_local -resume -with-trace -with-report
```

## Unit Tests

The workflow & process unit tests for this pipeline are written in the [nf-test](https://www.nf-test.com/) (v0.8.4) Nextflow testing framework. Nf-test will need to be installed to run the tests.

### Running Tests

The following command if entered from the repository top-level directory can be used to execute all of the per-process & per-workflow unit tests:

#### Run all tests

```{bash}
nf-test test
```

#### Run all module tests

```{bash}
nf-test test tests/modules/*.nf.test
```

#### Run all workflows tests

```{bash}
nf-test test tests/workflows/*.nf.test
```

#### Run whole pipeline test

```{bash}
nf-test test tests/main.nf.test
```

#### Run individual module/workflow test

```{bash}
nf-test test tests/<modules or workflows>/<module_to_test>.nf.test
```

## Input

### sort_reads

- kraken database path
- manifest for fastqs

### consensus_gen

- Manifest containing pair ended Reads and genome reference files for each pair
- manifest (required if starting from this entry point)
