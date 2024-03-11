# rvi_consensus_gen


## Installation

### build container

```
cd containers/
sudo singularity build base_container.sif baseContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
sudo singularity build pangolin.sif pangolinContainer.sing
# sudo singularity build kraken.sif kraken_container.singularity
```
### prepare reference genomes

```
cd virus_resources/
bwa index NC_045512.2.fasta
```

## Unit Tests
The workflow & process unit tests for this pipeline are written in the [nf-test](https://www.nf-test.com/) Nextflow testing framework. [nf-test](https://www.nf-test.com/) will need to be installed to run the tests.

### Running Tests
The following command if entered from the repository top-level directory can be used to execute all of the per-process & per-workflow unit tests:

#### Run all tests
```
nf-test test
```

#### Run all module tests
```
nf-test test tests/modules/*.nf.test
```

#### Run all workflows tests
```
nf-test test tests/workflows/*.nf.test
```

#### Run whole pipeline test 
```
nf-test test tests/main.nf.test
```

#### Run individual module/workflow test
```
nf-test test tests/<modules or workflows>/<module_to_test>.nf.test
```

## Usage

#### Generate manifest

For convenience, a script to generate the manifest is provided on this repo

```
python write_manifest.py ./path/to/my/fastqs_dir/ -fq1_ext my_r1_ext -fq2_ext my_r2_ext
```
#### Run pipeline

```
nextflow run /path/to/rvi_consensus_gen/main.nf --manifest /path/to/my/manifest.csv \
        --db_path /path/to/my/kraken_db \
        --results_dir outputs/ \
        --containers_dir /path/to/my/containers_dir/ \
        -profile sanger_local -resume -with-trace -with-report
```

> [TODO] add manifest description

#### Running from **consensus_gen**
```
nextflow run ${CHECKOUT}/main.nf --entry_point consensus_gen \
        --consensus_mnf sorted_manifest.csv     
        --results_dir $LAUNCHDIR/outputs/ \
        --containers_dir /path/to/containers_dir/ \
        -profile sanger_local -resume -with-trace -with-report
```

> [TODO] add manifest description

## Input

### sort_reads

- kraken database path
- manifest for fastqs

### consensus_gen
- Pair ended Reads
- virus resources json file
- manifest (required if starting from this entry point) 
