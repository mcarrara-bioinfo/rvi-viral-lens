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
        --outdir outputs/ \
        --containers_dir /path/to/my/containers_dir/ \
        -profile sanger_local -resume -with-trace -with-report
```

> [TODO] add manifest description

#### Running from **consensus_gen**
```
nextflow run ${CHECKOUT}/main.nf --entry_point consensus_gen \
        --consensus_mnf sorted_manifest.csv     
        --outdir $LAUNCHDIR/outputs/ \
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
