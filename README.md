# rvi_consensus_gen


## Installation

### build container

```
cd containers/
sudo singularity build base_container.sif baseContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
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


[TODO] add manifest description

```

nextflow run /path/to/rvi_consensus_gen/main.nf --manifest /path/to/my/manifest.csv \
        --db_path /path/to/my/kraken_db \
        --outdir outputs/ \
        --containers_dir /path/to/my/containers_dir/ \
        --virus_resources_json /path/to/my/virus_resources.json \
        -profile sanger_local -resume
```

## Input

- Pair ended Reads
- virus_resources_json
```{json}
{
    "SC2": {
        "reference_genome": ["GCF_009858895.2."] // https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009858895.2/
    },
    "H1N1A": {
        "reference_genome": ["GCF_001343785.1"] // https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/
    },
    "H3N2A":{
        "reference_genome": ["GCF_000865085.1"] // https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000865085.1/
    },
    "fluB":{
        "reference_genome": ["GCF_000820495.2"] //https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000820495.2/
    },
    "RSV":{
        "reference_genome":["GCF_000856445.1"] // https://www.ncbi.nlm.nih.gov/datasets/taxonomy/12814/
    }
}
```


## Support


## License
For open source projects, say how it is licensed.
