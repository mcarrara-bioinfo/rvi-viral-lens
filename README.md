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

```
nextflow run main.nf --virus_resources_json ./virus_resources.json --reads_dir ./test_data
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
