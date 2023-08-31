# rvi_consensus_gen


## Installation

### build container

```
cd containers/
sudo singularity build base_container.sif baseContainer.sing 
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

## Support


## License
For open source projects, say how it is licensed.
