# PIPELINE INPUT DOCUMENTATION

This pipeline relies on two **main inputs**:

- **`manifest`** : CSV Manifest of input fastq file pairs.
  - Must have `sample_id`,`reads_1` and `reads_2` collumns
  - If you have your set of fastq pairs in a single dir, a script (`write_manifest.py`) is provided to facilitate this process.

- **`db_path`** : Path of a valid [kraken2 database](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#kraken-2-databases)

## Manifest


The pipeline require as input a manifest containing a unique sample id (`sample_id`) and paths to each of the fastq pair file (`reads_1` and `reads_3`)

```csv
sample_id,reads_1,reads_2
sample1,/path/to/output/sample1_R1.fq,/path/to/output/sample1_R2.fq
sample2,/path/to/output/sample2_R1.fq,/path/to/output/sample2_R2.fq
sample3,/path/to/output/sample3_R1.fq,/path/to/output/sample3_R2.fq
```

**Write Manifest Script**

For user convenience, a script to write a manifest (`write_manifest.py`). This script generates a CSV manifest file from a directory of FASTQ files.

1. It takes a glob pattern to locate files, checks for paired-end FASTQ files (based on provided extensions)
2. extracts metadata such as sample IDs (optionaly taxonomic IDs) from the filenames.
    1. The user can specify if taxonomic IDs are included in the filenames, along with the filename separator. 
    2. If reference files are provided for taxonomic IDs, they are linked to the samples in the output. 

The script handles errors such as missing paired-end files or improperly formatted filenames, and generates a CSV manifest that includes sample information, file paths, and reference data (if applicable).

Example Command:

```bash
python write_manifest.py "output/*/reads_by_taxon/*.extracted_{1,2}.fq" \
    -fq1_ext "_R1.fq" \
    -fq2_ext "_R2.fq" \
    -filename_sep "." \
    -manifest_out "/path/to/output/manifest.csv"
```

**Breakdown**:

- `"output/*/reads_by_taxon/*.extracted_{1,2}.fq"` - A glob pattern to locate the FASTQ files in the specified directory structure.
- `-fq1_ext "_R1.fq"` - Specifies the naming pattern for the forward FASTQ files (e.g., `_R1.fq`).
- `-fq2_ext "_R2.fq"` - Specifies the naming pattern for the reverse FASTQ files (e.g., `_R2.fq`).
- `-filename_sep "."` - The separator used in the filename to extract sample id (e.g., `.`), anything before the separator will be considered as the value which should be the `sample_id`
- `-manifest_out "/path/to/output/manifest.csv"` - The output path for the generated CSV manifest file.

## Kraken Database

The pipeline only requirement assumes a valid Kraken Database.
However, this pipeline was developed under the RVI project and the following Database were developed to be used on this pipeline

**[ADD DESCRIPTION OF THE KRAKEN DATABASE SPECIFICS]**