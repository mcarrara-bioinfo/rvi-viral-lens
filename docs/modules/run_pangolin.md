# Nextflow Process Documentation: `run_pangolin`

## Overview

The `run_pangolin` process in this Nextflow pipeline is designed to classify SARS-CoV-2 sequences into lineages using the Pangolin tool. [Pangolin (Phylogenetic Assignment of Named Global Outbreak LINeages)](https://github.com/cov-lineages/pangolin) is a tool widely used in the analysis of SARS-CoV-2 genomes to determine the likely lineage of the virus based on the consensus sequence.

## Process: `run_pangolin`

The process runs the Pangolin tool on a consensus FASTA file to determine the SARS-CoV-2 lineage and extracts relevant metadata from the output. Below is a detailed breakdown of the process components:

### Tags and Labels

- Tag: `${meta.id}` – This tag identifies the process run by the sample ID, helping with tracking and logging within the pipeline.
- Label: `"pangolin"` – This label is used for resource configuration, at current version this label only sets which container to be used. For more information check the [Labels documentation[TODO]]().

### Input

- Input Tuple: `tuple val(meta), path(mapped_fasta)`
  - `meta`: Metadata associated with the sample, including identifiers like sample ID and taxonomic ID.
  - `mapped_fasta`: The path to the consensus FASTA file generated for the sample.

### Output

- Output Tuple: `tuple val(meta), path(mapped_fasta), env(lineage)`

  - `mapped_fasta`: The path to the input consensus FASTA file.
  - `env(lineage)`: Environment variables capturing the lineage assignment information, which are extracted from the Pangolin output.

### Publish Directory

- Publish Directory: `${params.results_dir}/${meta.sample_id}/${meta.taxid}/`
  - The results are published in a structured directory organized by sample ID and taxonomic ID, making it easy to locate outputs.

### Shell Script

```bash
lineage_report="${meta.id}_lineage.csv"
consensus_fasta=mapped_fasta

# Run Pangolin to classify the lineage
pangolin "${consensus_fasta}" --outfile "${lineage_report}"

# Extract lineage information from the Pangolin output
# This assumes the report has only one row (excluding the header)
IFS=, read -r taxon lineage conflict ambiguity_score scorpio_call \
scorpio_support scorpio_conflict scorpio_notes version \
pangolin_version scorpio_version constellation_version \
is_designated qc_status qc_notes note < <(tail -n +2 "${lineage_report}")

# Output lineage assignment to logs for verification
echo "Taxon: $taxon"
echo "Lineage: $lineage"
echo "Conflict: $conflict"
```

### Script Breakdown

1. **Run Pangolin**:

- Command: `pangolin "${consensus_fasta}" --outfile "${lineage_report}"`
- This command runs Pangolin on the input consensus FASTA file and writes the results to a CSV file named with the sample ID and `_lineage.csv` suffix.

2. **Extract Lineage Information**:

- The script uses the IFS (Internal Field Separator) command to parse the CSV output file from Pangolin, extracting values from each column into environment variables.
- The script assumes that there is only one row of data (excluding the header) in the lineage report and reads this data into variables such as taxon, lineage, conflict, ambiguity_score, and others.

3. **Log Lineage Details**:

- Lineage details such as taxon, lineage, and conflict are echoed to in order to be captured as strings on the output channel.
  - Curretly, only `$lineage` is used for downstream analysis.

> Note: The `Taxon` and `Conflict` values are echoed in order to have those value within easy access on the `.command.log` files, it was usefull for debugging during development, but may be not necessary to keep it anymore. We should consider to remove it.
