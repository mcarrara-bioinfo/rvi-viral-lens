# Kraken2ref JSON to TSV Report Script

## Overview

This Python script is designed to convert a Kraken2ref JSON output into a tab-separated values (TSV) report file. The script reads two main input files: the Kraken2ref JSON file and a corresponding Kraken2 taxonomic report. It then extracts relevant information about selected reference taxa, including virus subtypes and the number of reads per taxon, and generates a TSV report summarizing this data. The report provides detailed information on the sample, viruses, selected taxonomic IDs, flu segments, and subtyping data, specifically for influenza viruses if present.

## Inputs

1. **Kraken2ref JSON File** (`--input_json`): The JSON output produced by Kraken2ref that contains metadata and output information about selected references and associated viruses.
2. **Kraken2 Taxonomic Report** (`--report`): The Kraken2 report, which provides the taxonomic breakdown of the sample.

3. **Output Suffix** (`--out_suffix`, optional): A suffix for the output file. Defaults to `.viral_pipe.report.tsv`.

## Outputs

The script generates a TSV report file summarizing the viral taxonomic information for each selected reference. The output file is named after the sample ID followed by the specified suffix (e.g., `sample123.viral_pipe.report.tsv`).

## Key Features

- **Taxonomic and Virus Information Extraction**: The script identifies the virus taxonomic IDs and names from both the Kraken2ref JSON and Kraken2 report, allowing detailed annotation of selected viruses and reference sequences.
- **Influenza Subtyping**: If influenza is detected, the script extracts subtype information (e.g., H and N subtypes) and flu segment numbers. These details are specifically captured for influenza viruses, helping to identify the sample subtype.
- **Reads Per Taxon**: The script reports the number of reads assigned to each selected taxonomic ID, providing insights into the abundance of each virus in the sample.
Customizable Output: Users can specify the output file suffix, giving flexibility in naming the report files.

## Workflow

1. **Input Parsing**: The script takes input arguments using argparse for the input JSON, Kraken2 report, and optional output suffix.

2. **Data Extraction**:

   - Reads and parses both the Kraken2ref JSON and Kraken2 taxonomic report.
   - Extracts virus taxonomic IDs and their human-readable names.
   - Identifies influenza viruses and retrieves subtyping information based on flu segments (H/N subtypes for flu segments 4 and 6).

3. **TSV Report Generation**: The extracted data is organized into a structured table and written to a TSV file.

## Example Command

```bash
./kraken2ref_to_tsv.py -i kraken2ref_output.json -r kraken2_report.txt --out_suffix ".custom_report.tsv"
```

This command reads the JSON and report files, processes the data, and outputs a report named `<sample_id>.custom_report.tsv`.
