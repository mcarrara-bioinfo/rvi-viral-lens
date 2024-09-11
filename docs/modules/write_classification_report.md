# Nextflow Process Documentation: `write_classification_report`

## Overview

The process in this Nextflow pipeline generates a classification report in CSV format by writing provided lines of data into the report file. This report summarizes various metrics and classifications related to the sequencing analysis, including information on sample identification, taxonomic classification, genome coverage, and more.

## Process

This process creates a classification report from provided report lines, ensuring the output is correctly formatted and free from common formatting issues such as unexpected character encodings. Below is a detailed breakdown of the process components:

### Tags and Labels

- Tag: No explicit tag is defined in this process.
- Labels: No explicit tag is defined in this process


### Input

- Input Value: `val(list_of_report_lines)`
  - A list of lines that contain data to be written into the classification report. Each line corresponds to a data entry related to the classification results of samples.

### Output

- Output Path: `path(output_report_file)`
  - Outputs the path to the generated classification report file (`classification_report.csv`).

### Publish Directory

- Publish Directory: `${params.results_dir}/`
  - The report is published in the results directory specified by the pipeline parameters, making the output accessible in a centralized location.

### Script

The script section handles the generation of the CSV report:

Groovy:
```groovy
output_report_file="classification_report.csv"

// Replace double quotes with single quotes to avoid formatting issues
report_lines=list_of_report_lines.join("").replaceAll(/"/, "'")
```

Bash:
```bash
# Write header to the output report file
echo "Sample_ID,Virus_Taxon_ID,Virus,Species,Reference_Taxon_ID,Selected_Reference,Flu_Segment,Reference_Subtype,Sample_Subtype,Percentage_of_Genome_Covered,Total_Mapped_Reads,Longest_non_N_segment,Percentage_of_N_bases" > ${output_report_file}_pre

# Write data lines to the report file
echo "${report_lines}" >> ${output_report_file}_pre

# Remove any unwanted carriage return characters that may have been introduced
sed -e "s/\r//g" ${output_report_file}_pre > ${output_report_file}
# NOTE: The sed command is used to remove "^M" characters often introduced by different line endings
"""
```

### Script Breakdown

1. **Output File Name**:
The output report file is named classification_report.csv.

2. **Replace Double Quotes**: `report_lines = list_of_report_lines.join("").replaceAll(/"/, "'")`
This line of Groovy script joins all the report lines into a single string and replaces double quotes (") with single quotes ('). This substitution helps prevent potential formatting issues when writing to the CSV file.

3. **Write Header**:
Command: `echo "Sample_ID,...,Percentage_of_N_bases" > ${output_report_file}_pre`
Writes the header line to the pre-output file (${output_report_file}_pre). The header defines the columns in the CSV file, which include various identifiers and metrics related to the samples and their classifications.

4. **Write Data Lines**:
Command: `echo "${report_lines}" >> ${output_report_file}_pre`
Appends the processed data lines (${report_lines}) to the pre-output file.

5. **Remove Carriage Return Characters**:
Command: `sed -e "s/\r//g" ${output_report_file}_pre > ${output_report_file}`
This command uses sed to remove any carriage return (\r) characters that may be present in the file. These characters can be introduced when handling files across different operating systems (e.g., Windows vs. Unix-based systems), and their removal ensures consistent file formatting.
