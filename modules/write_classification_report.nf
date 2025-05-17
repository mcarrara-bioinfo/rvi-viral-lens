// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
process write_classification_report {
    /*
    * ---------------------------------------------------------------
    * Write classification report

    The process in this Nextflow pipeline generates a classification
    report in CSV format by writing provided lines of data into the
    report file. This report summarizes various metrics and
    classifications related to the sequencing analysis, including
    information on sample identification, taxonomic classification,
    genome coverage, and more.

    * ---------------------------------------------------------------
    * Input
        - A list of lines that contain data to be written into the
        classification report. Each line corresponds to a data entry
        related to the classification results of samples.
        - A minimum threshold of Total_Mapped_Reads to include the
        classification line in the report (Default: 100)
    * Output
        - Outputs the path to the generated classification report file
        (`classification_report.csv`).

    * ---------------------------------------------------------------
    */

    publishDir "${params.outdir}/", mode: 'copy'

    input:
        val(list_of_report_lines)
        val(min_total_mapped_reads)

    output:
        path(output_report_file)

    script:
        output_report_file = "classification_report.csv"
        // Replace " with ' to prevent issues with writing lines to file 
        report_lines = list_of_report_lines.join("").replaceAll(/"/, "'")
        """
        # Write header to output report file
        echo "Sample_ID,Virus_Taxon_ID,Virus,Species,Reference_Taxon_ID,Selected_Reference,Flu_Segment,Reference_Subtype,Sample_Subtype,Percentage_of_Genome_Covered,Total_Mapped_Reads,Longest_non_N_segment,Percentage_of_N_bases" > ${output_report_file}_pre

        # Write data lines to report file
        echo "${report_lines}" >> ${output_report_file}_pre
        sed -e "s/\r//g" ${output_report_file}_pre | \
            awk -v min="$min_total_mapped_reads" -F ',' ' \
                NR==1 { print; next } \    # Always print the header
                \$11 >= min { print } \    # Print report line only if column 11 >= min
            ' > ${output_report_file}
        # NOTE: the sed expression is there to remove "^M" added characteres
        """
/*
    * ---------------------------------------------------------------
# Script Breakdown

1. **Output File Name**:
The output report file is named classification_report.csv.

2. **Replace Double Quotes**: 
    - `report_lines = list_of_report_lines.join("").replaceAll(/"/, "'")`
    - This line of Groovy script joins all the report lines into a
    single string and replaces double quotes (") with single quotes
    ('). This substitution helps prevent potential formatting issues
    when writing to the CSV file.

3. **Write Header**:

    - Command: 
    ```
    echo "Sample_ID,...,Percentage_of_N_bases" > ${output_report_file}_pre`
    ```
    - Writes the header line to the pre-output file 
    (${output_report_file}_pre). The header defines the columns in the
    CSV file, which include various identifiers and metrics related to
    the samples and their classifications.

4. **Write Data Lines**:

    - Command: `echo "${report_lines}" >> ${output_report_file}_pre`
    - Appends the processed data lines (${report_lines}) to the pre-output file.

5. **Remove Carriage Return Characters**:

    - Command: 
    ```
    sed -e "s/\r//g" ${output_report_file}_pre | \
    ```
    - This command uses sed to remove any carriage return (\r) characters
    that may be present in the file. These characters can be introduced 
    when handling files across different operating systems (e.g., 
    Windows vs. Unix-based systems), and their removal ensures 
    consistent file formatting.
    The output is directly piped to the next command.

6. **Filter out any report line with Total_Mapped_Reads less than `min_total_mapped_reads`**

    - Command:
    ```
    awk -v min="$min_total_mapped_reads" -F ',' ' \
        NR==1 { print; next } \    # Always print the header
        \$11 >= min { print } \    # Print report line only if column 11 >= min
    ' > ${output_report_file}

    - This command uses awk to remove any report line where the value of the
    column `Total_Mapped_Reads` is less than the provided `min_total_mapped_reads`
    threshold. The code ensures that the header is always included in the output.
    The input is a direct pipe from the sed command at point 5.

*/
}
