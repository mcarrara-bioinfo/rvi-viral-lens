process write_classification_report {
    /*
    * Write classification report
    *
    * The process in this Nextflow pipeline generates
    * a classification report in CSV format by
    * writing provided lines of data into the report
    * file. This report summarizes various metrics
    * and classifications related to the sequencing
    * analysis, including information on sample
    * identification, taxonomic classification,
    * genome coverage, and more.
    *
    * check docs/modules/write_classification_report.md for more extensive documentation
    */

    publishDir "${params.results_dir}/", mode: 'copy'

    input:
        val(list_of_report_lines)

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
        sed -e "s/\r//g" ${output_report_file}_pre > ${output_report_file}
        # NOTE: the sed expression is there to remove "^M" added characteres
        """
}
