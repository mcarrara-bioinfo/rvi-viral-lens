process write_classification_report {

    publishDir "${params.results_dir}/", mode: 'copy'

    input:
        val(list_of_report_files)

    output:
        path(output_report_file)

    script:
        output_report_file = "classification_report.csv"
        report_files = list_of_report_files.join("")
        """
        # Write header to output report file
        echo "Sample_ID,Taxon_ID,Virus,Type,Flu_Segment,Percentage_of_Genome_Covered,Total_Mapped_Reads,Longest_non_N_segment,Percentage_of_N_bases" > ${output_report_file}_pre

        # Write data lines to report file
        echo '${report_files}' >> ${output_report_file}_pre
        sed -e "s/\r//g" ${output_report_file}_pre > ${output_report_file}
        # NOTE: the sed expression is there to remove "^M" added characteres
        """
}
