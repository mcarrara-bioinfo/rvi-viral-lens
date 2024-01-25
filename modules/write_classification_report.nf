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
        echo "Sample ID,Virus Species,Type,Flu Segment","Percentage of genome coverage","Read Depth" > ${output_report_file}

        # Write data lines to report file
        echo '${report_files}' >> ${output_report_file}
        """
}
