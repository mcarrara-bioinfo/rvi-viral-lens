include {write_classification_report} from '../modules/write_classification_report.nf'

workflow GENERATE_CLASSIFICATION_REPORT {
    take:
        /*
        meta must have the following keys:
            - sample_id
            - taxid
            - ref_selected
            - virus_name
            - virus_subtype
            - flu_segment
            - percentage_genome_coverage
            - total_mapped_reads
            - longest_no_N_segment
            - percentage_of_N_bases
        */
        meta_ch // meta

    main:
    
        // Create a report line for every sample and then aggregate them
        report_lines_ch = meta_ch.map{it ->
            // convert null values for type and segments to empty strings
            if (it[0].virus_subtype == null){
                virus_subtype='None'
            } else {
                virus_subtype = it[0].virus_subtype
            }

            if (it[0].flu_segment==null){
                flu_segment='None'
            } else {
                flu_segment = it[0].flu_segment
            }
            "${it[0].sample_id},${it[0].taxid},${it[0].ref_selected.replace(",","|")},${it[0].virus_subtype},${flu_segment},${it[0].percentage_genome_coverage},${it[0].total_mapped_reads.replace("^M", "")},${it[0].longest_no_N_segment},${it[0].percentage_of_N_bases}\n"
        }.collect()

        // Write all of the per-sample report lines to a report file
        write_classification_report(report_lines_ch)

}

workflow {
    manifest_channel = Channel.fromPath(params.manifest_file)
    | splitCsv(header: true, sep: ',')
    | map { row ->
        meta = [[id:row.sample_id,
            taxid:row.taxid,
            ref_selected:row.ref_selected,
            virus_name:row.virus_name,
            virus_subtype:row.virus_subtype,
            flu_segment:row.flu_segment,
            percentage_genome_coverage:row.percentage_genome_coverage,
            total_mapped_reads:row.total_mapped_reads,
            longest_no_N_segment:row.longest_no_N_segment,
            percentage_of_N_bases:row.percentage_of_N_bases
        ]]
    }

    GENERATE_CLASSIFICATION_REPORT(manifest_channel)
}
