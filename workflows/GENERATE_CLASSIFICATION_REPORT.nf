include {write_classification_report} from '../modules/write_classification_report.nf'

workflow GENERATE_CLASSIFICATION_REPORT {
    take:
        /*
        meta must have the following keys:
            - id
            - taxid
            - virus_species
            - type
            - segment
            - percentage_genome_coverage
            - read_depth
        */
        meta_ch // meta

    main:

        // Create a report line for every sample and then aggregate them
        report_lines_ch = meta_ch.map{meta ->
        // convert null values for type and segments to empty strings
        if (meta[0].virus_subtype == null){
            virus_subtype=''
        } else {
            virus_subtype = meta[0].virus_subtype
        }

        if (meta[0].flu_segment==null){
            flu_segment=''
        } else {
            flu_segment = meta[0].flu_segment
        }

        "${meta[0].sample_id},${meta[0].taxid},${meta[0].virus_name.replace(",","|")},${meta[0].virus_subtype},${flu_segment},${meta[0].percentage_genome_coverage},${meta[0].total_mapped_reads.replace("^M", "")},${meta[0].longest_no_N_segment},${meta[0].percentage_of_N_bases}\n"
        }.collect()

        // Write all of the per-sample report lines to a report file
        write_classification_report(report_lines_ch)

}

workflow {
    manifest_channel = Channel.fromPath(params.manifest_file)
    | splitCsv(header: true, sep: ',')
    | map { row ->
    meta = [id:row.sample_id,
            taxid:row.taxid,
            virus_species:row.virus_species,
            type:row.type,
            flu_segment:row.flu_segment,
            percentage_genome_coverage:row.percentage_genome_coverage,
            read_depth:row.read_depth
    ]}

    GENERATE_CLASSIFICATION_REPORT(manifest_channel)
}
