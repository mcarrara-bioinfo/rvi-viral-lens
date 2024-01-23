include {generate_report_line} from '../modules/generate_report_line.nf'
include {write_report} from '../modules/write_report.nf'

workflow GENERATE_REPORT {
    take:
        /*
        meta must have the following keys:
            - id
            - taxid
            - virus_species
            - type
            - segment
            - ref_files
            - bam_file_path
        */
        meta_ch // meta

    main:

        // Generate per-sample report lines
        reports_input_ch = meta_ch.map{meta -> tuple(meta, meta.ref_files)}
        generate_report_line(reports_input_ch)

        // Aggregate all of per-sample report lines
        report_lines_ch = generate_report_line.out.collect()

        // Write all of the per-sample report lines to a report file
        write_report(report_lines_ch)

    emit:
        write_report.out
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
            ref_files:row.ref_files, 
            bam_file:row.bam_file_path
        ]}

    GENERATE_REPORT(manifest_channel)
}
