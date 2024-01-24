//include {generate_classification_report_line} from '../modules/generate_classification_report_line.nf'
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
        "${meta.id},${meta.virus_species},${meta.type},${meta.flu_segment},${meta.percentage_genome_coverage},${meta.read_depth}\n"
        }.collect()

        // Write all of the per-sample report lines to a report file
        write_classification_report(report_lines_ch)

    emit:
        write_classification_report.out
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
