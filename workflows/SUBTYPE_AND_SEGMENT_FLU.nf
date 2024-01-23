include {retrieve_flu_subtype_and_segment} from '../modules/retrieve_flu_subtype_and_segment.nf'

workflow SUBTYPE_AND_SEGMENT_FLU {
    take:
        /*
        meta must have the following keys:
            - taxid
        */
        sample_kraken_report_channel // tuple (meta, kraken_report)

    main:

        // Where possible, retrieve the flu type and flu segment from the kraken report file
        retrieve_flu_subtype_and_segment(sample_kraken_report_channel)

        // Add the flu type and segment to the metadata map
        subtyped_sample_kraken_report_channel = retrieve_flu_subtype_and_segment.out
        | map{meta, kraken_report, type, flu_segment -> 
            meta.type = type
            meta.flu_segment = flu_segment
            tuple (meta, kraken_report)
        }

    emit:
        subtyped_sample_kraken_report_channel // tuple (meta, kraken_report)
}   

workflow {
    manifest_channel = Channel.fromPath(params.manifest_file)
    | splitCsv(header: true, sep: ',')
    | map { row -> 
    meta = [taxid:row.taxon_id]
    tuple (meta, row.kraken_report_file_path)}

    SUBTYPE_AND_SEGMENT_FLU(manifest_channel)
}
