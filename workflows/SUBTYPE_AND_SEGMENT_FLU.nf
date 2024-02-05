include {retrieve_flu_subtype_and_segment} from '../modules/retrieve_flu_subtype_and_segment.nf'

workflow SUBTYPE_AND_SEGMENT_FLU {
    take:
        /*
        meta must have the following key:
            - taxid_name : "NC_002020.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 8, complete sequence"
        */
        sample_meta_channel // meta

    main:

        // Where possible, retrieve the flu type and flu segment from taxid_name string in the meta
        retrieve_flu_subtype_and_segment(sample_meta_channel)

        // Add the flu type and segment to the metadata map
        subtyped_sample_meta_channel = retrieve_flu_subtype_and_segment.out
        | map{meta, type, flu_segment -> 
            meta.type = type
            meta.flu_segment = flu_segment
            meta
        }

    emit:
        subtyped_sample_meta_channel // meta
}   

workflow {
    manifest_channel = Channel.fromPath(params.manifest_file)
    | splitCsv(header: true, sep: '\t')
    | map { row -> 
    meta = [taxid_name:row.taxid_name]
    }

    SUBTYPE_AND_SEGMENT_FLU(manifest_channel)
}
