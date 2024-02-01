include {retrieve_flu_subtype_and_segment} from '../modules/retrieve_flu_subtype_and_segment.nf'

workflow FLU_SUBTYPING {
    take:
        /*
        meta must have the following key:
            - taxid_name : "NC_002020.1 Influenza A virus (A/Puerto Rico/8/1934(H1N1)) segment 8, complete sequence"
        */
        sample_meta_channel // meta

    main:
        sample_meta_channel
            .map{it -> tuple(it[0], it[0].taxid_name)}
            .set{taxid_name_ch}
        // Where possible, retrieve the flu type and flu segment from taxid_name string in the meta
        retrieve_flu_subtype_and_segment(taxid_name_ch)

        // Add the flu type and segment to the metadata map
        subtyped_sample_meta_channel = retrieve_flu_subtype_and_segment.out
        | map{meta, type, flu_segment -> 
            meta.virus_subtype = type
            if (flu_segment == "Null"){
                flu_segment = null
            }
            meta.flu_segment = flu_segment
            tuple(meta)
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
