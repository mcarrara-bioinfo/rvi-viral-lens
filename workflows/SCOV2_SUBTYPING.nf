include {run_pangolin} from '../modules/run_pangolin.nf'

workflow SCOV2_SUBTYPING {
    take:
        /*
        meta must have the following keys:
            - id
            - taxid
            - sample_id
            - ref_files
        */
        consensus_seq_ch // tuple (meta, consensus_seq)

    main:

        // generate consensus
        run_pangolin(consensus_seq_ch)
        run_pangolin.out
            .map { meta, consensus_seq, lineage -> 
                meta.virus_subtype = lineage
                [meta, consensus_seq]
            }
            .set {scov2_subtype_out_ch}
        scov2_subtype_out_ch.view()
    emit:
        scov2_subtype_out_ch // tuple (meta, consensus)
}
