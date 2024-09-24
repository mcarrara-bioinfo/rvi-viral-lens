include {run_pangolin} from '../modules/run_pangolin.nf'

workflow SCOV2_SUBTYPING {
    take:
        /*
        Subtype SCOV2 sequences

        The SCOV2_SUBTYPING workflow is designed to 
        determine the SARS-CoV-2 lineage (subtype) of
        consensus sequences using the PANGOLIN tool.
        This workflow takes in a channel of consensus
        sequences along with their metadata, runs the
        PANGOLIN lineage classification, and outputs
        updated metadata with the assigned lineage.

        meta must have the following keys:
            - id
            - taxid
            - sample_id

        check docs/workflow/SCOV2_SUBTYPING.md for a more extensive documentation 
        */
        consensus_seq_ch // tuple (meta, consensus_seq)

    main:

        // get SCOV2 lineage classification
        run_pangolin(consensus_seq_ch)
        run_pangolin.out
            .map { meta, consensus_seq, lineage -> 
                meta.virus_subtype = lineage
                tuple(meta)
            }
            .set {scov2_subtype_out_ch}

    emit:
        scov2_subtype_out_ch // tuple (meta, consensus)
}
