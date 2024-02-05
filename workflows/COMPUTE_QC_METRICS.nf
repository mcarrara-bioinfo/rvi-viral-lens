include {compute_depth_and_coverage} from '../modules/compute_depth_and_coverage.nf'
//include {compute_valid_calls_percentage} from '../modules/compute_valid_calls_percentage.nf'

workflow COMPUTE_QC_METRICS {
    take:
        /*
        meta must have the following keys:
            - id
            - taxid
            - sample_id
            - ref_files
        */

        qc_metrics_In_ch // [meta, fasta_file, [quality_txt_files], variant_tsv]

    main:

        // store fasta_files at meta
        qc_metrics_In_ch
            | map {meta, fasta_file, quality_files, variant_tsv ->
                meta.consensus_fa = fasta_file 
                [meta, meta.ref_files[0], meta.bam_file]
            }
            | set {depth_and_coverage_In_ch}

        compute_depth_and_coverage(depth_and_coverage_In_ch)

        compute_depth_and_coverage.out
            | map {meta, read_depth, percentage_genome_coverage, genome_size -> 
                meta.read_depth = read_depth
                meta.percentage_genome_coverage = percentage_genome_coverage
                meta.genome_size = genome_size
                [meta]
            }
            | set {qc_Out_ch}
        //qc_Out_ch.view()
        //qc_Out_ch
        //    | map {meta -> tuple(meta, meta.consensus_fa[0], meta.ref_files[0][0])}
        //    | set {valid_calls_percentage_In_ch}

        //valid_calls_percentage_In_ch.view()
        //compute_valid_calls_percentage(valid_calls_percentage_In_ch)
        //compute_valid_calls_percentage.out.view()

    emit:
        qc_Out_ch // (meta)
}