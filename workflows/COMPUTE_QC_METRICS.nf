include {compute_depth_and_coverage} from '../modules/compute_depth_and_coverage.nf'
include {compute_percentage_covered_bases} from '../modules/compute_percentage_covered_bases.nf'

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

        // compute percentage of bases covered
        qc_metrics_In_ch
            | map {meta, fasta_file, quality_files, variant_tsv ->
                meta.consensus_fa = fasta_file 
                [meta, fasta_file, meta.ref_files[0]]
            }
            | set {percentage_covered_bases_In_ch}
        
        compute_percentage_covered_bases(percentage_covered_bases_In_ch)
        compute_percentage_covered_bases.out
            | map {meta, stdout_str -> 
                    tokens_lst = stdout_str.tokenize("|")
                    // store results on meta
                    meta.percentage_genome_coverage = tokens_lst[0]
                    meta.consensus_lenght_sequence = tokens_lst[1]
                    meta.consensus_total_N_calls = tokens_lst[2]
                    // perpare for depth computation
                    tuple(meta)
            }
            //| view()
        compute_depth_and_coverage(depth_and_coverage_In_ch)


        
        compute_depth_and_coverage.out
            | map {meta, read_depth, percentage_genome_coverage, genome_size -> 
                meta.read_depth = read_depth
                meta.percentage_genome_coverage = percentage_genome_coverage
                meta.genome_size = genome_size
                tuple(meta)
            }
            | set {qc_Out_ch}

    emit:
        qc_Out_ch // (meta)
}