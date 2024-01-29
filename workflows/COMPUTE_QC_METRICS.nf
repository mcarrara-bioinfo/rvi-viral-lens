include {compute_percentage_covered_bases} from '../modules/compute_percentage_covered_bases.nf'
include {compute_number_of_aligned_reads} from '../modules/compute_number_of_aligned_reads.nf'
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
        // compute percentage of bases covered
        qc_metrics_In_ch
            | map {meta, fasta_file, quality_files, variant_tsv ->
                // store fasta_files at meta
                meta.consensus_fa = fasta_file
                [meta, fasta_file, meta.ref_files[0]]
            }
            | set {percentage_covered_bases_In_ch}
        
        compute_percentage_covered_bases(percentage_covered_bases_In_ch)
        
        // compute total mapped reads
        compute_percentage_covered_bases.out
            | map {meta, stdout_str -> 
                    tokens_lst = stdout_str.tokenize("|")
                    // store results on meta
                    meta.percentage_genome_coverage = tokens_lst[0]
                    meta.consensus_lenght_sequence = tokens_lst[1]
                    meta.consensus_total_N_calls = tokens_lst[2]
                    // perpare for depth computation
                    tuple(meta, meta.bam_file)
            }
            | set {compute_n_aligned_read_In_ch}

        compute_number_of_aligned_reads(compute_n_aligned_read_In_ch)

        compute_number_of_aligned_reads.out
            | map {meta, stdout_str -> 
                stdout_lst = stdout_str.tokenize("|")
                meta.total_mapped_reads = stdout_lst[0]
                meta.total_unmapped_reads = stdout_lst[1]
                tuple(meta)
            }
            | set {qc_Out_ch}

    emit:
        qc_Out_ch // (meta)
}