include {run_qc_script} from '../modules/run_qc_script.nf'

workflow COMPUTE_QC_METRICS {
    take:
        /*

        Obtain QC Metrics

        The `COMPUTE_QC_METRICS` workflow is designed to compute
        quality control (QC) metrics for consensus sequences
        generated from sequencing data. The workflow processes
        each sample's data to evaluate the quality and coverage
        of the generated consensus sequences. The QC metrics 
        include the percentage of bases covered, the percentage
        of N bases, the longest segment without N bases, and 
        read alignment statistics (total reads aligned, unmapped
        and mapped).

        meta must have the following keys:
            - id
            - taxid
            - sample_id
            - bam_file
            - ref_files
            - mpileup_file

        check docs/workflow/COMPUTE_QC_METRICS.md for a more extensive documentation
        */

        qc_metrics_In_ch // [meta, fasta_file]

    main:
        qc_metrics_In_ch
            | map {meta, fasta_file ->
                // store fasta_files at meta
                meta.consensus_fa = fasta_file
                tuple(meta, meta.bam_file, fasta_file, meta.ref_files[0], meta.mpileup_file)
            }
            | set{qc_script_In_ch}

        // compute percentage of bases covered
        run_qc_script(qc_script_In_ch)
        
        // populate meta with the qc values
        run_qc_script.out
            | map {meta, qc_csv, stdout_str -> 
                tokens_lst = stdout_str.tokenize(",")
                meta.percentage_of_N_bases = tokens_lst[1]
                meta.percentage_genome_coverage = tokens_lst[2]
                meta.longest_no_N_segment = tokens_lst[3]
                meta.total_aligned_reads = tokens_lst[4]
                meta.total_unmapped_reads = tokens_lst[10].replace("\n","")
                meta.total_mapped_reads = tokens_lst[8]
                
                tuple(meta, meta.bam_file)
            }
            | set {qc_Out_ch}
   emit:
        qc_Out_ch // (meta)

}

// DEV NOTES: We should consider to make reference files on the input channels

workflow {
    manifest_channel = Channel.fromPath(params.manifest_file)
    | splitCsv(header: true, sep: ',')
    | map { row ->
        meta = [id:row.id,
            taxid:row.taxid,
            sample_id:row.sample_id,
            bam_file:row.bam_file,
            ref_files:[row.ref_files],
            mpileup_file:row.mpileup_file
        ]
        [meta, row.fasta_file]
    }
    COMPUTE_QC_METRICS(manifest_channel)
}
>>>>>>> workflows/COMPUTE_QC_METRICS.nf
