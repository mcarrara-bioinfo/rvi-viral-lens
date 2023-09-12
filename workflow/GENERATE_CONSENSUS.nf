
include {bwa_alignment_and_post_processing} from '../modules/bwa_alignment.nf'
include {run_ivar} from '../modules/run_ivar.nf'

workflow GENERATE_CONSENSUS {
    take:
        input_ch // tuple(file_id, virus_id, fastq_pairs, ref_fasta)
    main:
        // align reads to reference
        bwa_alignment_and_post_processing (input_ch)
        bams_ch = bwa_alignment_and_post_processing.out // tuple (file_id, [sorted_bam, bai])

        // set ivar input channel
        bams_ch
        | join(input_ch) // tuple(file_id. [sorted_bam, bai], virus_id, fastq_pairs, ref_fasta)
        | map { it -> tuple( it[0], it[1], it[4])} //tuple (file_id, [sorted_bam, bai], ref_fasta)
        | set {ivar_in_ch}

        // generate consensus
        run_ivar(ivar_in_ch)

    emit:
        run_ivar.out // tuple (file_id, [fasta_files], [quality_txt_files], variant_tsv)
}

// NOTE for now everything out of ivar is emitted, after the pipeline mature a bit
//     only whatever is used by the pipeline should be output by the process. 