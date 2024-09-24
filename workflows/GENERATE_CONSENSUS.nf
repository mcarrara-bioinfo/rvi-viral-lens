include {bwa_alignment_and_post_processing} from '../modules/bwa_alignment.nf'
include {run_ivar} from '../modules/run_ivar.nf'

workflow GENERATE_CONSENSUS {
    take:
        /*
        meta must have the following keys:
            - id
            - taxid
            - sample_id
            - ref_files
        */
        sample_taxid_ch // tuple (meta, reads)
        
    main:
        // prepare bwa input channel
        sample_taxid_ch
            | map {meta, reads -> tuple(meta,reads,meta.ref_files)}
            | set {bwa_input_ch} // tuple(meta, reads, ref_genome_paths)

        // align reads to reference
        bwa_alignment_and_post_processing (bwa_input_ch)
        bams_ch = bwa_alignment_and_post_processing.out // tuple (meta, [sorted_bam, bai])

        // set ivar input channel
        bams_ch
            | map {meta, bams -> 
                meta.bam_file = bams[0]
                tuple(meta, bams)}
            | set {ivar_in_ch}

        // generate consensus
        run_ivar(ivar_in_ch)

        // add mpileup output file to meta
        run_ivar.out // tuple (meta, fasta_file, mpileup_file)
            | map {meta, fasta_file, mpileup_file -> 
                meta.mpileup_file = mpileup_file
                tuple(meta, fasta_file)}
            | set {out_ch}

    emit:
        out_ch // tuple (meta, fasta_file)
}

def parse_consensus_mnf_meta(consensus_mnf) {
    // consensus_mnf <Channel.fromPath()>
    def mnf_ch =  Channel.fromPath(consensus_mnf)
                        | splitCsv(header: true, sep: ',')
                        | map {row -> 
                            // set meta
                            meta = [sample_id: row.sample_id,
                                    taxid: row.taxid,
                                    ref_files: row.ref_files.split(";").collect()]

                            meta.id = "${row.sample_id}.${row.taxid}"
                            
                            
                            // set files
                            reads = [row.reads_1, row.reads_2]

                            // declare channel shape
                            tuple(meta, reads)
                        }
    return mnf_ch // tuple(index, [fastq_pairs])
}

def check_generate_consensus_params(){
    // < placeholder >
    def errors = 0
    return errors
}
