include {run_trimmomatic} from "../modules/run_trimmomatic.nf"
include {run_fastq2fasta} from "../modules/run_fastq2fasta.nf"
include {run_trf} from "../modules/run_trf.nf"
include {run_rmRepeatFromFq} from "../modules/run_rmRepeatFromFq.nf"
include {run_sra_human_scrubber} from "../modules/run_scrubber.nf"

workflow PREPROCESSING {
    /*
    -----------------------------------------------------------------
    Preprocessing fastq files


    -----------------------------------------------------------------
    # Inputs
    - **Sample Taxid Channel **: A channel containing tuples of 

    -----------------------------------------------------------------
    # Key Processes

    -----------------------------------------------------------------
    # Outputs
    - `run_ivar.out`: A channel containing tuples of metadata and the generated consensus FASTA file.

    */

    take:

        sample_taxid_ch // tuple (meta, reads)

    main:
        reads_ch = sample_taxid_ch.map{meta, reads -> tuple (meta, reads[0], reads[1])}
        // run trimmomatic
        if (params.run_trimmomatic){
            run_trimmomatic(reads_ch)
            reads_ch = run_trimmomatic.out.paired_channel // tuple (meta, read_1, read_2)
        }

        // run trf
        if (params.run_trf){
            
            // preapare fastq channel to be join by id
            reads_ch.map{meta, fq_1, fq_2 -> 
                tuple (meta.id, meta, [fq_1, fq_2])
                }
                .set {fqs_ch}
            
            // convert
            run_fastq2fasta(reads_ch)
            run_fastq2fasta.out // tuple (meta, fasta_1, fasta_2)
                | set {trf_in_ch}
            
            run_trf(trf_in_ch)
            run_trf.out.paired_trf // tuple (meta, trf_out_1, trf_out_2)
                | map {meta, trf_out_1, trf_out_2 -> 
                    tuple (meta.id,[trf_out_1, trf_out_2])
                 }
                | set {trf_ch} // tuple (meta.id, [trfs_out])

            fqs_ch // tuple (meta.id, meta, [fqs])
                | join(trf_ch) // tuple (meta.id, meta, [fqs], [trfs_out])
                | map {id, meta, fqs, trfs -> tuple(meta, fqs[0], fqs[1], trfs[0], trfs[1])}
                | set {rmTRFfromFq_In_ch}
            run_rmRepeatFromFq(rmTRFfromFq_In_ch)
            run_rmRepeatFromFq.out.fastqs
                | set {reads_ch} // tuple (meta, trf_fq_1, trf_fq_2)
        }

        // run human-sra-scrubble
        if (params.run_scrubble){
            run_sra_human_scrubber(reads_ch)
            reads_ch = run_sra_human_scrubber.out
        }
    reads_ch
      .map{meta, reads_1, reads_2 -> tuple(meta, [reads_1, reads_2])}
      .set {out_ch}
    emit:
        out_ch // tuple (meta, fastq_pair)

}

def parse_mnf_meta(preprocessing_mnf) {
    // consensus_mnf <Channel.fromPath()>
    println(preprocessing_mnf)
    def mnf_ch =  Channel.fromPath(preprocessing_mnf)
                        | splitCsv(header: true, sep: ',')
                        | map {row -> 
                            // set meta
                            meta = [
                                id: row.sample_id,
                                sample_id: row.sample_id
                                ]

                            // set files
                            reads = [row.reads_1, row.reads_2]

                            // declare channel shape
                            tuple(meta, reads)
                        }
    return mnf_ch // tuple(index, [fastq_pairs])
}

// run preprocessing in isolation
workflow {
    mnf_ch = parse_mnf_meta(params.preprocessing_mnf)
    PREPROCESSING(mnf_ch)
}