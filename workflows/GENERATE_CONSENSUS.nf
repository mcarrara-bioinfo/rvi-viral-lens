
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
            | map {meta, reads -> [meta,reads,meta.ref_files]}
            | set {bwa_input_ch} // tuple(meta, reads, ref_genome_paths)

        // align reads to reference
        bwa_alignment_and_post_processing (bwa_input_ch)
        bams_ch = bwa_alignment_and_post_processing.out // tuple (meta, [sorted_bam, bai])

        // set ivar input channel
        bams_ch
            | map {meta, bams -> tuple(meta, bams, meta.ref_files)}
            | set {ivar_in_ch}

        // generate consensus
        run_ivar(ivar_in_ch)

        // NOTE for now everything out of ivar is emitted by the process, as the pipeline evolves
        //     we need to review to publish only whatever is used by the pipeline should be on the output channel.


    emit:
        run_ivar.out // tuple (meta, [fasta_files], [quality_txt_files], variant_tsv)
}


def parse_consensus_mnf_meta(consensus_mnf) {
    // consensus_mnf <Channel.fromPath()>
    def mnf_ch =  Channel.fromPath(consensus_mnf)
                        | splitCsv(header: true, sep: ',')
                        | map {row -> 
                            // set meta
                            meta = [sample_id: row.sample_id,
                                    taxid: row.taxid]
                            meta.id = "${meta.sample_id}.${meta.taxid}"

                            // set files
                            reads = [row.reads_1, row.reads_2]

                            // declare channel shape
                            [meta, reads]
                        }
    return mnf_ch // tuple(index, [fastq_pairs])
}

def check_generate_consensus_params(){

    def errors = 0
    // Is virus resource param set to something?
    if (param.virus_resources_json == null){
        log.error("No virus resource json file set")
        errors +=1
    }
    // if yes, is it a file which exists?
    if (params.virus_resources_json){
        virus_res_json = file(params.virus_resources_dir)
        if (!virus_res_json.exists()){
            log.error("The virus resource json provided (${params.irods_manifest}) does not exist.")
            errors += 1
        }
        //TODO
        //else {
        //    validate_virus_resource(params.irods_manifest, params.panels_settings)
        //}
    }
    return errors
}