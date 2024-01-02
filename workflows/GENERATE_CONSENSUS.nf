
include {bwa_alignment_and_post_processing} from '../modules/bwa_alignment.nf'
include {run_ivar} from '../modules/run_ivar.nf'

workflow GENERATE_CONSENSUS {
    take:
        //consensus_mnf // manifest path
        sample_taxid_ch // tuple (meta, reads)
        json_ch // virus_resources dict

    main:
        /*
        mnf_ch = parse_consensus_mnf(consensus_mnf) // tuple(index, sample_id, [fastq_pair], virus_id)
        */
        // assemble 
        sample_taxid_ch // tuple(meta, fastq_pair)
            | map {meta, reads ->
                // add the vresources to meta
                if (json_ch[meta.taxid].value!= null){
                    ref_gnm_path = "${params.virus_resources_dir}${json_ch[meta.taxid].reference_genome.value[0]}"
                    meta.reference_genome_path = ref_gnm_path
                }
                else {
                    meta.reference_genome_path = null
                }
                //meta.virus_resources_dct = json_ch[meta.taxid]
               
                [meta, reads]
            }
            | branch { meta, reads ->
                missing_taxids: meta.reference_genome_path == null //virus_resources_dct.value == null
                valid_taxids: meta.reference_genome_path != null //meta.virus_resources_dct.value != null
            }
            | set {input_to_virus_rsrc_ch}

        input_to_virus_rsrc_ch.valid_taxids // tuple (meta, reads)
            | map {meta, reads ->
                [meta,reads,meta.reference_genome_path]
                }
            | set {bwa_input_ch} // tuple(meta, reads, ref_genome_path)

        // align reads to reference
        bwa_alignment_and_post_processing (bwa_input_ch)
        bams_ch = bwa_alignment_and_post_processing.out // tuple (meta, [sorted_bam, bai])


        // set ivar input channel
        // prep for join on meta.id
        bams_ch
            | map {meta, bams -> tuple(meta, bams, meta.reference_genome_path)}
            | set {ivar_in_ch}
        // generate consensus
        run_ivar(ivar_in_ch)

    emit:
        run_ivar.out // tuple (meta, [fasta_files], [quality_txt_files], variant_tsv)
}

// NOTE for now everything out of ivar is emitted, after the pipeline mature a bit
//     only whatever is used by the pipeline should be output by the process.

def parse_consensus_mnf(consensus_mnf) {
    // consensus_mnf <Channel.fromPath()>
    def mnf_ch = consensus_mnf
                        | splitCsv(header: true, sep: ',')
                        | map { row ->
                            // TODO: validateMnf(row)
                            tuple(
                                row.idx,
                                row.sample_id,
                                [row.reads_1, row.reads_2],
                                row.tax_id
                            )
                        }
    return mnf_ch // tuple(index, sample_id, [fastq_pairs], virus_id)
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