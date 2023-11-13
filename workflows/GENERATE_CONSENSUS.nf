
include {bwa_alignment_and_post_processing} from '../modules/bwa_alignment.nf'
include {run_ivar} from '../modules/run_ivar.nf'

workflow GENERATE_CONSENSUS {
    take:
        consensus_mnf // manifest path
        json_ch // virus_resources dict

    main:
        mnf_ch = parse_consensus_mnf(consensus_mnf) // tuple(index, sample_id, [fastq_pair], virus_id)

        mnf_ch
            | combine(json_ch) // tuple(index, sample_id, [fastq_pair], virus_id, meta)
            | map {it -> tuple("${it[0]}_${it[1]}_${it[3]}", it[3], it[2], it[4][it[3]])} // tuple(index_sample_id, virus_id, [fastq_pair], meta[virus_id])
            | branch {
                missing_taxids: it[3] == null
                valid_taxids: it[3] != null
            }
            | set {input_to_virus_rsrc_ch}


        input_to_virus_rsrc_ch.valid_taxids
            | map {it ->
                ref_gnm_path = "${params.virus_resources_dir}${it[3]["reference_genome"][0]}" 
                tuple(it[0], it[1],it[2],ref_gnm_path)
                }
            | set {bwa_input_ch} // tuple(index_sample_id, virus_id, [fastq_pair], ref_fasta_path)
        
        // align reads to reference
        bwa_alignment_and_post_processing (bwa_input_ch)
        bams_ch = bwa_alignment_and_post_processing.out // tuple (file_id, [sorted_bam, bai])

        // set ivar input channel
        bams_ch
        | join(bwa_input_ch) // tuple(file_id. [sorted_bam, bai], virus_id, [fastq_pairs], ref_fasta)
        | map { it -> tuple(it[0], it[1], it[4])} //tuple (file_id, [sorted_bam, bai], ref_fasta)
        | set {ivar_in_ch}

        // generate consensus
        run_ivar(ivar_in_ch)

    emit:
        run_ivar.out // tuple (file_id, [fasta_files], [quality_txt_files], variant_tsv)
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