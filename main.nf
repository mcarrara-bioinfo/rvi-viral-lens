#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include {load_json_file} from './workflow/load_json.nf'
include {bwa_alignment_and_post_processing} from './modules/bwa_alignment.nf'
class PipelineParameters {
    
    // Function that parses json output file 
    public static Map readParams(json_file) {
        def jsonSlurper = new groovy.json.JsonSlurper()
        String pipelineparameter = new File(json_file).text
        def Map configparam = (Map) jsonSlurper.parseText(pipelineparameter)
        return configparam
    }

    // Function that writes a map to a json file
    public static void writeParams(params, filename) {
        def json = new groovy.json.JsonBuilder(params)
        def myFile = new File(filename)
        myFile.write(groovy.json.JsonOutput.prettyPrint(json.toString()))
    }
}
// Main entry-point workflow
workflow {

    // === 1 - Process input ===
    // 1.1 - Get reads files 
    reads_channel_fqgz = channel.fromFilePairs("${params.reads_dir}/*_R{1,2}*.fq.gz")
    reads_channel_fagz = channel.fromFilePairs("${params.reads_dir}/*_R{1,2}*.fastq.gz")
    reads_channel_fagz
        | concat(reads_channel_fqgz) // tuple (sample_id, fq_pairs)
        | map{it -> tuple(it[0], it[0].split("_")[-1], it[1])} // tuple (sample_id, virus_code, reads)
        | set {reads_ch}

    // 1.2 - load virus settings
    //virus_resources = load_json_file(params.virus_resources_json)
    //println(virus_resources.sarscov2.reference_genome[0])

    json_params = PipelineParameters.readParams(params.virus_resources_json)
    json_ch = Channel.value(json_params)

    // add meta
    reads_ch
        | combine(json_ch) // tuple (sample_id, viral_id, read_pairs, meta)
        | set {input_ch}

    // ==========================

    // Do whole reference genome mapping

    // generate bwa in channel
    input_ch
        |map{it -> 
            // NOTE remember we need to handle multiple references, this is a temporary fix
            def reference_gnm = "${params.virus_resources_dir}${it[3][it[1]]['reference_genome'][0]}.fasta"
                tuple(
                    it[0], it[2], reference_gnm//"${params.virus_resources_dir}/${it[3][v_id]}"
                )
            }
        | set {bwa_input_ch}
    bwa_input_ch.view()
    bwa_alignment_and_post_processing (bwa_input_ch)
    bwa_alignment_and_post_processing.out.view()
    // "${reference_gnm}.amb", "${reference_gnm}.ann",
    //                               "${reference_gnm}.bwt", "${reference_gnm}.pac", "${reference_gnm}.sa"]
    //bwa-mem
    //samtools-sort
    //samtools mpileup
    
    // Do genome consensus generation
    //ivar - > consensus
    
    // Do alignment of consensus to genome to ref

    // Do species specific

    // general

    // SARS-CoV-2

    // Flu


}
