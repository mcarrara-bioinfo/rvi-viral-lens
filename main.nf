#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include {SORT_READS_BY_REF} from './workflows/SORT_READS_BY_REF.nf'
include {GENERATE_CONSENSUS} from './workflows/GENERATE_CONSENSUS.nf'

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
    // 1.0 - load virus settings
    json_params = PipelineParameters.readParams(params.virus_resources_json)
    json_ch = Channel.value(json_params)

    // ==========================
    // Map reads to virus

    // SORT_READS_BY_REF
    SORT_READS_BY_REF(params.manifest)
    consensus_mnf = SORT_READS_BY_REF.out

    // Gen consensus
    GENERATE_CONSENSUS(consensus_mnf, json_ch)
    // Do consensus sequence analysis

    // Do virus specific analysis

    // SARS-CoV-2

    // Flu

}
