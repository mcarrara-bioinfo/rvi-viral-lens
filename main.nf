#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include {check_generate_consensus_params; parse_consensus_mnf_meta} from './workflows/GENERATE_CONSENSUS.nf'
include {check_sort_reads_params} from './workflows/SORT_READS_BY_REF.nf'

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
    check_main_params()
    // ==========================
    // Map reads to virus

    // sort reads by taxon
    // input = per-sample fastq manifest; output = per-sample, per-taxon fastq manifest
    
    if (params.entry_point == "sort_reads"){

        SORT_READS_BY_REF(params.manifest)
        sample_taxid_ch = SORT_READS_BY_REF.out // tuple (meta, reads)
    }
    
    // generate consensus
    if (params.entry_point == "consensus_gen"){
        // process manifest
        sample_taxid_ch = parse_consensus_mnf_meta(params.consensus_mnf)
    }

    // load virus settings
    json_params = PipelineParameters.readParams(params.virus_resources_json)
    json_ch = Channel.value(json_params)

    GENERATE_CONSENSUS(sample_taxid_ch, json_ch)
    // TO DO: //
    // Do consensus sequence analysis
    // Do virus specific analysis
    // SARS-CoV-2
    // Flu
}

def __check_if_params_file_exist(param_name, param_value){

  def error = 0

  if (!(param_value==null)){
    param_file = file(param_value)
    if (!param_file.exists()){
      log.error("${param_file} does not exist")
      error +=1
    }
  }

  if (param_value==null){
    log.error("${param_name} must be provided")
    error +=1
  }
  // ----------------------
  return error
}

def check_main_params(){

    def errors = 0
    def valid_entry_points = ["sort_reads", "consensus_gen"]
    
    // check if execution mode is valid
    if (!valid_entry_points.contains(params.entry_point)){
        log.error("The execution mode provided (${params.entry_point}) is not valid. valid modes = ${valid_entry_points}")
        errors += 1
    }

    if (params.entry_point == "sort_reads"){
        errors += check_sort_reads_params()
    }

    if (params.entry_point=="consensus_gen"){
        // check if manifest was provided
        errors += __check_if_params_file_exist("consensus_mnf", params.consensus_mnf)
    }
    errors += __check_if_params_file_exist("virus_resources", params.virus_resources_json)

    //errors += check_generate_consensus_params()

    if (errors > 0) {
        log.error("Parameter errors were found, the pipeline will not run.")
        exit 1
    }
}
