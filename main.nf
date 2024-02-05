#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include {check_generate_consensus_params} from './workflows/GENERATE_CONSENSUS.nf'
include {check_sort_reads_params} from './workflows/SORT_READS_BY_REF.nf'

include {SORT_READS_BY_REF} from './workflows/SORT_READS_BY_REF.nf'
include {GENERATE_CONSENSUS} from './workflows/GENERATE_CONSENSUS.nf'
include {SCOV2_SUBTYPING} from './workflows/SCOV2_SUBTYPING.nf'
include {COMPUTE_QC_METRICS} from './workflows/COMPUTE_QC_METRICS.nf'
include {FLU_SUBTYPING} from './workflows/FLU_SUBTYPING.nf'
include {GENERATE_CLASSIFICATION_REPORT} from './workflows/GENERATE_CLASSIFICATION_REPORT.nf'
/*
* ANSI escape codes to color output messages
*/
ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"
ANSI_BOLD = "\033[1m"
  
log.info """${ANSI_RESET}
  ===========================================
  Viral Pipeline [Dev - Prototype]
  Used parameters:
  -------------------------------------------
  --> general pipeline parameters:

    --entry_point              : ${params.entry_point}
    --containers_dir           : ${params.containers_dir}
    --results_dir              : ${params.results_dir}

  --> SORT_READS_BY_REF workflow parameters:

    --db_path                  : ${params.db_path}
    --db_library_fa_path       : ${params.db_library_fa_path}
    --min_reads_for_taxid      : ${params.min_reads_for_taxid}

  --> GENERATE_CONSENSUS workflow parameters:
    --consensus_mnf            : ${params.consensus_mnf}
    --depth_treshold           : ${params.depth_treshold}
    --mapping_quality_treshold : ${params.mapping_quality_treshold}

  --> viral subtyping branching parameters:
    --scv2_keyword             : ${params.scv2_keyword}
    --flu_keyword              : ${params.flu_keyword}

  ------------------------------------------
  Runtime data:
  -------------------------------------------
  Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
  Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
  Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
  Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
  ------------------------------------------
""".stripIndent()

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
        consensus_mnf = SORT_READS_BY_REF.out
    }
    
    // generate consensus
    if (params.entry_point == "consensus_gen"){
        // process manifest
        consensus_mnf = Channel.fromPath(params.consensus_mnf, checkIfExists: true)
    }
    // 1.0 - load virus settings
    json_params = PipelineParameters.readParams(params.virus_resources_json)
    json_ch = Channel.value(json_params)

    GENERATE_CONSENSUS(sample_taxid_ch)

    // branching output from consensus for subtyping

    COMPUTE_QC_METRICS(GENERATE_CONSENSUS.out)

    COMPUTE_QC_METRICS.out
      .branch{ it -> 
        flu_subtyping_workflow_in_ch: it[0].taxid_name.contains("${params.flu_keyword}")
        scv2_subtyping_workflow_in_ch: it[0].taxid_name.contains("${params.scv2_keyword}")
        no_subtyping_ch: true
      }
      .set {qc_metrics_out_ch}
    

    if (params.do_scov2_subtyping == true){
      qc_metrics_out_ch.scv2_subtyping_workflow_in_ch
        .map {it -> [it[0], it[0].consensus_fa]}
        .set {scov_2_subt_In_ch}
      SCOV2_SUBTYPING(scov_2_subt_In_ch)
      //report_in_ch.concat(SCOV2_SUBTYPING.out)
      SCOV2_SUBTYPING.out.set{scov2_subtyped_ch}
    }

    if (params.do_flu_subtyping == true){
      FLU_SUBTYPING(qc_metrics_out_ch.flu_subtyping_workflow_in_ch)
      FLU_SUBTYPING.out.set{flu_subtyped_ch}
    }
  
  // TODO handle with no subtypinh is requested
  if (!params.do_flu_subtyping == true){
    flu_subtyped_ch = Channel.empty()
  }

  if (!params.do_scov2_subtyping == true){
    scov2_subtyped_ch = Channel.empty()
  }
  qc_metrics_out_ch.no_subtyping_ch.concat(scov2_subtyped_ch, flu_subtyped_ch)
    .set{report_in_ch}

  report_in_ch
  GENERATE_CLASSIFICATION_REPORT(report_in_ch)
}

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
