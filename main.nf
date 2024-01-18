#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include {check_generate_consensus_params; parse_consensus_mnf_meta} from './workflows/GENERATE_CONSENSUS.nf'
include {check_sort_reads_params} from './workflows/SORT_READS_BY_REF.nf'

include {SORT_READS_BY_REF} from './workflows/SORT_READS_BY_REF.nf'
include {GENERATE_CONSENSUS} from './workflows/GENERATE_CONSENSUS.nf'

// Main entry-point workflow
workflow {
    // === 1 - Process input ===
    check_main_params()
    // ==========================
    // Map reads to virus

    // sort reads by taxon
    // input = per-sample fastq manifest; output = per-sample, per-taxon fastq manifest
    
    if (params.entry_point == "sort_reads"){
        // check if 
        SORT_READS_BY_REF(params.manifest)
        sample_taxid_ch = SORT_READS_BY_REF.out // tuple (meta, reads)
    }
    
    // generate consensus
    if (params.entry_point == "consensus_gen"){
        // process manifest
        sample_taxid_ch = parse_consensus_mnf_meta(params.consensus_mnf)
    }


    GENERATE_CONSENSUS(sample_taxid_ch)

    // branching output from consensus for subtyping
    GENERATE_CONSENSUS.out // [meta, [fasta_files], [quality_txt_files], variant_tsv]
      .branch { meta, fasta_files, quality_files, variant_tsv ->
        flu_subtyping_workflow_in_ch: meta.taxid_name.contains("${params.flu_keyword}")
        scv2_subtyping_workflow_in_ch: meta.taxid_name.contains("${params.scv2_keyword}")
        no_subtyping_ch: true
      }
      .set {consensus_seq_out_ch}

    //consensus_seq_out_ch.flu_subtyping_workflow_in_ch.view()
    //consensus_seq_out_ch.scv2_subtyping_workflow_in_ch.view()
    //consensus_seq_out_ch.no_subtyping_ch.view()

    // TO DO:
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

    //errors += check_generate_consensus_params()

    if (errors > 0) {
        log.error("Parameter errors were found, the pipeline will not run.")
        exit 1
    }
}
/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html
 */
workflow.onComplete {
  // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_bold = params.monochrome_logs ? '' : "\033[1m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  println """
  Pipeline execution summary
  ---------------------------
  Completed at : ${workflow.complete}
  Duration     : ${workflow.duration}
  Success      : ${c_bold}${workflow.success ? c_green : c_red}${workflow.success}${c_reset}
  Results Dir  : ${file(params.outdir)}
  Work Dir     : ${workflow.workDir}
  Exit status  : ${workflow.exitStatus}
  Error report : ${workflow.errorReport ?: '-'}
  """.stripIndent()
}
/*
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: \n ${workflow.errorMessage}"
}
*/